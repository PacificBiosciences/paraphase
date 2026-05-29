use crate::assembly::assembly_result::AssembledPaths;
use crate::assembly::variant_graph::Graph as VariantGraph;
use crate::io::json::{ReadAlignmentId, ReadFingerprintMap};
use crate::phaser;
use crate::phaser::base_qual;
use crate::toolkit::range;
use crate::toolkit::util::{raw_qpos, DError, HashMap};

use vstr::{VStr, VString};

use itertools::{iproduct, Itertools};
use rust_htslib::bam::Read;

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

pub mod defaults {
    pub const FLANKING_SPURIOUS_BP: i32 = 10;
    pub const MIN_BASE_QUALITY: u8 = 25; // In Python: paraphase.phaser.Phaser.MEAN_BASE_QUAL
}

pub type ReadSupport<'a> = (
    BTreeMap<VStr<'a>, Vec<VStr<'a>>>,
    BTreeMap<VStr<'a>, Vec<u32>>,
    (BTreeMap<VStr<'a>, i32>, BTreeMap<VStr<'a>, f64>),
);

impl phaser::Phaser {
    /// Generate hap: read map from read->hap map.
    #[must_use]
    pub fn simplify_read_haps(
        reads: &ReadFingerprintMap,
    ) -> (
        BTreeMap<VStr<'_>, Vec<&ReadAlignmentId>>,
        BTreeMap<&ReadAlignmentId, VStr<'_>>,
    ) {
        let mut haps_to_reads = BTreeMap::<VStr<'_>, Vec<&ReadAlignmentId>>::new();
        let mut reads_to_haps = BTreeMap::<&ReadAlignmentId, VStr<'_>>::new();
        for (read, hap) in reads {
            haps_to_reads.entry(hap.vstr()).or_default().push(read);
            reads_to_haps.insert(read, hap.vstr());
        }
        (haps_to_reads, reads_to_haps)
    }

    /// Identify likely spurious haplotypes and remove them from the candidate set.
    ///
    /// Compares near-identical hap pairs and uses local pileup evidence around the
    /// mismatch site to decide whether one hap is likely an artifact.
    pub(crate) fn adjust_spurious_haplotypes<'a>(
        &self,
        current_hap_asm: &'a BTreeMap<VStr<'a>, Vec<VStr<'a>>>,
        flanking_bp: Option<i32>,
        min_base_quality: Option<u8>,
    ) -> Result<AssembledPaths, DError> {
        let min_base_quality = min_base_quality.unwrap_or(defaults::MIN_BASE_QUALITY);
        let flanking_bp = i64::from(flanking_bp.unwrap_or(defaults::FLANKING_SPURIOUS_BP));
        let mut passing = AssembledPaths::from_seqs(current_hap_asm.keys());
        let mut suspicion = Vec::new();
        for (hap1, hap2) in
            iproduct!(current_hap_asm.keys(), current_hap_asm.keys()).filter(|x| x.0 != x.1)
        {
            let (matches, mismatches, sites) = hap1
                .iter()
                .copied()
                .zip(hap2.iter().copied())
                .enumerate()
                .filter(|(_id, (x, y))| *x != b'x' && *y != b'x')
                .fold(
                    (0i32, 0i32, vec![]),
                    |(mut matches, mut mismatches, mut sites), (idx, (x, y))| {
                        let is_1_or_2 = |x: u8| -> bool { matches!(x, b'1' | b'2') };
                        if x == y {
                            matches += 1;
                        } else if is_1_or_2(x) && is_1_or_2(y) {
                            mismatches += 1;
                            sites.push(idx as i32);
                        }
                        (matches, mismatches, sites)
                    },
                );
            if matches >= 5 && mismatches == 1 && sites.len() == 1 {
                let mismatch_site = &self.het_sites[sites[0] as usize];
                let mismatch_pos = mismatch_site.pos;
                let Some(hap1_reads) = current_hap_asm.get(hap1) else {
                    continue;
                };
                let Some(hap2_reads) = current_hap_asm.get(hap2) else {
                    continue;
                };
                if let Some(pair) = if hap1_reads.len() <= 5 && hap2_reads.len() >= 6 {
                    Some((hap2, hap1))
                } else if hap2_reads.len() <= 5 && hap1_reads.len() >= 6 {
                    Some((hap1, hap2))
                } else {
                    None
                } {
                    suspicion.push((pair, mismatch_pos));
                }
            }
        }

        suspicion.sort();
        suspicion.dedup();
        for ((hap1, hap2), mismatch_pos) in &suspicion {
            let Some(hap1_reads) = current_hap_asm.get(*hap1) else {
                continue;
            };
            let Some(hap2_reads) = current_hap_asm.get(*hap2) else {
                continue;
            };
            let mut hap1_at_pos = BTreeSet::<VString>::new();
            let mut hap2_at_pos = BTreeSet::<VString>::new();
            let mut bam = self.try_realigned_bam()?;
            let tid = self
                .genome_tid()
                .map(|x| x as i32)
                .ok_or("adjust_spurious_haplotypes: missing chr tid")?;
            bam.fetch((tid, *mismatch_pos - 1, *mismatch_pos))?;
            for x in bam.pileup() {
                let x = x?;
                match i64::from(x.pos()).cmp(mismatch_pos) {
                    Ordering::Less => continue,
                    Ordering::Greater => break,
                    Ordering::Equal => {}
                };
                for aln in x.alignments().filter(|x| !x.is_del() && !x.is_refskip()) {
                    let record = aln.record();
                    if base_qual(&aln) < min_base_quality {
                        continue;
                    }
                    let qpos = raw_qpos(&aln) as i64;
                    let read_name = self.get_read_name(&record);
                    let has_hap1 = hap1_reads.contains(&VStr::from(&read_name[..]));
                    let has_hap2 = hap2_reads.contains(&VStr::from(&read_name[..]));
                    if (has_hap1 || has_hap2)
                        && qpos >= flanking_bp
                        && (qpos + flanking_bp < record.seq_len() as i64)
                    {
                        let seq = record.seq().as_bytes();
                        let start = qpos - flanking_bp;
                        let end = qpos + flanking_bp;
                        let slice = &seq[start as usize..end as usize];
                        if has_hap1 {
                            hap1_at_pos.insert(slice.into());
                        }
                        if has_hap2 {
                            hap2_at_pos.insert(slice.into());
                        }
                    }
                }
            }
            if hap1_at_pos.intersection(&hap2_at_pos).next().is_some() {
                passing.remove(&VString::from(*hap2));
            }
        }
        Ok(passing)
    }

    /// Compute per-haplotype read support counts in stable high-support windows.
    ///
    /// This mirrors Paraphase behavior by looking for the longest contiguous variant
    /// range where all haplotypes have enough support, then counting reads in a
    /// narrow center slice of that range.
    fn get_read_count_window<'a>(
        uniquely_supporting_haps: &BTreeMap<VStr<'a>, Vec<&'a VString>>,
    ) -> Option<(usize, usize)> {
        if uniquely_supporting_haps.is_empty() {
            return None;
        }
        let first_hap_reads = uniquely_supporting_haps.values().next()?;
        let first_read_hap = first_hap_reads.first()?;
        let nvar = first_read_hap.len();
        let nhap = uniquely_supporting_haps.len();
        let mut hap_base_counts = vec![Vec::<i32>::new(); nhap];
        for (hap_counts, matches) in hap_base_counts
            .iter_mut()
            .zip(uniquely_supporting_haps.values())
        {
            for i in 0..nvar {
                let counted_bases = matches
                    .iter()
                    .filter(|x| !matches!(x[i], b'x' | b'0'))
                    .count();
                hap_counts.push(counted_bases as i32);
            }
        }
        let mut ranges = Vec::new();
        assert!(hap_base_counts.iter().all(|x| x.len() == nvar));
        for fingerprint_idx in 0..nvar {
            if hap_base_counts
                .iter()
                .map(|x| x[fingerprint_idx])
                .min()
                .is_some_and(|x| x >= 5)
            {
                for j_idx in (fingerprint_idx + 1)..nvar {
                    if j_idx == nvar - 1
                        || hap_base_counts
                            .iter()
                            .map(|x| x[j_idx])
                            .min()
                            .map_or(true, |x| x < 5)
                    {
                        ranges.push(range::I64::new(fingerprint_idx as i64, j_idx as i64));
                        break;
                    }
                }
            }
        }
        let longest_range = ranges
            .iter()
            .sorted_by(|x, y| y.len().cmp(&x.len()))
            .next()?;

        let mid = (longest_range.end + longest_range.start) / 2;
        let nstart = std::cmp::max(mid - 1, longest_range.start);
        let nend = std::cmp::min(mid + 1, longest_range.end);
        if nend == nstart {
            return None;
        }
        Some((nstart as usize, nend as usize))
    }

    /// Compute per-haplotype unique/non-unique read support counts in stable high-support windows.
    ///
    /// This mirrors Paraphase behavior by looking for the longest contiguous variant
    /// range where all haplotypes have enough support, then counting reads in a
    /// narrow center slice of that range.
    fn get_read_counts<'a>(
        uniquely_supporting_haps: &BTreeMap<VStr<'a>, Vec<&'a VString>>,
        raw_read_haps: &HashMap<&str, &VString>,
        nonuniquely_supporting_reads: &HashMap<&str, Vec<VStr<'a>>>,
    ) -> (BTreeMap<VStr<'a>, i32>, BTreeMap<VStr<'a>, f64>) {
        let mut read_count_unique = BTreeMap::<VStr<'_>, i32>::new();
        let mut read_count_nonunique = BTreeMap::<VStr<'_>, f64>::new();
        let Some((nstart, nend)) = Self::get_read_count_window(uniquely_supporting_haps) else {
            return (read_count_unique, read_count_nonunique);
        };

        for (hap, matches) in uniquely_supporting_haps {
            let l_reads = matches
                .iter()
                .filter(|x| x[nstart..nend].iter().any(|b| *b != b'x'))
                .count();
            read_count_unique.insert(*hap, l_reads as i32);
            read_count_nonunique.insert(*hap, l_reads as f64);
        }

        // evenly distribute nonunique reads across all haplotypes supported.
        for (read_name, haps_supported) in nonuniquely_supporting_reads
            .iter()
            .filter(|(_read, haps)| haps.len() > 1)
        {
            let Some(read_hap) = raw_read_haps.get(read_name) else {
                continue;
            };
            if !read_hap[nstart..nend].iter().any(|b| *b != b'x') {
                continue;
            }
            let read_weight = 1.0 / haps_supported.len() as f64;
            for hap in haps_supported {
                *read_count_nonunique.entry(*hap).or_default() += read_weight;
            }
        }
        (read_count_unique, read_count_nonunique)
    }

    /// Build unique/non-unique read support maps for assembled haplotypes.
    ///
    /// Returns:
    /// - uniquely supporting reads per assembled haplotype,
    /// - non-uniquely supporting read-to-haplotype-id links,
    /// - read-count summary from [`Self::get_read_counts`].
    pub(crate) fn get_read_support<'a>(
        reads: &'a ReadFingerprintMap,
        haps_to_reads: &'a BTreeMap<VStr<'a>, Vec<&ReadAlignmentId>>,
        assembled_haps: &'a [VString],
    ) -> Result<ReadSupport<'a>, DError> {
        let support = VariantGraph::match_reads_and_haps(reads, assembled_haps, None)?;
        let uniquely_supporting_haps = &support.support_reads;
        let mut uniquely_supporting_reads = assembled_haps
            .iter()
            .map(|hap| (VStr::from(hap), vec![]))
            .collect::<BTreeMap<_, _>>();
        let mut uniquely_supporting_haps_by_seq = BTreeMap::<VStr<'a>, Vec<&'a VString>>::new();
        let mut nonuniquely_supporting_reads_by_name = HashMap::<&str, Vec<VStr<'a>>>::default();
        let raw_read_haps = reads
            .iter()
            .map(|(read_id, hap)| (read_id.read_name.as_str(), hap))
            .collect::<HashMap<&str, &VString>>();
        for (haplotype_index, reads) in uniquely_supporting_haps {
            let ass_hap = VStr::from(&assembled_haps[*haplotype_index as usize]);
            uniquely_supporting_haps_by_seq.insert(ass_hap, reads.iter().map(|x| x.path).collect());
        }
        for (read, hap_ids) in support
            .supporting_haps_per_read
            .iter()
            .filter(|(_read, haps)| haps.len() > 1)
        {
            let haps_supported = hap_ids
                .iter()
                .map(|hap_id| VStr::from(&assembled_haps[*hap_id as usize]))
                .collect::<Vec<_>>();
            nonuniquely_supporting_reads_by_name.insert(read.name, haps_supported);
        }

        let read_counts = Self::get_read_counts(
            &uniquely_supporting_haps_by_seq,
            &raw_read_haps,
            &nonuniquely_supporting_reads_by_name,
        );
        for (haplotype_index, reads) in uniquely_supporting_haps.iter() {
            let ass_hap = VStr::from(&assembled_haps[*haplotype_index as usize]);
            let read_support = uniquely_supporting_reads.entry(ass_hap).or_default();
            for read in reads.iter() {
                let read_hap = read.path.vstr();
                for name in &haps_to_reads[&read_hap] {
                    read_support.push(VStr::from(&name.read_name[..]));
                }
            }
            read_support.sort();
            read_support.dedup();
        }

        let mut nonuniquely_supporting_reads = BTreeMap::<VStr<'a>, Vec<u32>>::new();
        for (read, hap_ids) in support
            .supporting_haps_per_read
            .iter()
            .filter(|(_read, haps)| haps.len() > 1)
        {
            nonuniquely_supporting_reads.insert(VStr::from(&read.name), hap_ids.clone());
        }
        Ok((
            uniquely_supporting_reads,
            nonuniquely_supporting_reads,
            read_counts,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::json::ReadAlignmentId;

    #[test]
    fn simplify_read_haps_matches_python_case() {
        let reads = BTreeMap::from([
            (ReadAlignmentId::from_name("r1"), VString::from("12")),
            (ReadAlignmentId::from_name("r2"), VString::from("21")),
        ]);

        let (haps_to_reads, reads_to_haps) = phaser::Phaser::simplify_read_haps(&reads);

        let haps_to_reads = haps_to_reads
            .into_iter()
            .map(|(hap, reads)| {
                (
                    hap.to_string(),
                    reads
                        .into_iter()
                        .map(std::string::ToString::to_string)
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<BTreeMap<_, _>>();
        let reads_to_haps = reads_to_haps
            .into_iter()
            .map(|(read, hap)| (read.to_string(), hap.to_string()))
            .collect::<BTreeMap<_, _>>();

        assert_eq!(
            haps_to_reads,
            BTreeMap::from([
                (String::from("12"), vec![String::from("r1")]),
                (String::from("21"), vec![String::from("r2")]),
            ])
        );
        assert_eq!(
            reads_to_haps,
            BTreeMap::from([
                (String::from("r1"), String::from("12")),
                (String::from("r2"), String::from("21")),
            ])
        );
    }
}
