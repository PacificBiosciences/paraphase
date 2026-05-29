use crate::phaser;
use crate::phaser::base_qual;
use crate::phaser::PhasedResult;
use crate::toolkit::range;
use crate::toolkit::site_selection::{self, CandidateSite};
use crate::toolkit::util::{raw_qpos, DError};

use vstr::{VStr, VString};

use itertools::{enumerate, Itertools};
use rust_htslib::bam::{self, pileup::Indel, Read};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

/// Find first and last non-gap (`x`) character indices.
#[must_use]
pub fn get_start_end(x: &[u8]) -> (i32, i32) {
    let start = x
        .iter()
        .enumerate()
        .find(|x| *x.1 != b'x')
        .map_or(x.len() - 1, |x| x.0) as i32;
    let end = x
        .iter()
        .rev()
        .enumerate()
        .find(|x| *x.1 != b'x')
        .map_or(0, |(count, _item)| x.len() - 1 - count) as i32;
    (start, end)
}

/// Assignment for token - reference, alt, or other (missing or special)>
#[derive(Debug, PartialEq, Copy, Clone, Eq, PartialOrd, Ord, Hash)]
pub enum Assignment {
    Ref,
    Alt,
    Dot,
}

#[derive(Debug, PartialEq, Copy, Clone, Eq)]
pub enum Genotype {
    One,
    Zero,
    Dot,
}

/// Haplotype summary information.
/// Set of variants, reference coordinate range, and whether or not it is truncated.
/// `boundary_gene2` marks gene2 locations.
#[derive(Clone, Debug)]
pub struct HapInfo {
    pub variants: Vec<CandidateSite>,
    pub boundary: range::I64,
    pub boundary_gene2: Option<range::I64>,
    pub is_truncated: Vec<String>,
}

#[derive(Clone, Debug, serde::Deserialize, serde::Serialize)]
pub struct HapInfoForJson {
    pub variants: Vec<String>,
    pub boundary: [i64; 2],
    pub boundary_gene2: Option<[Option<i64>; 2]>,
    pub is_truncated: Vec<String>,
}

impl std::convert::From<&HapInfo> for HapInfoForJson {
    fn from(info: &HapInfo) -> Self {
        let variants = info
            .variants
            .iter()
            .map(std::string::ToString::to_string)
            .collect::<Vec<_>>();
        let bound_to_output = |x: &range::I64| -> [i64; 2] { [x.start + 1, x.end + 1] };
        let boundary = bound_to_output(&info.boundary);
        let boundary_gene2 = info.boundary_gene2.as_ref().map(|x| {
            let to_nullable = |v: i64| if v < 0 { None } else { Some(v + 1) };
            [to_nullable(x.start), to_nullable(x.end)]
        });
        let is_truncated = info.is_truncated.clone();
        Self {
            variants,
            boundary,
            boundary_gene2,
            is_truncated,
        }
    }
}

impl phaser::Phaser {
    /// Infer a haplotype genotype at one site from supporting read assignments.
    pub(crate) fn get_genotype_in_hap(
        &self,
        var_reads: &BTreeMap<String, Assignment>,
        hap_reads: &[String],
        hap_reads_nonunique: &BTreeMap<&String, &Vec<VString>>,
    ) -> Result<Genotype, DError> {
        let mut hap_reads_contain_var = hap_reads
            .iter()
            .filter(|x| var_reads.contains_key(&x[..]))
            .collect::<Vec<_>>();
        if hap_reads_contain_var.len() < 3 {
            for name in hap_reads_nonunique
                .keys()
                .filter(|x| var_reads.contains_key(&x[..]))
            {
                hap_reads_contain_var.push(name);
            }
        }
        let hap_reads_contain_var = hap_reads_contain_var
            .iter()
            .filter_map(|x| var_reads.get(&x[..]).copied())
            .collect::<Vec<Assignment>>();

        if hap_reads_contain_var.len() >= 3 {
            let ctr = hap_reads_contain_var
                .into_iter()
                .collect::<counter::Counter<_, i32>>();
            let most_common = ctr.k_most_common_ordered(2);
            let top = most_common[0];
            let top_threshold = top.1 as f32 * 0.15;
            let threshold = if top_threshold > 2. {
                2.
            } else {
                top_threshold
            };
            if most_common.len() == 1 || most_common[1].1 as f32 <= threshold {
                match top.0 {
                    Assignment::Alt => return Ok(Genotype::One),
                    Assignment::Ref => return Ok(Genotype::Zero),
                    _ => {}
                }
            }
        }
        Ok(Genotype::Dot)
    }

    /// Classify reads at a candidate variant as ref/alt/other for haplotype checks.
    pub fn check_variants_in_haplotypes(
        &self,
        site: &CandidateSite,
        ref_seq: &[u8],
        min_bq: Option<u8>,
    ) -> Result<BTreeMap<String, Assignment>, DError> {
        let mut dreads = BTreeMap::<String, Assignment>::new();
        let var_size = site.var_seq.len() as i64 - site.ref_seq.len() as i64;
        let indel_base_in_read = match var_size.cmp(&0) {
            Ordering::Less => format!(
                "{}{}{}",
                std::str::from_utf8(&[site.ref_seq[0]])?,
                var_size,
                VStr::from(&site.ref_seq[1..])
            ),
            Ordering::Greater => format!(
                "{}+{}{}",
                std::str::from_utf8(&[site.ref_seq[0]])?,
                var_size,
                VStr::from(&site.var_seq[1..])
            ),
            _ => String::default(),
        };
        let ref_base = site.ref_seq[0];
        let alt_base = if var_size == 0 {
            VStr::from(&site.var_seq[..1])
        } else {
            VStr::from(&indel_base_in_read[..])
        };
        let mut reader = self.try_realigned_bam()?;
        let tid = self
            .genome_tid()
            .map(|x| x as i32)
            .ok_or("get_indel_haplotypes: missing chr tid")?;
        reader.fetch((tid, site.pos - 1, site.pos + 1))?;
        for pile in reader.pileup() {
            let pile = pile?;
            match pile.pos().cmp(&(site.pos as u32)) {
                Ordering::Less => continue,
                Ordering::Equal => {}
                Ordering::Greater => break,
            }
            let min_base_quality = min_bq.unwrap_or(
                self.settings
                    .site_selection_settings
                    .min_candidate_base_quality,
            );
            let offset = self.offset() as usize;
            for aln in pile.alignments() {
                let query_pos_raw = raw_qpos(&aln);
                let record = aln.record();
                let is_reverse = record.is_reverse();
                let refskip_char = if is_reverse { b'<' } else { b'>' };
                let bq = base_qual(&aln);
                if bq < min_base_quality {
                    continue;
                }
                let record = aln.record();
                let seq = record.seq().as_bytes();
                let mut query_seq = VString::default();
                let base = if !aln.is_del() {
                    let base = seq.get(query_pos_raw).copied().unwrap_or(b'N');
                    site_selection::maybe_strand_mark_char(base, is_reverse, false)
                } else if aln.is_refskip() {
                    refskip_char
                } else {
                    b'*'
                };
                query_seq.push(base);
                let pos = pile.pos() as usize;
                match aln.indel() {
                    Indel::Ins(x) => {
                        query_seq.push(b'+');
                        query_seq.extend_from_slice(x.to_string().as_bytes());
                        for j in 1..=(x as usize) {
                            query_seq.push(site_selection::maybe_strand_mark_char(
                                seq[j + query_pos_raw],
                                is_reverse,
                                false,
                            ));
                        }
                    }
                    Indel::Del(x) => {
                        query_seq.push(b'-');
                        query_seq.extend_from_slice(x.to_string().as_bytes());
                        for j in 1..=(x as usize) {
                            query_seq.push(site_selection::maybe_strand_mark_char(
                                ref_seq[j + pos - offset],
                                is_reverse,
                                false,
                            ));
                        }
                    }
                    bam::pileup::Indel::None => {}
                }
                let asn = match (query_seq == ref_base, query_seq == alt_base) {
                    (true, _) => Assignment::Ref,
                    (_, true) => Assignment::Alt,
                    _ => Assignment::Dot,
                };
                let this_read_names = self.get_read_names(&record, None);
                let this_read_name = this_read_names
                    .first()
                    .ok_or_else(|| {
                        phaser::Exception::new(format!(
                            "No read names extracted for record '{}'",
                            String::from_utf8_lossy(record.qname())
                        ))
                    })?
                    .to_string();
                dreads.insert(this_read_name, asn);
            }
        }

        Ok(dreads)
    }

    /// Compute genomic boundary span for a (possibly partial) haplotype string.
    ///
    /// Leading/trailing `x` segments are expanded to neighboring site intervals.
    fn get_hap_variant_ranges(&self, hap: VStr<'_>) -> range::I64 {
        let (start, end) = get_start_end(&hap);
        let range = range::I64::new(start.into(), end.into());
        let nstart_previous_pos: i64 = if range.start == 0 {
            self.left_boundary_0based()
        } else {
            self.het_sites[(range.start - 1) as usize].pos + 1
        };
        let nend_next_pos = if end == (hap.len() - 1) as i32 {
            self.right_boundary_0based()
        } else {
            self.het_sites[(end + 1) as usize].pos - 1
        };
        range::I64::new(nstart_previous_pos, nend_next_pos)
    }

    /// Find corresponding coordinates in the secondary region.
    pub fn get_range_in_other_gene(&self, pos: i64, search_range: Option<i64>) -> Option<i64> {
        let search_range = search_range.unwrap_or(200);
        self.matches.get(&pos).copied().or_else(|| {
            (pos..pos + search_range)
                .filter_map(|x| self.matches.get(&x))
                .copied()
                .next()
        })
    }

    /// Infer the 5' clip breakpoint represented by a haplotype fingerprint.
    pub fn get_5pclip_from_hap(&self, hap: &VStr) -> Result<Option<i64>, DError> {
        let het_sites = &self.het_sites;
        let hap_len = hap.iter().len();
        assert_eq!(het_sites.len(), hap_len);
        let mut clips_not_present = Vec::new();
        for (index, base) in enumerate(hap.iter()) {
            if index < hap_len - 1 {
                let next_base = hap.get(index + 1).ok_or_else(|| {
                    phaser::Exception::new(format!(
                        "Haplotype indexing failed at index {} while scanning 5p clips (hap_len={})",
                        index + 1,
                        hap_len
                    ))
                })?;
                let site_before = het_sites
                    .get(index)
                    .ok_or_else(|| {
                        phaser::Exception::new(format!(
                            "Missing het_site at index {} while scanning 5p clips",
                            index
                        ))
                    })?
                    .pos;
                let site_after = het_sites
                    .get(index + 1)
                    .ok_or_else(|| {
                        phaser::Exception::new(format!(
                            "Missing het_site at index {} while scanning 5p clips",
                            index + 1
                        ))
                    })?
                    .pos;
                for clip_position in &self.clip_5p_positions {
                    if *clip_position > site_before && *clip_position < site_after {
                        if *base == b'0' && *next_base != b'0' && *next_base != b'x' {
                            return Ok(Some(*clip_position));
                        }
                        if *next_base != b'0'
                            && *next_base != b'x'
                            && ((*base != b'0' && *base != b'x')
                                || (*base == b'x' && site_after - *clip_position < 5000))
                        {
                            clips_not_present.push(*clip_position);
                        }
                    }
                }
            }
        }
        if clips_not_present == self.clip_5p_positions {
            return Ok(Some(0));
        }
        Ok(None)
    }

    /// Given a haplotype, get its 3p clip position
    /// Infer the 3' clip breakpoint represented by a haplotype fingerprint.
    pub fn get_3pclip_from_hap(&self, hap: &VStr) -> Result<Option<i64>, DError> {
        let het_sites = &self.het_sites;
        let hap_len = hap.iter().len();
        assert_eq!(het_sites.len(), hap_len);
        let mut clips_not_present = Vec::new();
        for (index, base) in enumerate(hap.iter()) {
            if index < hap_len - 1 {
                let next_base = hap.get(index + 1).ok_or_else(|| {
                    phaser::Exception::new(format!(
                        "Haplotype indexing failed at index {} while scanning 3p clips (hap_len={})",
                        index + 1,
                        hap_len
                    ))
                })?;
                let site_before = het_sites
                    .get(index)
                    .ok_or_else(|| {
                        phaser::Exception::new(format!(
                            "Missing het_site at index {} while scanning 3p clips",
                            index
                        ))
                    })?
                    .pos;
                let site_after = het_sites
                    .get(index + 1)
                    .ok_or_else(|| {
                        phaser::Exception::new(format!(
                            "Missing het_site at index {} while scanning 3p clips",
                            index + 1
                        ))
                    })?
                    .pos;
                for clip_position in &self.clip_3p_positions {
                    if *clip_position > site_before && *clip_position < site_after {
                        if *next_base == b'0' && *base != b'0' && *base != b'x' {
                            return Ok(Some(*clip_position));
                        }
                        if *base != b'0'
                            && *base != b'x'
                            && ((*next_base != b'0' && *next_base != b'x')
                                || (*next_base == b'x' && *clip_position - site_before < 5000))
                        {
                            clips_not_present.push(*clip_position);
                        }
                    }
                }
            }
        }
        if clips_not_present == self.clip_3p_positions {
            return Ok(Some(0));
        }
        Ok(None)
    }

    /// Summarize variants per hap.
    /// Build per-haplotype variant summaries and boundary metadata.
    ///
    /// Returns map keyed by haplotype label for downstream JSON/VCF reporting.
    pub fn output_variants_in_haps(
        &mut self,
        result: &PhasedResult,
        known_del: &BTreeMap<char, String>,
        assembled_haps: BTreeMap<VStr<'_>, String>,
    ) -> Result<BTreeMap<String, HapInfo>, DError> {
        use std::str::FromStr;

        if result.assemblies.main_haps.is_empty() {
            return Ok(BTreeMap::new());
        }
        let het_sites = self.het_sites.clone();
        let no_phasing_sites = self.het_sites_no_phasing.clone();
        let mut hap_info = BTreeMap::<String, HapInfo>::new();
        let mut hap_variants = assembled_haps
            .values()
            .cloned()
            .map(|x| (x, BTreeSet::new()))
            .collect::<BTreeMap<String, BTreeSet<String>>>();

        let ref_seq = {
            let faidx = self.make_faidx()?;
            let (chrom, start, stop) = self.parsed_nchr_0based().ok_or_else(|| {
                phaser::Exception::new(format!(
                    "Malformed realign region '{}'; expected 'chr:start-end' format",
                    self.realign_region
                ))
            })?;

            faidx
                .fetch_seq(chrom, start as usize, stop as usize)?
                .to_owned()
        };
        if !result.uniquely_supporting_reads.is_empty() {
            for var in &no_phasing_sites {
                let var_reads = self.check_variants_in_haplotypes(var, &ref_seq, None)?;
                let mut haps_with_variant = Vec::new();
                for (hap, hap_name) in &assembled_haps {
                    let hap_vstring: VString = hap.into();
                    let Some(hap_reads) = result.uniquely_supporting_reads.get(&hap_vstring) else {
                        continue;
                    };
                    let hap_reads = &hap_reads[..];
                    let hap_reads_nonunique = result
                        .nonuniquely_supporting_reads
                        .iter()
                        .filter(|(_read, hap_set)| hap_set.contains(&hap.into()))
                        .collect::<BTreeMap<_, _>>();
                    let genotype =
                        self.get_genotype_in_hap(&var_reads, hap_reads, &hap_reads_nonunique)?;
                    if genotype == Genotype::One {
                        haps_with_variant.push(hap_name);
                    }
                }
                if haps_with_variant.is_empty() {
                    self.het_sites_no_phasing.retain(|x| x != var);
                } else {
                    for hap_name in haps_with_variant {
                        if let Some(v) = hap_variants.get_mut(hap_name) {
                            v.insert(var.to_string());
                        }
                    }
                }
            }
        }

        for (hap, hap_name) in &assembled_haps {
            log::debug!("Finding variants in haplotype {hap_name}");
            let mut hap_bounds = self.get_hap_variant_ranges(*hap);
            log::debug!("Haplotype bounds: {hap_bounds}");
            let mut is_truncated = Vec::new();
            for (base, het_site) in hap.iter().zip(het_sites.iter()) {
                if *base == b'2' {
                    hap_variants
                        .entry(hap_name.clone())
                        .or_default()
                        .insert(het_site.to_string());
                } else if let Some(del_name) = known_del.get(&(*base as char)) {
                    hap_variants.entry(hap_name.clone()).or_default();
                    if !hap_variants[hap_name].contains(del_name) {
                        hap_variants
                            .entry(hap_name.clone())
                            .or_default()
                            .insert(del_name.clone());
                    }
                }
            }
            let mut filtered_hom = self.hom_sites.clone();
            for del_data in self.del_data.iter() {
                let del_name = &del_data.name();
                if hap_variants[hap_name].contains(del_name) {
                    log::debug!("This haplotype has deletion {del_name}.");
                    log::debug!(
                        "Removing homozygous sites within the deletion range, between {} and {}",
                        del_data.threep().start,
                        del_data.fivep().end
                    );
                    filtered_hom.retain(|x| {
                        x.pos < del_data.threep().start || x.pos > del_data.fivep().end
                    });
                    log::debug!("Filtered homozygous sites after checking against long deletion {del_name}: {filtered_hom:?}");
                }
            }
            let clip_position_5p = self.get_5pclip_from_hap(hap)?;
            if let Some(clip_position_5p_value) = clip_position_5p {
                if clip_position_5p_value != 0 {
                    filtered_hom.retain(|x| x.pos > clip_position_5p_value);
                    hap_bounds.start = std::cmp::max(hap_bounds.start, clip_position_5p_value);
                    hap_variants
                        .entry(hap_name.clone())
                        .or_default()
                        .insert(format!("{}_clip_5p", clip_position_5p_value + 1));
                    if clip_position_5p_value + 1 > self.gene_start() {
                        is_truncated.push(String::from("5p"));
                    }
                }
            }
            let clip_position_3p = self.get_3pclip_from_hap(hap)?;
            if let Some(clip_position_3p_value) = clip_position_3p {
                if clip_position_3p_value != 0 {
                    filtered_hom.retain(|x| x.pos < clip_position_3p_value);
                    hap_bounds.end = std::cmp::min(hap_bounds.end, clip_position_3p_value);
                    hap_variants
                        .entry(hap_name.clone())
                        .or_default()
                        .insert(format!("{}_clip_3p", clip_position_3p_value + 1));
                    if clip_position_3p_value + 1 < self.gene_end() {
                        is_truncated.push(String::from("3p"));
                    }
                }
            }
            for site in &filtered_hom {
                hap_variants
                    .entry(hap_name.clone())
                    .or_default()
                    .insert(site.to_string());
            }
            hap_bounds.start = std::cmp::max(hap_bounds.start, self.left_boundary_0based());
            hap_bounds.end = std::cmp::min(hap_bounds.end, self.right_boundary_0based());
            let boundary_gene2 = if self
                .locus_config()
                .gene2_region(self.settings.genome == "37")
                .is_some()
            {
                let (start, end) = [hap_bounds.start, hap_bounds.end]
                    .into_iter()
                    .map(|x| self.get_range_in_other_gene(x, None))
                    .next_tuple()
                    .ok_or_else(|| {
                        phaser::Exception::new(format!(
                            "Failed to compute boundary_gene2 coordinate pair for hap '{}' (start={}, end={})",
                            hap_name, hap_bounds.start, hap_bounds.end
                        ))
                    })?;
                match (start, end) {
                    (Some(start), Some(end)) => Some(range::I64::new(
                        std::cmp::min(start, end),
                        std::cmp::max(start, end),
                    )),
                    (None, Some(end)) => Some(range::I64::new(-2, end)),
                    (Some(start), None) => Some(range::I64::new(start, -2)),
                    (None, None) => Some(range::I64::new(-2, -2)),
                }
            } else {
                None
            };
            let var_tmp = hap_variants.get(hap_name).ok_or_else(|| {
                phaser::Exception::new(format!(
                    "Haplotype '{}' missing from hap_variants map",
                    hap_name
                ))
            })?;
            let getpos = |x: &str| {
                x.split_terminator('_')
                    .next()
                    .and_then(|x| x.parse::<i64>().ok())
            };
            let var_tmp = var_tmp
                .iter()
                .filter_map(|x| {
                    let pos = getpos(x)?;
                    if pos > hap_bounds.start && pos <= hap_bounds.end + 1 {
                        Some(x)
                    } else {
                        None
                    }
                })
                .sorted_by_cached_key(|x| getpos(&x[..]))
                .collect::<Vec<_>>();
            let mut variants = Vec::with_capacity(var_tmp.len());
            for var in var_tmp
                .into_iter()
                .map(|s| CandidateSite::from_str(s.as_str()))
            {
                variants.push(var?);
            }
            let info = HapInfo {
                variants,
                boundary: hap_bounds,
                boundary_gene2,
                is_truncated,
            };
            hap_info.insert(hap_name.into(), info);
        }

        Ok(hap_info)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config;
    use crate::phaser;
    use crate::phaser::Phaser;
    use crate::toolkit::site_selection::CandidateSite;
    use crate::toolkit::util;
    use std::str::FromStr;

    fn build_smn1_phaser() -> Option<Phaser> {
        let outdir = tempfile::TempDir::new().expect("tempdir should build");
        let Some(genome_path) = std::env::var("HG38")
            .ok()
            .map(|x| x.trim_end_matches(".mmi").to_string())
        else {
            log::warn!("Skipping hap_variants test because the HG38 environment variable is not configured.");
            return None;
        };
        let settings = phaser::Settings::new(
            "HG00733",
            (genome_path, util::test_file("bams/HG00733.smn1.bam")),
            outdir.path(),
            "smn1",
            &config::Region::try_load(None).expect("region config should load"),
            None,
            None,
            String::from("38"),
            None,
            0.03,
            false,
        );
        let gene_config = config::Gene::try_load(None).expect("gene config should load");
        Some(Phaser::new(settings, Some(gene_config), None, None).expect("phaser should build"))
    }

    #[test]
    fn get_start_end_matches_python_case() {
        let (start, end) = get_start_end(b"xx1211xxx");
        assert_eq!((start, end), (2, 5));
    }

    #[test]
    fn get_hap_variant_ranges_matches_python_partial_case() {
        let Some(mut phaser) = build_smn1_phaser() else {
            return;
        };
        phaser.het_sites = vec![
            CandidateSite::from_str("70917101_A_C").unwrap(),
            CandidateSite::from_str("70917111_A_C").unwrap(),
            CandidateSite::from_str("70917150_A_C").unwrap(),
            CandidateSite::from_str("70917200_A_C").unwrap(),
            CandidateSite::from_str("70917300_A_C").unwrap(),
            CandidateSite::from_str("70917400_A_C").unwrap(),
        ];

        let bounds = phaser.get_hap_variant_ranges(VStr::from("x121xx"));
        assert_eq!(bounds.start, 70917101);
        assert_eq!(bounds.end, 70917298);
    }
}
