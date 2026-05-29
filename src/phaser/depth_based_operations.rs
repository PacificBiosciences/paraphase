use crate::depth::{self};
use crate::phaser::{Exception, Phaser};
use crate::phaser::{HapInfo, PhasedResult};
use crate::toolkit::math::depth_prob;
use crate::toolkit::util;
use crate::toolkit::util::DError;

use rust_htslib::{bam, bam::Read};
use vstr::{VStr, VString};

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

pub(crate) fn counts(
    x: &bam::pileup::Pileup,
    ref_seq: &[u8],
    min_base_quality: u8,
    offset: usize,
) -> counter::Counter<VString, i32> {
    use crate::phaser::base_qual;
    use crate::toolkit::site_selection::maybe_strand_mark_char;

    let mut ret = counter::Counter::new();
    let mark_strand = false;
    for aln in x.alignments() {
        let record = aln.record();
        let query_pos = util::raw_qpos(&aln);
        let seq = record.seq().as_bytes();
        if base_qual(&aln) < min_base_quality {
            continue;
        }
        let is_reverse = record.is_reverse();
        let refskip_char = if is_reverse { b'<' } else { b'>' };
        let base = if !aln.is_del() && !aln.is_refskip() {
            seq.get(query_pos).copied().unwrap_or(b'N')
        } else if aln.is_refskip() {
            refskip_char
        } else {
            b'*'
        };
        let pos = record.pos() as usize;
        let mut query_seq = VString::from(vec![base]);
        match aln.indel() {
            bam::pileup::Indel::Ins(x) => {
                query_seq.push(b'+');
                query_seq.extend_from_slice(x.to_string().as_bytes());
                for j in 1..=(x as usize) {
                    query_seq.push(maybe_strand_mark_char(
                        seq[j + query_pos],
                        is_reverse,
                        mark_strand,
                    ));
                }
            }
            bam::pileup::Indel::Del(x) => {
                query_seq.push(b'-');
                query_seq.extend_from_slice(x.to_string().as_bytes());
                for j in 1..=(x as usize) {
                    query_seq.push(maybe_strand_mark_char(
                        ref_seq[j + pos - offset],
                        is_reverse,
                        mark_strand,
                    ));
                }
            }
            bam::pileup::Indel::None => {}
        }
        *ret.entry(query_seq).or_default() += 1;
    }
    ret
}

impl Phaser {
    /// Identify which haplotypes may have two copies based on depth.
    /// For each haplotype, identify the variants where it's different
    /// from other haplotypes. Check depth at those variant sites and
    /// see if the depth suggests twice coverage.
    pub(crate) fn compare_depth(
        &self,
        haps: &BTreeMap<String, HapInfo>,
        assembled_haps: &BTreeMap<VStr<'_>, String>,
        loose: bool,
        stringent: bool,
    ) -> Result<Vec<String>, DError> {
        if haps.len() <= 1 {
            return Ok(vec![]);
        }
        let mut two_cp_haps = Vec::<String>::new();

        let mut hap_name_to_seq = BTreeMap::new();
        for (hap_seq, hap_name) in assembled_haps {
            hap_name_to_seq.insert(hap_name, hap_seq);
        }
        let mut bound_start = Vec::new();
        let mut bound_end = Vec::new();
        for (hap, hap_info) in haps {
            let bound = &hap_info.boundary;
            let hap_seq = hap_name_to_seq.get(hap).ok_or_else(|| {
                Exception::new(format!(
                    "Haplotype '{}' missing from assembled haplotype name->sequence map",
                    hap
                ))
            })?;
            let hap_clip_5p = self.get_5pclip_from_hap(hap_seq)?;
            let hap_clip_3p = self.get_3pclip_from_hap(hap_seq)?;
            let mut true_start = bound.start;
            if let Some(hap_clip_5p_value) = hap_clip_5p {
                if hap_clip_5p_value != 0 {
                    true_start = hap_clip_5p_value;
                }
            }
            bound_start.push(true_start);
            let mut true_end = bound.end;
            if let Some(hap_clip_3p_value) = hap_clip_3p {
                if hap_clip_3p_value != 0 {
                    true_end = hap_clip_3p_value;
                }
            }
            bound_end.push(true_end);
        }
        let start = bound_start
            .iter()
            .max()
            .ok_or_else(|| Exception::new("compare_depth: no hap start bounds were computed"))?;
        let end = bound_end
            .iter()
            .min()
            .ok_or_else(|| Exception::new("compare_depth: no hap end bounds were computed"))?;
        let range = crate::toolkit::range::I64::new(*start, *end);
        let vars = haps
            .values()
            .flat_map(|info| {
                info.variants.iter().filter(|x| {
                    range.contains(&x.pos) && self.het_sites.contains(x) && x.is_unit_length()
                })
            })
            .collect::<BTreeSet<_>>();
        log::trace!("[compare_depth] candidate comparison variants: vars={vars:?}");
        let faidx = self.make_faidx()?;
        let (chrom, start, stop) = self
            .parsed_nchr_0based()
            .ok_or_else(|| Exception::new("compare_depth: malformatted region string"))?;
        let ref_seq = faidx.fetch_seq(chrom, start as usize, stop as usize)?;
        let offset = self.try_offset()? as usize;
        let threshold: f64 = match (stringent, loose) {
            (true, _) => 0.8,
            (_, true) => 0.5,
            _ => 0.6,
        };
        let mut reader = self.try_realigned_bam()?;
        for (hap_name, info) in haps {
            let mut sites = BTreeMap::new();
            let other_haps = haps.iter().filter(|x| x.0 != hap_name).collect::<Vec<_>>();
            let other_cn = other_haps.len();
            let this_hap_var = &info.variants[..];
            let other_haps_var = other_haps
                .iter()
                .flat_map(|(_name, info)| info.variants.iter())
                .collect::<Vec<_>>();
            for var in &vars {
                let in_this = this_hap_var.contains(var);
                let in_other = other_haps_var.contains(var);
                if in_this && !in_other {
                    sites.insert(var.pos, &var.var_seq);
                } else if !in_this
                    && other_haps_var.iter().filter(|x| *x == var).count() == other_cn
                {
                    sites.insert(var.pos, &var.ref_seq);
                }
            }
            log::trace!(
                "[compare_depth] per-haplotype discriminator sites: hap_name={hap_name}, sites={sites:?}"
            );
            let mut double_prob_count = 0usize;
            let tid = self
                .genome_tid()
                .map(|x| x as i32)
                .ok_or_else(|| Exception::new("compare_depth: missing chr tid in BAM header"))?;
            for (site, base) in &sites {
                reader.fetch((tid, *site, site + 1)).map_err(|e| {
                    Exception::new(format!(
                        "compare_depth: failed to fetch {tid}:{site}-{}: {e}",
                        site + 1
                    ))
                })?;
                for pileup in reader.pileup() {
                    let pileup = match pileup {
                        Ok(pileup) => pileup,
                        Err(e) => {
                            log::warn!(
                                "compare_depth: pileup iteration failed at site {site}: {e}"
                            );
                            continue;
                        }
                    };
                    let pos = i64::from(pileup.pos());
                    match pos.cmp(site) {
                        Ordering::Less => continue,
                        Ordering::Greater => break,
                        Ordering::Equal => {}
                    }
                    if self.settings.min_base_quality != 25 {
                        log::warn!(
                            "compare_depth: expected min_base_quality=25 for calibrated depth heuristic; observed {}",
                            self.settings.min_base_quality
                        );
                    }
                    let got_counts =
                        counts(&pileup, &ref_seq, self.settings.min_base_quality, offset);
                    log::trace!(
                        "Observed pileup base counts at site {site}: got_counts={got_counts:?}"
                    );
                    debug_assert!(
                        got_counts.keys().all(|x| **x == x.to_ascii_uppercase()),
                        "counts had non uppercase keys: {got_counts:?}"
                    );
                    let mut base_count = 0;
                    let mut total = 0;
                    for (each_base, count) in &got_counts {
                        total += count;
                        let Some(each_base_first) = each_base.first().copied() else {
                            continue;
                        };
                        if **base == each_base_first {
                            base_count += count;
                        }
                    }
                    let other_count = total - base_count;
                    let double_prob =
                        depth_prob(base_count, f64::from(other_count) / other_cn as f64);
                    log::trace!(
                        "[compare_depth] depth-probability input/output: hap_name={hap_name}, site={site}, total={total}, base_count={base_count}, double_prob={double_prob:?}"
                    );
                    if let Some(double_prob_value) = double_prob {
                        if double_prob_value[0] < 0.25 {
                            double_prob_count += 1;
                        }
                    }
                }
            }
            log::trace!(
                "[compare_depth] per-haplotype double-copy evidence summary: hap_name={hap_name}, nsite={}, double_prob_count={}",
                sites.len(),
                double_prob_count
            );
            if double_prob_count > 0
                && sites.len() >= 5
                && double_prob_count as f64 >= sites.len() as f64 * threshold
            {
                two_cp_haps.push(hap_name.clone());
            }
        }

        if two_cp_haps.len() > 1 {
            two_cp_haps.clear();
        }

        log::debug!(
            "Differentiating variant depth-based two-copy hap candidates: haps={two_cp_haps:?}"
        );
        Ok(two_cp_haps)
    }

    /// Check if the haplotype with the highest depth has twice the reads
    /// of the haplotype with the second highest depth
    /// Equivalent to get_cn2_haplotype in python
    pub fn compare_depth_by_read_count(
        &mut self,
        assembled_haps: &BTreeMap<VStr<'_>, String>,
        phase_results: &PhasedResult,
        prob_threshold: f32,
        excluded_haplotypes: &[VString],
    ) -> Vec<String> {
        let min_read_count = if self.settings.targeted { 15 } else { 10 };
        let prob_threshold_nonunique = prob_threshold.max(0.5);
        let mut two_cp = Vec::new();
        let read_counts_unique = &phase_results.read_counts.0;
        let read_counts_nonunique = &phase_results.read_counts.1;
        if assembled_haps.len() < 2 || read_counts_unique.len() < 2 {
            return two_cp;
        }
        let excluded_haplotypes = excluded_haplotypes.iter().collect::<BTreeSet<_>>();

        let haps = read_counts_unique
            .keys()
            .filter(|hap| !excluded_haplotypes.contains(hap))
            .collect::<Vec<_>>();
        if haps.len() < 2 {
            return two_cp;
        }
        let counts_unique = haps
            .iter()
            .map(|hap| *read_counts_unique.get(*hap).unwrap_or(&0))
            .collect::<Vec<_>>();
        let Some(max_count) = counts_unique.iter().copied().max() else {
            return two_cp;
        };
        let Some(cp2_hap_candidate_index) = counts_unique.iter().position(|x| *x == max_count)
        else {
            return two_cp;
        };
        let mut sorted_counts_unique = counts_unique.clone();
        sorted_counts_unique.sort_by(|a, b| b.cmp(a));
        let Some(others_max) = sorted_counts_unique.get(1).copied() else {
            return two_cp;
        };
        let Some(others_max_index) = counts_unique.iter().position(|x| *x == others_max) else {
            return two_cp;
        };
        let probs_unique = depth_prob(max_count, f64::from(others_max));

        let counts_nonunique = haps
            .iter()
            .map(|hap| *read_counts_nonunique.get(*hap).unwrap_or(&0.0))
            .collect::<Vec<_>>();
        let cp2_hap_nonunique = counts_nonunique[cp2_hap_candidate_index];
        let other_nonunique_counts = counts_nonunique
            .iter()
            .enumerate()
            .filter(|(idx, _)| *idx != cp2_hap_candidate_index)
            .map(|(_, count)| *count)
            .collect::<Vec<_>>();
        // Take the mean of the nonunique counts of all remaining haplotypes
        let other_nonunique_mean =
            other_nonunique_counts.iter().sum::<f64>() / other_nonunique_counts.len() as f64;
        let probs_nonunique = depth_prob(
            cp2_hap_nonunique.round() as i32,
            other_nonunique_mean.round(),
        );

        log::debug!(
            "Read-count two-copy check: hap_count={}, candidate_idx={}, second_idx={}, unique_counts={counts_unique:?}, unique_probs={probs_unique:?}, nonunique_counts={counts_nonunique:?}, nonunique_other_mean={other_nonunique_mean}, nonunique_probs={probs_nonunique:?}",
            read_counts_unique.len(),
            cp2_hap_candidate_index,
            others_max_index
        );

        if let Some(probs_value_unique) = probs_unique {
            if probs_value_unique[0] < prob_threshold && others_max >= min_read_count {
                if let Some(probs_value_nonunique) = probs_nonunique {
                    // check that when considering nonunique reads, the haplotype still has more reads than expected
                    // use a very lenient cutoff for now
                    if !(probs_value_nonunique[0] < prob_threshold_nonunique
                        && other_nonunique_mean >= f64::from(min_read_count))
                    {
                        return two_cp;
                    }
                } else {
                    return two_cp;
                }

                let cp2_hap_candidate = haps[cp2_hap_candidate_index];
                if let Some(hap_name) = assembled_haps.get(&cp2_hap_candidate.vstr()) {
                    two_cp.push(hap_name.clone());
                } else {
                    log::warn!(
                        "Skipping two-copy assignment because the top read-count haplotype was not found in the assembled haplotype map."
                    );
                }
            }
        }
        two_cp
    }

    #[allow(clippy::too_many_arguments)]
    pub fn adjust_depth(
        &mut self,
        assembled_haps: BTreeMap<VStr<'_>, String>,
        haps: BTreeMap<String, HapInfo>,
        phase_results: PhasedResult,
        expect_cn4: bool,
        stringent: bool,
        loose: bool,
        prob_threshold: f32,
    ) -> Result<(Vec<String>, usize), DError> {
        let mut two_cp_haps = Vec::<String>::new();
        if assembled_haps.len() == 1 && self.init_het_sites.is_empty() {
            two_cp_haps = assembled_haps.values().cloned().collect::<Vec<_>>();
        } else if self.settings.targeted {
            two_cp_haps = self.compare_depth_by_read_count(
                &assembled_haps,
                &phase_results,
                prob_threshold,
                &[],
            );
        } else if (assembled_haps.len() == 3 && !self.expect_cn2() && self.gene_name() != "BPY2")
            || (self.gene_name() == "BPY2" && assembled_haps.len() < 3)
        {
            two_cp_haps = self.compare_depth(&haps, &assembled_haps, loose, stringent)?;
            if two_cp_haps.is_empty() && !phase_results.read_counts.0.is_empty() {
                two_cp_haps = self.compare_depth_by_read_count(
                    &assembled_haps,
                    &phase_results,
                    prob_threshold,
                    &[],
                );
            }
        }

        if self.settings.targeted && two_cp_haps.is_empty() {
            let gene1_cn2 = self.locus_config().get("gene1_cn2");
            if gene1_cn2.is_some() {
                let region_length = self.right_boundary() - self.left_boundary();
                let snp_count = self.het_sites.len();
                log::debug!(
                    "Targeted CN heuristic inputs: region_length={region_length}, snp_count={snp_count}"
                );
                if snp_count as f64 > region_length as f64 * 0.008 {
                    let mut gene1_haps = Vec::new();
                    let mut gene2_haps = Vec::new();
                    for (hap_seq, hap_name) in &assembled_haps {
                        let count2 = hap_seq.iter().filter(|x| **x == b'2').count();
                        if count2 as f64 > hap_seq.len() as f64 * 0.7 {
                            gene2_haps.push(hap_name);
                        } else {
                            gene1_haps.push(hap_name);
                        }
                    }
                    log::debug!(
                        "Targeted CN heuristic hap buckets: gene1_haps={gene1_haps:?}, gene2_haps={gene2_haps:?}"
                    );
                    if gene1_haps.len() == 1 {
                        if let Some(hap) = gene1_haps.first() {
                            two_cp_haps.push((*hap).to_string());
                        }
                    }
                }
            }
        }

        let mut total_cn = assembled_haps.len() + two_cp_haps.len();
        if assembled_haps.is_empty() && self.init_het_sites.is_empty() {
            total_cn = 2;
        }
        if total_cn == 2 && !self.expect_cn2() && self.gene_name() != "BPY2" {
            if expect_cn4 {
                two_cp_haps = assembled_haps.values().cloned().collect::<Vec<_>>();
                total_cn = 4;
            } else if let Some(depth) = self.settings.depth.as_ref() {
                if self.gene_config().uses_depth(self.gene_name())
                    && depth.status() == depth::MedianDepthStatus::Passing
                {
                    let copy_number_probs =
                        depth_prob(self.region_avg_depth[0].0 as i32, depth.median as f32);
                    if let Some(copy_number_probs_value) = copy_number_probs {
                        let same_cn_prob = copy_number_probs_value[0];
                        if same_cn_prob < 0.75 {
                            total_cn = 4;
                            if two_cp_haps.is_empty() && !assembled_haps.is_empty() {
                                two_cp_haps = assembled_haps.values().cloned().collect::<Vec<_>>();
                            }
                        }
                    }
                }
            }
        }
        Ok((two_cp_haps, total_cn))
    }
}
