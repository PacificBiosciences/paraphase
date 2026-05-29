// SMN1 specific caller
use crate::io::json::GeneCall;
use crate::phaser::Phaser;
use crate::phaser::{HapInfo, HapInfoForJson};
use crate::toolkit::math::depth_prob;
use crate::toolkit::range;
use crate::toolkit::site_selection::CandidateSite;
use crate::toolkit::util::DError;
use itertools::intersperse;
use rust_htslib::bam::Read;
use serde_json::Value;
use std::cmp;
use std::collections::{BTreeMap, HashSet};
use std::ops::Index;
use vstr::VString;

#[derive(Clone, Debug)]
pub struct KnownHapMatch {
    group: String,
    nmismatch: Vec<usize>,
}

impl Phaser {
    /// check splice site
    fn check_smn1_smn2_presence(&mut self) -> Result<(usize, usize), DError> {
        let mut smn1_read_splice = HashSet::new();
        let mut smn2_read_splice = HashSet::new();
        let mut bam = self.try_realigned_bam()?;
        let tid = self.genome_tid().map(|x| x as i32).ok_or_else(|| {
            crate::phaser::Exception::new(format!(
                "Missing chromosome tid while checking SMN1/SMN2 splice support for gene '{}'",
                self.gene_name()
            ))
        })?;
        let pivot_site = self.pivot_site_0based().ok_or_else(|| {
            crate::phaser::Exception::new(format!(
                "Missing pivot site for gene '{}' while checking SMN1/SMN2 presence",
                self.gene_name()
            ))
        })?;
        bam.fetch((tid, pivot_site, pivot_site + 1))?;
        for p in bam.pileup() {
            let p = p?;
            if i64::from(p.pos()) == pivot_site {
                for pileup_read in p.alignments() {
                    if !pileup_read.is_del() && !pileup_read.is_refskip() {
                        let record = pileup_read.record();
                        let Some(qpos) = pileup_read.qpos() else {
                            continue;
                        };
                        let this_base = record.seq().index(qpos).to_ascii_uppercase();
                        let base_q = record.qual().get(qpos).copied().unwrap_or(0);
                        let qname = std::str::from_utf8(record.qname())?.to_string();
                        if base_q >= 13 {
                            if this_base == b'C' {
                                smn1_read_splice.insert(qname.clone());
                            } else if this_base == b'T' {
                                // For reverse alignments, "next" genomic base is previous query base.
                                let next_qpos = if record.is_reverse() {
                                    qpos.checked_sub(1)
                                } else if qpos + 1 < record.seq_len() {
                                    Some(qpos + 1)
                                } else {
                                    None
                                };
                                if let Some(next_qpos) = next_qpos {
                                    let next_base =
                                        record.seq().index(next_qpos).to_ascii_uppercase();
                                    if next_base != b'C' {
                                        smn2_read_splice.insert(qname.clone());
                                    }
                                }
                            }
                        }
                        log::trace!(
                            "{qname} {base_q} this_base {this_base} {}",
                            std::str::from_utf8(&[this_base])?
                        );
                    }
                }
            }
        }

        Ok((smn1_read_splice.len(), smn2_read_splice.len()))
    }

    /// Adjust SMN1 CN in various scenarios
    fn update_smn1_cn(
        &mut self,
        smn1_cn: Option<i32>,
        smn2_cn: Option<i32>,
        smn1_nread: usize,
        smn2_nread: usize,
        has_smn2: bool,
        smn1_haps: &BTreeMap<String, String>,
        two_cp_haps_candidate: &Vec<String>,
        hcn_higher: bool,
    ) -> Result<(Option<i32>, Vec<String>), DError> {
        let mut new_smn1_cn = smn1_cn;
        let mut two_cp_haps: Vec<String> = Vec::new();
        let Some(smn1_cn_value) = smn1_cn else {
            return Ok((None, two_cp_haps));
        };
        if let Some(depth) = self.settings.depth.as_ref() {
            let genome_depth = depth.median;
            if genome_depth < 20.0 && smn1_cn_value == 1 {
                return Ok((None, two_cp_haps));
            }
            let copy_number_probs = depth_prob(smn1_nread as i32, genome_depth / 2.0_f64);
            if let Some(copy_number_probs_value) = copy_number_probs {
                if smn1_cn_value == 1 {
                    // when there is only one haplotype, check if depth is
                    // consistent with haploid depth
                    let copy_one_prob = copy_number_probs_value[0];
                    if copy_one_prob < 0.05 {
                        two_cp_haps = smn1_haps.clone().into_values().collect::<Vec<String>>();
                        return Ok((Some(2), two_cp_haps));
                    }
                    if copy_one_prob < 0.25 {
                        new_smn1_cn = None;
                    }
                } else if smn1_cn_value == 2 {
                    //  scenario where two two-copy alleles are identical
                    let copy_four_prob = copy_number_probs_value[3];
                    if copy_four_prob > 0.75 {
                        two_cp_haps = smn1_haps.clone().into_values().collect::<Vec<String>>();
                        return Ok((Some(4), two_cp_haps));
                    }
                }
            }
            if (smn1_cn_value == 2 || smn1_cn_value == 3) && !self.settings.targeted {
                // check if one smn1 haplotype has more reads than others
                if two_cp_haps_candidate.len() == 1 {
                    if let Some(two_cp_hap) = two_cp_haps_candidate.first() {
                        let two_cp_hap = two_cp_hap.to_string();
                        if two_cp_hap.contains("smn1hap") {
                            two_cp_haps.push(two_cp_hap);
                            return Ok((Some(smn1_cn_value + 1), two_cp_haps));
                        }
                    }
                }
            }
        }
        if smn1_cn_value == 1 {
            let Some(smn2_cn_value) = smn2_cn else {
                return Ok((new_smn1_cn, two_cp_haps));
            };
            if smn2_cn_value >= 1 {
                let haploid_depth = (smn2_nread as f32) / (smn2_cn_value as f32);
                if haploid_depth >= 10.0 {
                    let copy_number_probs = depth_prob(smn1_nread as i32, haploid_depth);
                    // when there is only one haplotype, check if depth is
                    // consistent with haploid depth
                    if let Some(copy_number_probs_value) = copy_number_probs {
                        let copy_one_prob = copy_number_probs_value[0];
                        if copy_one_prob < 0.25 {
                            new_smn1_cn = None;
                            if smn2_cn_value > 1 && copy_one_prob < 0.05 {
                                two_cp_haps =
                                    smn1_haps.clone().into_values().collect::<Vec<String>>();
                                return Ok((Some(2), two_cp_haps));
                            }
                        }
                    }
                }
            }
        }
        // if we see more haplotypes in other regions of the gene
        if smn1_cn_value == 1 && hcn_higher {
            if !has_smn2 {
                two_cp_haps = smn1_haps.clone().into_values().collect::<Vec<String>>();
                return Ok((Some(2), two_cp_haps));
            } else {
                new_smn1_cn = None;
            }
        }
        Ok((new_smn1_cn, two_cp_haps))
    }

    /// Adjust SMN2 CN in various scenarios
    fn update_smn2_cn(
        &mut self,
        smn1_cn: Option<i32>,
        smn2_cn: Option<i32>,
        smn1_nread: usize,
        smn2_nread: usize,
        exon78_del_nreads: usize,
        smn2_haps: &BTreeMap<String, String>,
        two_cp_haps_candidate: &Vec<String>,
    ) -> Result<(Option<i32>, Vec<String>), DError> {
        let mut new_smn2_cn = smn2_cn;
        let mut two_cp_haps: Vec<String> = Vec::new();
        let Some(smn2_cn_value) = smn2_cn else {
            return Ok((None, two_cp_haps));
        };
        if let Some(depth) = self.settings.depth.as_ref() {
            let genome_depth = depth.median;
            if genome_depth < 20.0 && smn2_cn_value == 1 {
                return Ok((None, two_cp_haps));
            }
            // if smn1 cn is 3, then no need to adjust smn2 cn
            if smn2_cn_value == 1 && smn1_cn.map_or(true, |x| x <= 2) && exon78_del_nreads <= 1 {
                let copy_number_probs = depth_prob(smn2_nread as i32, genome_depth / 2.0_f64);
                if let Some(copy_number_probs_value) = copy_number_probs {
                    //when there is only one haplotype, check if depth is
                    //consistent with haploid depth
                    let copy_one_prob = copy_number_probs_value[0];
                    if copy_one_prob < 0.25 {
                        new_smn2_cn = None;
                        if copy_one_prob < 0.05 {
                            two_cp_haps = smn2_haps.clone().into_values().collect::<Vec<String>>();
                            return Ok((Some(2), two_cp_haps));
                        }
                    }
                }
            } else if smn2_cn_value == 2 && smn1_cn == Some(0) {
                let copy_number_probs = depth_prob(smn2_nread as i32, genome_depth / 2.0_f64);
                if let Some(copy_number_probs_value) = copy_number_probs {
                    //  smn1 = 0, smn2 = 4, but only see two haplotypes
                    let copy_four_prob = copy_number_probs_value[3];
                    if copy_four_prob > 0.75 {
                        two_cp_haps = smn2_haps.clone().into_values().collect::<Vec<String>>();
                        return Ok((Some(4), two_cp_haps));
                    }
                }
            }
        }
        if let Some(smn1_cn_value) = smn1_cn {
            if [0, 1].contains(&smn1_cn_value)
                && [2, 3].contains(&smn2_cn_value)
                && !self.settings.targeted
            {
                if two_cp_haps_candidate.len() == 1 {
                    if let Some(two_cp_hap) = two_cp_haps_candidate.first() {
                        let two_cp_hap = two_cp_hap.to_string();
                        if two_cp_hap.contains("smn2hap") {
                            two_cp_haps.push(two_cp_hap);
                            return Ok((Some(smn2_cn_value + 1), two_cp_haps));
                        }
                    }
                }
            }
            if smn2_cn_value == 1 && smn1_cn_value == 2 && exon78_del_nreads <= 1 {
                let haploid_depth = (smn1_nread as f32) / (smn1_cn_value as f32);
                if haploid_depth >= 10.0 {
                    let copy_number_probs = depth_prob(smn2_nread as i32, haploid_depth);
                    // when there is only one haplotype, check if depth is
                    // consistent with haploid depth
                    if let Some(copy_number_probs_value) = copy_number_probs {
                        let copy_one_prob = copy_number_probs_value[0];
                        if copy_one_prob < 0.25 {
                            new_smn2_cn = None;
                            if copy_one_prob < 0.05 {
                                two_cp_haps =
                                    smn2_haps.clone().into_values().collect::<Vec<String>>();
                                return Ok((Some(2), two_cp_haps));
                            }
                        }
                    }
                }
            }
        }
        Ok((new_smn2_cn, two_cp_haps))
    }

    fn assign_haps_to_gene(
        &self,
        ass_haps: &[String],
        has_smn1: bool,
        smn2_reads_splice: usize,
    ) -> (Vec<String>, Vec<String>, Vec<String>) {
        let found_splice = self.get_pivot_index();
        let mut smn1_haps = Vec::new();
        let mut smn2_haps = Vec::new();
        let mut smn2_del_haps = Vec::new();

        for hap in ass_haps {
            if hap.contains('3') {
                smn2_del_haps.push(hap.to_string());
            }
        }

        if let Some(splice_index) = found_splice {
            let splice_index = splice_index as usize;
            for hap in ass_haps {
                let hap_bytes = hap.as_bytes();
                if hap_bytes.get(splice_index) == Some(&b'1') {
                    smn1_haps.push(hap.to_string());
                } else if !hap.contains('3') {
                    smn2_haps.push(hap.to_string());
                }
            }
        } else if !has_smn1 || smn2_reads_splice < 2 {
            for hap in ass_haps {
                if !hap.contains('3') {
                    if has_smn1 {
                        smn1_haps.push(hap.to_string());
                    } else {
                        smn2_haps.push(hap.to_string());
                    }
                }
            }
        }
        (smn1_haps, smn2_haps, smn2_del_haps)
    }

    /// Assign a phased haplotype to the closest known SMN1/SMN2 haplogroup label.
    ///
    /// Uses two matching passes (full and `strip_c` region) to mirror Python logic.
    pub fn assign_hap_to_group(
        &mut self,
        hap_info: &HapInfo,
        known_haps: &Vec<(HapInfo, String)>,
    ) -> Result<Option<String>, DError> {
        let mut best_match = None;
        let best_match1 = self.get_best_match(hap_info, known_haps, false, 15)?;
        let best_match2 = self.get_best_match(hap_info, known_haps, true, 10)?;
        if best_match1.is_some() {
            best_match = best_match1;
        } else if best_match2.is_some() {
            best_match = best_match2;
        }
        if let Some(best_match_value) = best_match.clone() {
            if best_match_value == "S1-9" {
                let deletion = self
                    .locus_config()
                    .get("deletion1_name")
                    .and_then(|x| x.as_str())
                    .ok_or_else(|| {
                        crate::phaser::Exception::new(format!(
                            "Missing deletion1_name in region config for gene '{}'",
                            self.gene_name()
                        ))
                    })?
                    .to_string();
                let fields = deletion.split_terminator('_').collect::<Vec<_>>();
                let var_pos = fields[0].parse::<i64>()? - 1;
                let ref_base = fields[1].to_string();
                let alt_base = fields[2].to_string();
                let deletion_var = CandidateSite::new(var_pos, ref_base, alt_base);
                if hap_info.variants.contains(&deletion_var) {
                    best_match = Some(String::from("S1-9d"));
                }
            }
        }
        Ok(best_match)
    }

    /// Return the best haplogroup match under overlap and mismatch thresholds.
    ///
    /// When `strip_c` is true, matching is restricted to the configured C-rich
    /// subregion and trailing `c` suffixes are normalized.
    pub fn get_best_match(
        &mut self,
        hap_info: &HapInfo,
        known_haps: &Vec<(HapInfo, String)>,
        strip_c: bool,
        max_mismatch: usize,
    ) -> Result<Option<String>, DError> {
        let strip_c_region_start = self
            .locus_config()
            .get("strip_c_region_start")
            .and_then(|x| x.as_i64())
            .ok_or_else(|| {
                crate::phaser::Exception::new(format!(
                    "Missing strip_c_region_start in region config for gene '{}'",
                    self.gene_name()
                ))
            })?;
        let strip_c_region_end = self
            .locus_config()
            .get("strip_c_region_end")
            .and_then(|x| x.as_i64())
            .ok_or_else(|| {
                crate::phaser::Exception::new(format!(
                    "Missing strip_c_region_end in region config for gene '{}'",
                    self.gene_name()
                ))
            })?;
        let min_overlap_len: i64 = 12000;
        let vars_all = &hap_info.variants;
        let boundary = &hap_info.boundary;
        let nstart = boundary.start;
        let nend = boundary.end;
        let mut haplogroup_candidates: BTreeMap<String, Vec<usize>> = BTreeMap::new();
        for (known_hap_info, known_haplogroup) in known_haps {
            let known_hap_vars_all = &known_hap_info.variants;
            let known_bound = &known_hap_info.boundary;
            let known_hap_nstart = known_bound.start;
            let known_hap_nend = known_bound.end;
            let mut new_nstart = cmp::max(known_hap_nstart, nstart);
            let mut new_nend = cmp::min(known_hap_nend, nend);
            if strip_c {
                new_nstart = cmp::max(new_nstart, strip_c_region_start);
                new_nend = cmp::min(new_nend, strip_c_region_end);
            }
            let mut new_known_haplogroup = known_haplogroup.to_string();
            if new_nend - new_nstart >= min_overlap_len {
                if strip_c {
                    let tmp = known_haplogroup.strip_suffix("c");
                    if let Some(tmp) = tmp {
                        new_known_haplogroup = tmp.to_string();
                    }
                }
                let var1 = vars_all
                    .iter()
                    .filter(|x| x.pos >= new_nstart && x.pos <= new_nend)
                    .collect::<Vec<_>>();
                let var2 = known_hap_vars_all
                    .iter()
                    .filter(|x| x.pos >= new_nstart && x.pos <= new_nend)
                    .collect::<Vec<_>>();
                //log::debug!("var1 {:?}", var1);
                //log::debug!("group {:?} var2 {:?}", new_known_haplogroup.clone(), var2);
                let nmismatch = var1.iter().filter(|x| !var2.contains(*x)).count()
                    + var2.iter().filter(|x| !var1.contains(*x)).count();
                //log::debug!("nmismatch {:?}", nmismatch);
                haplogroup_candidates
                    .entry(new_known_haplogroup)
                    .or_default()
                    .push(nmismatch);
            }
        }
        if haplogroup_candidates.is_empty() {
            return Ok(None);
        }
        let mut sort_candidates = Vec::new();
        for (a, b) in haplogroup_candidates {
            let mut c = b.clone();
            c.sort();
            sort_candidates.push(KnownHapMatch {
                group: a,
                nmismatch: c,
            });
        }
        sort_candidates.sort_by(|a, b| a.nmismatch.iter().min().cmp(&(b.nmismatch.iter().min())));
        log::debug!(
            "SMN1 haplogroup candidates (strip_c={strip_c}, max_mismatch={max_mismatch}): {sort_candidates:?}"
        );
        let first_match = &sort_candidates[0].group;
        let first_match_nmismatch = &sort_candidates[0].nmismatch;
        let second_match_nmismatch = &sort_candidates[1].nmismatch;
        let mut best_match = None;
        if first_match_nmismatch.len() >= 2 {
            if (first_match_nmismatch[0] <= 2 || first_match_nmismatch[1] <= max_mismatch)
                && first_match_nmismatch[0] < second_match_nmismatch[0] - 2
                && first_match_nmismatch[1] < second_match_nmismatch[0]
            {
                best_match = Some(first_match.to_string());
            }
        }
        Ok(best_match)
    }

    /// Run SMN1/SMN2-specific workflow, including splice-site-based gene assignment,
    /// deletion handling, CN estimation, and haplogroup annotation.
    pub fn run_smn1(&mut self) -> Result<GeneCall, DError> {
        // Initial setup:
        // S1: Get local region.
        let seq = self.realign()?;
        // Check coverage after aligning to local reference.
        let coverage_passes = self.coverage_passes();
        let mut call = self.get_default_call();
        if !coverage_passes {
            log::debug!(
                "Coverage check failed after read alignment for gene {}; returning default call with failed_for_coverage=true",
                self.gene_name()
            );
            call.failed_for_coverage = true;
            return Ok(call);
        }

        const KNOWNHAPS_38: &[u8] = std::include_bytes!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/data/38/smn1/known_haplotypes.json"
        ));
        const KNOWNHAPS_19: &[u8] = std::include_bytes!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/data/19/smn1/known_haplotypes.json"
        ));
        const KNOWNHAPS_13: &[u8] = std::include_bytes!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/data/chm13/smn1/known_haplotypes.json"
        ));
        let tmp = if self.settings.genome == "19" || self.settings.genome == "37" {
            std::str::from_utf8(KNOWNHAPS_19)?
        } else if self.settings.genome == "chm13" {
            std::str::from_utf8(KNOWNHAPS_13)?
        } else {
            std::str::from_utf8(KNOWNHAPS_38)?
        };
        let known_haplotypes: Value = serde_json::from_str(tmp)?;
        let mut known_haps = BTreeMap::new();
        for gene_name in vec!["smn1", "smn2"] {
            let mut tmp = Vec::new();
            let hap_obj = known_haplotypes[gene_name].as_object().ok_or_else(|| {
                crate::phaser::Exception::new(format!(
                    "known_haplotypes['{}'] should be an object",
                    gene_name
                ))
            })?;
            for (_a, b) in hap_obj {
                let c = b.as_array().ok_or_else(|| {
                    crate::phaser::Exception::new("known hap entry should be an array")
                })?;
                let nlen = c.len();
                if nlen < 2 {
                    return Err("known hap entry must have at least bounds+haplogroup".into());
                }
                let known_hap_vars_all = c[..(nlen - 2)]
                    .iter()
                    .map(|x| {
                        x.as_str().ok_or_else(|| {
                            crate::phaser::Exception::new("known hap variant should be string")
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()?
                    .into_iter()
                    .map(|x| x.to_string())
                    .map(|x| {
                        let fields = x.split_terminator('_').collect::<Vec<_>>();
                        if fields.len() < 3 {
                            return Err("known hap variant should be pos_ref_alt".into());
                        }
                        let var_pos = fields[0].parse::<i64>()? - 1;
                        let ref_base = fields[1].to_string();
                        let alt_base = fields[2].to_string();
                        Ok(CandidateSite::new(var_pos, ref_base, alt_base))
                    })
                    .collect::<Result<Vec<CandidateSite>, DError>>()?;
                let known_haplogroup = c
                    .last()
                    .and_then(|x| x.as_str())
                    .ok_or_else(|| {
                        crate::phaser::Exception::new("known haplogroup should be string")
                    })?
                    .to_string();
                let known_bound = c[nlen - 2]
                    .as_array()
                    .ok_or_else(|| {
                        crate::phaser::Exception::new("known hap bounds should be array")
                    })?
                    .iter()
                    .map(|x| {
                        x.as_number().and_then(|n| n.as_i64()).ok_or_else(|| {
                            crate::phaser::Exception::new("known hap bound should be integer")
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                if known_bound.len() < 2 {
                    return Err("known hap bounds should include start and end".into());
                }
                let known_hap_nstart = known_bound[0];
                let known_hap_nend = known_bound[1];
                let known_hap_info = HapInfo {
                    variants: known_hap_vars_all,
                    boundary: range::I64::new(known_hap_nstart, known_hap_nend),
                    boundary_gene2: None,
                    is_truncated: vec![],
                };
                tmp.push((known_hap_info, known_haplogroup));
            }
            known_haps.insert(String::from(gene_name), tmp);
        }
        log::debug!(
            "Loaded known SMN haplotype references: smn1={}, smn2={}",
            known_haps.get("smn1").map_or(0, std::vec::Vec::len),
            known_haps.get("smn2").map_or(0, std::vec::Vec::len)
        );

        let (hom_sites_to_add, add_sites) = self.get_sites(&seq, None, None)?;
        let (smn1_nread_splice, smn2_nread_splice) = self.check_smn1_smn2_presence()?;
        let mut has_smn1 = false;
        let mut has_smn2 = false;
        if smn1_nread_splice >= 2 {
            has_smn1 = true;
        }
        // reverse del_data
        self.del_data.reverse();

        let exon78_del_reads = self.del_data[0].del_reads_partial.clone();
        let exon78_del_nreads = exon78_del_reads.len();
        if smn2_nread_splice + exon78_del_nreads >= 2 {
            has_smn2 = true;
        }

        let tid = self.genome_tid().map(|x| x as i32).ok_or_else(|| {
            crate::phaser::Exception::new(format!(
                "Missing chromosome tid while running '{}' workflow",
                self.gene_name()
            ))
        })?;
        let init_read_hap_map = self.haplotypes_from_reads(
            None,
            /* kept_sites */ &hom_sites_to_add,
            Some(&add_sites),
            /* partial_deletion_reads */ None,
            (
                /* min_mapq= */ 5,
                /* check_clip= */ true,
                /* min_clip_len */ Some(50u32),
            ),
            tid,
            None,
            &hom_sites_to_add,
        )?;
        let (phase_results, known_del) =
            self.update_indel_and_phase(init_read_hap_map.clone(), &mut call)?;

        // rename haplotypes
        let mut assembled_haps = BTreeMap::new();
        let main_haps_clone = phase_results.assemblies.main_haps.clone();
        let mod_gene_name =
            intersperse(self.gene_name().split_terminator('-'), ",").collect::<String>();
        let mut smn1_haps = BTreeMap::new();
        let mut smn2_haps = BTreeMap::new();
        let mut smn_del_haps = BTreeMap::new();
        let mut counter_smn1 = 0;
        let mut counter_smn2 = 0;
        let mut counter_smndel = 0;
        let main_hap_strings = main_haps_clone
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>();
        let (smn1_assigned, smn2_assigned, smn_del_assigned) =
            self.assign_haps_to_gene(&main_hap_strings, has_smn1, smn2_nread_splice);
        let smn1_assigned = smn1_assigned.into_iter().collect::<HashSet<_>>();
        let smn2_assigned = smn2_assigned.into_iter().collect::<HashSet<_>>();
        let smn_del_assigned = smn_del_assigned.into_iter().collect::<HashSet<_>>();
        for hap in main_haps_clone.iter() {
            let hap_string = hap.to_string();
            if smn_del_assigned.contains(&hap_string) {
                counter_smndel += 1;
                let hap_name = format!("{mod_gene_name}_smndel78hap{}", counter_smndel);
                assembled_haps.insert(hap.vstr(), hap_name.clone());
                smn_del_haps.insert(hap_string, hap_name);
            } else if smn1_assigned.contains(&hap_string) {
                counter_smn1 += 1;
                let hap_name = format!("{mod_gene_name}_smn1hap{}", counter_smn1);
                assembled_haps.insert(hap.vstr(), hap_name.clone());
                smn1_haps.insert(hap_string, hap_name);
            } else if smn2_assigned.contains(&hap_string) {
                counter_smn2 += 1;
                let hap_name = format!("{mod_gene_name}_smn2hap{}", counter_smn2);
                assembled_haps.insert(hap.vstr(), hap_name.clone());
                smn2_haps.insert(hap_string, hap_name);
            }
        }
        log::debug!(
            "SMN assembled haplotype summary: total={}, smn1={}, smn2={}, smn_del78={}",
            assembled_haps.len(),
            smn1_haps.len(),
            smn2_haps.len(),
            smn_del_haps.len()
        );
        call.final_haplotypes = assembled_haps
            .clone()
            .into_iter()
            .map(|(k, v)| (k.to_string(), v))
            .collect::<BTreeMap<_, _>>();

        let mut smn1_cn: Option<i32> = None;
        let mut smn2_cn: Option<i32> = None;
        let mut smn_del_cn: i32 = counter_smndel;
        // homozygous case
        if self.init_het_sites.is_empty() && assembled_haps.is_empty() {
            if has_smn1 && !has_smn2 {
                smn1_cn = Some(2);
                smn2_cn = Some(0);
            } else if !has_smn1 && has_smn2 {
                smn1_cn = Some(0);
                smn2_cn = Some(2);
            }
        } else {
            smn1_cn = Some(counter_smn1);
            smn2_cn = Some(counter_smn2);
        }

        let mut two_cp_haps = Vec::new();
        if self.settings.targeted {
            let excluded_haplotypes = assembled_haps
                .keys()
                .filter(|hap| hap.contains(&b'3'))
                .map(|hap| VString::from(*hap))
                .collect::<Vec<_>>();
            two_cp_haps = self.compare_depth_by_read_count(
                &assembled_haps,
                &phase_results,
                0.05,
                &excluded_haplotypes,
            );
            for hap in &two_cp_haps {
                if hap.contains("smn1hap") {
                    if let Some(smn1_cn_value) = smn1_cn {
                        smn1_cn = Some(smn1_cn_value + 1);
                    }
                } else if hap.contains("smn2hap") {
                    if let Some(smn2_cn_value) = smn2_cn {
                        smn2_cn = Some(smn2_cn_value + 1);
                    }
                }
            }
            call.two_copy_haplotypes = two_cp_haps.clone();
        }
        if two_cp_haps.is_empty() {
            let two_cp_haps_candidate =
                self.compare_depth_by_read_count(&assembled_haps, &phase_results, 0.0005, &[]);
            let (updated_smn1_cn, mut two_cp_haps) = self.update_smn1_cn(
                smn1_cn,
                smn2_cn,
                smn1_nread_splice,
                smn2_nread_splice,
                has_smn2,
                &smn1_haps,
                &two_cp_haps_candidate,
                phase_results.assemblies.highest_cn > assembled_haps.len(),
            )?;
            smn1_cn = updated_smn1_cn;
            let (updated_smn2_cn, two_cp_haps_smn2) = self.update_smn2_cn(
                updated_smn1_cn,
                smn2_cn,
                smn1_nread_splice,
                smn2_nread_splice,
                exon78_del_nreads,
                &smn2_haps,
                &two_cp_haps_candidate,
            )?;
            smn2_cn = updated_smn2_cn;
            for hap in two_cp_haps_smn2 {
                two_cp_haps.push(hap);
            }
            call.two_copy_haplotypes = two_cp_haps;
        }

        // Output variants
        let haps =
            self.output_variants_in_haps(&phase_results, &known_del, assembled_haps.clone())?;
        // assign haplogroups
        let mut haplogroups = BTreeMap::new();
        for (hap_name, hap_info) in &haps {
            if hap_name.contains("del") {
                let haplogroup = Some(String::from("smn_del_exon78"));
                haplogroups.insert(hap_name, haplogroup);
            } else {
                let gene_name = &hap_name
                    .split_terminator('_')
                    .map(std::borrow::ToOwned::to_owned)
                    .collect::<Vec<_>>()[1]
                    .split_terminator('h')
                    .map(std::borrow::ToOwned::to_owned)
                    .collect::<Vec<_>>()[0];
                let haplogroup = self.assign_hap_to_group(
                    hap_info,
                    known_haps.get(gene_name).ok_or_else(|| {
                        crate::phaser::Exception::new(format!(
                            "Known haplotypes missing key '{}'",
                            gene_name
                        ))
                    })?,
                )?;
                haplogroups.insert(hap_name, haplogroup.clone());
                //log::trace!("hap {hap_name} haplogroup {haplogroup:?}");
            }
        }

        // homozygous case
        if assembled_haps.len() == 1 && self.init_het_sites.is_empty() {
            if has_smn1 && !has_smn2 {
                smn1_cn = Some(2);
                call.two_copy_haplotypes = smn1_haps.clone().into_values().collect::<Vec<String>>();
            } else if has_smn2 && !has_smn1 {
                if smn2_nread_splice > 0 {
                    smn2_cn = Some(2);
                    call.two_copy_haplotypes =
                        smn2_haps.clone().into_values().collect::<Vec<String>>();
                } else {
                    smn_del_cn = 2;
                    call.two_copy_haplotypes =
                        smn_del_haps.clone().into_values().collect::<Vec<String>>();
                }
            }
        }

        call.haplotype_details = haps
            .iter()
            .map(|(key, val)| (key.clone(), HapInfoForJson::from(val)))
            .collect::<BTreeMap<_, _>>();

        self.fill_in_call(phase_results, &mut call);
        call.region_specific_info
            .insert(String::from("smn1_cn"), smn1_cn.into());
        call.region_specific_info
            .insert(String::from("smn2_cn"), smn2_cn.into());
        call.region_specific_info
            .insert(String::from("smn_del78_cn"), smn_del_cn.into());
        call.region_specific_info
            .insert(String::from("smn1_read_number"), smn1_nread_splice.into());
        call.region_specific_info
            .insert(String::from("smn2_read_number"), smn2_nread_splice.into());
        call.region_specific_info.insert(
            String::from("smn_del78_read_number"),
            exon78_del_nreads.into(),
        );
        call.region_specific_info.insert(
            String::from("smn1_haplotypes"),
            serde_json::to_value(&smn1_haps)?,
        );
        call.region_specific_info.insert(
            String::from("smn2_haplotypes"),
            serde_json::to_value(&smn2_haps)?,
        );
        call.region_specific_info.insert(
            String::from("smn_del78_haplotypes"),
            serde_json::to_value(&smn_del_haps)?,
        );
        call.region_specific_info.insert(
            String::from("haplogroup"),
            serde_json::to_value(&haplogroups)?,
        );
        Ok(call)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config;
    use crate::depth;
    use crate::phaser;
    use crate::toolkit::util;
    use std::path::{Path, PathBuf};

    fn hg38_reference() -> Option<std::path::PathBuf> {
        std::env::var("HG38").ok().map(std::path::PathBuf::from)
    }

    fn stage_realigned_fixture(outdir: &Path, sample_id: &str, gene: &str) {
        std::fs::create_dir_all(outdir).expect("test outdir should be creatable");
        let src_bam = util::test_file("bams/HG00733.smn1.bam");
        let src_bai = util::test_file("bams/HG00733.smn1.bam.bai");
        let dest_bam = outdir.join(format!("{sample_id}_{gene}_realigned.bam"));
        let dest_bai = outdir.join(format!("{sample_id}_{gene}_realigned.bam.bai"));
        if !dest_bam.exists() {
            std::fs::copy(&src_bam, &dest_bam).expect("realigned BAM fixture should copy");
        }
        if !dest_bai.exists() {
            std::fs::copy(&src_bai, &dest_bai).expect("realigned BAM index should copy");
        }
    }

    fn temp_outdir() -> PathBuf {
        tempfile::TempDir::new()
            .expect("tempdir should build")
            .into_path()
    }

    fn build_smn1_phaser(outdir: PathBuf, depth: Option<depth::Result>) -> Option<Phaser> {
        let Some(hg38_reference) = hg38_reference() else {
            log::warn!(
                "Skipping SMN1 test because the HG38 environment variable is not configured."
            );
            return None;
        };
        stage_realigned_fixture(&outdir, "HG00733", "smn1");
        let settings = phaser::Settings::new(
            "HG00733",
            (hg38_reference, util::test_file("bams/HG00733.smn1.bam")),
            outdir,
            "smn1",
            &config::Region::try_load(None).expect("region config should load"),
            depth,
            None,
            String::from("38"),
            None,
            0.03,
            false,
        );
        let gene_config = config::Gene::try_load(None).expect("gene config should load");
        Some(Phaser::new(settings, Some(gene_config), None, None).expect("phaser should build"))
    }

    fn label_smn1_deletions(phaser: &mut Phaser) {
        phaser
            .parse_deletions_from_config()
            .expect("deletions should parse from config");
        phaser
            .label_big_dels()
            .expect("deletion labeling should succeed");
    }

    #[test]
    fn check_smn1_smn2_presence_reports_both() {
        let Some(mut phaser) = build_smn1_phaser(temp_outdir(), None) else {
            return;
        };
        let (smn1_reads, smn2_reads) = phaser
            .check_smn1_smn2_presence()
            .expect("splice-support check should succeed");
        assert!(smn1_reads > 0, "expected some SMN1 splice-supporting reads");
        assert!(smn2_reads > 0, "expected some SMN2 splice-supporting reads");
    }

    #[test]
    fn get_long_del_reads_matches_python_fixture_expectations() {
        let Some(mut phaser) = build_smn1_phaser(temp_outdir(), None) else {
            return;
        };
        label_smn1_deletions(&mut phaser);

        let smn1_del = phaser
            .del_data
            .iter()
            .find(|x| x.name() == "70917700_del_3559")
            .expect("smn1 configured deletion should exist");
        let smn2_del = phaser
            .del_data
            .iter()
            .find(|x| x.name() == "70948286_del_6310")
            .expect("smn2 configured deletion should exist");

        assert!(
            !smn2_del.del_reads.is_empty(),
            "expected some full-supporting SMN2 deletion reads"
        );
        assert!(
            !smn2_del.del_reads_partial.is_empty(),
            "expected some partial-supporting SMN2 deletion reads"
        );
        assert!(
            !smn2_del.del_negative_reads.is_empty(),
            "expected some negative SMN2 deletion reads"
        );
        assert!(
            smn1_del.del_reads.is_empty(),
            "expected no full-supporting SMN1 deletion reads in fixture"
        );
        assert!(
            smn1_del.del_reads_partial.is_empty(),
            "expected no partial-supporting SMN1 deletion reads in fixture"
        );
        assert!(
            !smn1_del.del_negative_reads.is_empty(),
            "expected some negative SMN1 deletion reads"
        );
    }

    #[test]
    fn allow_del_bases_matches_python_cases() {
        let Some(mut phaser) = build_smn1_phaser(temp_outdir(), None) else {
            return;
        };
        label_smn1_deletions(&mut phaser);

        assert!(phaser.allow_del_bases(70948287));
        assert!(!phaser.allow_del_bases(70948285));
    }

    #[test]
    fn assign_haps_to_gene_matches_python_cases() {
        let Some(phaser) = build_smn1_phaser(temp_outdir(), None) else {
            return;
        };

        let mut with_splice = phaser;
        with_splice.het_sites = vec![
            CandidateSite::new(70951939, String::from("A"), String::from("C")),
            CandidateSite::new(70951945, String::from("T"), String::from("G")),
            CandidateSite::new(70951957, String::from("T"), String::from("G")),
        ];
        let haps = vec![
            String::from("111"),
            String::from("121"),
            String::from("131"),
        ];
        let (smn1_haps, smn2_haps, smn2_del_haps) =
            with_splice.assign_haps_to_gene(&haps, true, 25);
        assert_eq!(smn1_haps, vec![String::from("111")]);
        assert_eq!(smn2_haps, vec![String::from("121")]);
        assert_eq!(smn2_del_haps, vec![String::from("131")]);

        let Some(mut no_splice_has_smn1) = build_smn1_phaser(temp_outdir(), None) else {
            return;
        };
        no_splice_has_smn1.het_sites = vec![
            CandidateSite::new(70951946, String::from("A"), String::from("C")),
            CandidateSite::new(70951948, String::from("T"), String::from("G")),
            CandidateSite::new(70951957, String::from("T"), String::from("G")),
        ];
        let (smn1_haps, smn2_haps, smn2_del_haps) =
            no_splice_has_smn1.assign_haps_to_gene(&haps, true, 0);
        assert_eq!(smn1_haps, vec![String::from("111"), String::from("121")]);
        assert!(smn2_haps.is_empty());
        assert_eq!(smn2_del_haps, vec![String::from("131")]);

        let Some(mut no_splice_has_smn2) = build_smn1_phaser(temp_outdir(), None) else {
            return;
        };
        no_splice_has_smn2.het_sites = vec![
            CandidateSite::new(70951946, String::from("A"), String::from("C")),
            CandidateSite::new(70951948, String::from("T"), String::from("G")),
            CandidateSite::new(70951957, String::from("T"), String::from("G")),
        ];
        let (smn1_haps, smn2_haps, smn2_del_haps) =
            no_splice_has_smn2.assign_haps_to_gene(&haps, false, 25);
        assert!(smn1_haps.is_empty());
        assert_eq!(smn2_haps, vec![String::from("111"), String::from("121")]);
        assert_eq!(smn2_del_haps, vec![String::from("131")]);

        let Some(mut no_splice_has_both) = build_smn1_phaser(temp_outdir(), None) else {
            return;
        };
        no_splice_has_both.het_sites = vec![
            CandidateSite::new(70951946, String::from("A"), String::from("C")),
            CandidateSite::new(70951948, String::from("T"), String::from("G")),
            CandidateSite::new(70951957, String::from("T"), String::from("G")),
        ];
        let (smn1_haps, smn2_haps, smn2_del_haps) =
            no_splice_has_both.assign_haps_to_gene(&haps, true, 25);
        assert!(smn1_haps.is_empty());
        assert!(smn2_haps.is_empty());
        assert_eq!(smn2_del_haps, vec![String::from("131")]);
    }

    #[test]
    fn update_smn1_cn_matches_key_python_cases() {
        let Some(mut phaser) = build_smn1_phaser(
            temp_outdir(),
            Some(depth::Result {
                median: 15.0,
                median_absolute_difference: 0.1,
                sex: crate::depth::Sex::Other,
            }),
        ) else {
            return;
        };
        let smn1_haps = BTreeMap::from([(String::from("11"), String::from("smn1hap1"))]);
        let (cn, two_cp_haps) = phaser
            .update_smn1_cn(Some(1), Some(1), 2, 2, true, &smn1_haps, &vec![], false)
            .expect("smn1 cn update should succeed");
        assert_eq!(cn, None);
        assert!(two_cp_haps.is_empty());

        let Some(mut phaser) = build_smn1_phaser(temp_outdir(), None) else {
            return;
        };
        let (cn, two_cp_haps) = phaser
            .update_smn1_cn(Some(1), Some(1), 3, 1, false, &smn1_haps, &vec![], true)
            .expect("smn1 cn update should succeed");
        assert_eq!(cn, Some(2));
        assert_eq!(two_cp_haps, vec![String::from("smn1hap1")]);

        let Some(mut phaser) = build_smn1_phaser(
            tempfile::TempDir::new()
                .expect("tempdir should build")
                .into_path(),
            None,
        ) else {
            return;
        };
        let (cn, two_cp_haps) = phaser
            .update_smn1_cn(Some(1), Some(2), 30, 32, true, &smn1_haps, &vec![], false)
            .expect("smn1 cn update should succeed");
        assert_eq!(cn, Some(2));
        assert_eq!(two_cp_haps, vec![String::from("smn1hap1")]);
    }

    #[test]
    fn update_smn2_cn_matches_key_python_cases() {
        let Some(mut phaser) = build_smn1_phaser(
            temp_outdir(),
            Some(depth::Result {
                median: 15.0,
                median_absolute_difference: 0.1,
                sex: crate::depth::Sex::Other,
            }),
        ) else {
            return;
        };
        let smn2_haps = BTreeMap::from([(String::from("22"), String::from("smn2hap1"))]);
        let (cn, two_cp_haps) = phaser
            .update_smn2_cn(Some(1), Some(1), 1, 1, 0, &smn2_haps, &vec![])
            .expect("smn2 cn update should succeed");
        assert_eq!(cn, None);
        assert!(two_cp_haps.is_empty());

        let Some(mut phaser) = build_smn1_phaser(
            tempfile::TempDir::new()
                .expect("tempdir should build")
                .into_path(),
            Some(depth::Result {
                median: 30.0,
                median_absolute_difference: 0.1,
                sex: crate::depth::Sex::Other,
            }),
        ) else {
            return;
        };
        let smn2_haps = BTreeMap::from([
            (String::from("22"), String::from("smn2hap1")),
            (String::from("21"), String::from("smn2hap2")),
        ]);
        let (cn, two_cp_haps) = phaser
            .update_smn2_cn(Some(0), Some(2), 0, 60, 0, &smn2_haps, &vec![])
            .expect("smn2 cn update should succeed");
        assert_eq!(cn, Some(4));
        let mut expected = vec![String::from("smn2hap1"), String::from("smn2hap2")];
        let mut found = two_cp_haps;
        expected.sort();
        found.sort();
        assert_eq!(found, expected);

        let Some(mut phaser) = build_smn1_phaser(
            tempfile::TempDir::new()
                .expect("tempdir should build")
                .into_path(),
            None,
        ) else {
            return;
        };
        let (cn, two_cp_haps) = phaser
            .update_smn2_cn(
                Some(1),
                Some(2),
                20,
                20,
                0,
                &smn2_haps,
                &vec![String::from("smn2hap2")],
            )
            .expect("smn2 cn update should succeed");
        assert_eq!(cn, Some(3));
        assert_eq!(two_cp_haps, vec![String::from("smn2hap2")]);
    }
}
