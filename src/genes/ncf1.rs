// NCF1 specific caller
use crate::io::json::GeneCall;
use crate::phaser::Phaser;
use crate::phaser::{Assignment, HapInfoForJson};
use crate::toolkit::math::depth_prob;
use crate::toolkit::site_selection::query_seq_counter;
use crate::toolkit::site_selection::CandidateSite;
use crate::toolkit::util::{DError, HashMap};
use itertools::intersperse;
use rust_htslib::bam::Read;
use std::collections::BTreeMap;
use std::str::FromStr;
use vstr::VString;

impl Phaser {
    pub fn run_ncf1(&mut self) -> Result<GeneCall, DError> {
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

        let (hom_sites_to_add, add_sites) = self.get_sites(&seq, None, None)?;
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
        let uniquely_supporting_reads = &phase_results.uniquely_supporting_reads;

        // main variant is 74777266_GGT_G
        let pivot_var = self
            .locus_config()
            .get("pivot_var")
            .and_then(|x| x.as_str())
            .ok_or_else(|| {
                crate::phaser::Exception::new(format!(
                    "Missing pivot_var in region config for gene '{}'",
                    self.gene_name()
                ))
            })?
            .to_string();
        let pivot_var_site = CandidateSite::from_str(&pivot_var)?;
        let ref_seq = {
            let faidx = self.make_faidx()?;
            let (chrom, start, stop) = self.parsed_nchr_0based().ok_or_else(|| {
                crate::phaser::Exception::new(format!(
                    "Malformed realign region '{}'; expected 'chr:start-end' format",
                    self.realign_region
                ))
            })?;

            faidx
                .fetch_seq(chrom, start as usize, stop as usize)?
                .to_owned()
        };
        let mut gene_reads = 0;
        let mut pseudo_reads = 0;
        let var_size = pivot_var_site.var_seq.len() as i64 - pivot_var_site.ref_seq.len() as i64;
        let pseudo_token = if var_size < 0 {
            format!(
                "{}{}{}",
                std::str::from_utf8(&[pivot_var_site.ref_seq[0]])?,
                var_size,
                VString::from(&pivot_var_site.ref_seq[1..])
            )
        } else {
            pivot_var_site.var_seq.to_string()
        };
        if let Some(pivot_site) = self.pivot_site {
            let mut bam = self.try_realigned_bam()?;
            let tid = self.genome_tid().map(|x| x as i32).ok_or_else(|| {
                crate::phaser::Exception::new(format!(
                    "Missing chromosome tid while counting NCF1 pivot pileup reads for gene '{}'",
                    self.gene_name()
                ))
            })?;
            bam.fetch((tid, pivot_site, pivot_site + 1))?;
            let mut aln2seq = HashMap::default();
            for p in bam.pileup() {
                let p = p?;
                if i64::from(p.pos()) != pivot_site {
                    continue;
                }
                let (counts, _first_seen) = query_seq_counter(
                    &p,
                    &ref_seq,
                    &self.settings.site_selection_settings,
                    self.offset(),
                    &mut aln2seq,
                );
                for (seq, count) in counts {
                    if seq.as_slice() == b"G" {
                        gene_reads += count;
                    } else if seq.as_slice() == pseudo_token.as_bytes() {
                        pseudo_reads += count;
                    }
                }
                break;
            }
        }
        let pivot_var_reads =
            self.check_variants_in_haplotypes(&pivot_var_site, &ref_seq, Some(13))?;
        log::debug!(
            "NCF1 pivot variant read assignments: site={pivot_var_site:?}, total_reads={}",
            pivot_var_reads.len()
        );

        // first name haplotypes
        let mut assembled_haps = BTreeMap::new();
        let main_haps_clone = phase_results.assemblies.main_haps.clone();
        let mod_gene_name =
            intersperse(self.gene_name().split_terminator('-'), ",").collect::<String>();
        for (idx, hap) in main_haps_clone.iter().enumerate() {
            assembled_haps.insert(hap.vstr(), format!("{mod_gene_name}_hap{}", idx + 1));
        }

        // Output variants
        let haps =
            self.output_variants_in_haps(&phase_results, &known_del, assembled_haps.clone())?;

        // rename haplotypes
        let mut renamed_haps = BTreeMap::new();
        let mut counter_gene = 0;
        let mut counter_pseudo = 0;
        for (hap_seq, hap_name) in &assembled_haps {
            let hap_vars = &haps
                .get(hap_name)
                .ok_or_else(|| {
                    crate::phaser::Exception::new(format!(
                        "Haplotype '{}' missing from output_variants_in_haps map",
                        hap_name
                    ))
                })?
                .variants;
            if !hap_vars.contains(&pivot_var_site) {
                let hap_reads = uniquely_supporting_reads
                    .get(&VString::from(hap_seq))
                    .ok_or_else(|| {
                        crate::phaser::Exception::new(format!(
                            "Assembled hap sequence '{}' missing from uniquely_supporting_reads",
                            hap_seq
                        ))
                    })?;
                let mut this_hap_reads_check_var = Vec::new();
                for read in hap_reads {
                    if pivot_var_reads.contains_key(read) {
                        this_hap_reads_check_var.push(pivot_var_reads[read]);
                    }
                }
                let total_count = this_hap_reads_check_var.len();
                let alt_count = this_hap_reads_check_var
                    .iter()
                    .filter(|x| **x == Assignment::Alt)
                    .collect::<Vec<_>>()
                    .len();
                if alt_count as f32 > (total_count as f32) * 0.7 {
                    counter_pseudo += 1;
                    renamed_haps.insert(
                        hap_name,
                        format!("{mod_gene_name}_pseudohap{}", counter_pseudo),
                    );
                } else {
                    counter_gene += 1;
                    renamed_haps.insert(hap_name, format!("{mod_gene_name}_hap{}", counter_gene));
                }
            } else {
                counter_pseudo += 1;
                renamed_haps.insert(
                    hap_name,
                    format!("{mod_gene_name}_pseudohap{}", counter_pseudo),
                );
            }
        }
        let mut assembled_haps_renamed = BTreeMap::new();
        for (hap_seq, hap_name) in &assembled_haps {
            assembled_haps_renamed.insert(
                *hap_seq,
                renamed_haps
                    .get(hap_name)
                    .ok_or_else(|| {
                        crate::phaser::Exception::new(format!(
                            "Haplotype '{}' missing from renamed_haps map",
                            hap_name
                        ))
                    })?
                    .to_string(),
            );
        }
        call.final_haplotypes = assembled_haps_renamed
            .clone()
            .into_iter()
            .map(|(k, v)| (k.to_string(), v))
            .collect::<BTreeMap<_, _>>();

        let mut hap_details_renamed = BTreeMap::new();
        for (hap_name, hap_info) in &haps {
            hap_details_renamed.insert(
                renamed_haps
                    .get(hap_name)
                    .ok_or_else(|| {
                        crate::phaser::Exception::new(format!(
                            "Haplotype '{}' missing from renamed_haps map",
                            hap_name
                        ))
                    })?
                    .to_string(),
                hap_info.clone(),
            );
        }
        call.haplotype_details = hap_details_renamed
            .iter()
            .map(|(key, val)| (key.clone(), HapInfoForJson::from(val)))
            .collect::<BTreeMap<_, _>>();

        let mut total_cn = assembled_haps_renamed.len();
        let mut two_cp_haps: Vec<String> = Vec::new();
        let mut gene_cn: Option<i32> = Some(counter_gene);
        let mut total_cn_is_no_call = false;
        // scenario where only three haplotypes are found, possibly each at CN2
        if total_cn == 3 {
            if let Some(depth) = self.settings.depth.as_ref() {
                let region_median_depth = self.region_avg_depth[0].0.ceil() as i32;
                let prob = depth_prob(region_median_depth, depth.median);
                if let Some(prob_value) = prob {
                    if prob_value[2] + prob_value[3] > 0.95 {
                        two_cp_haps = assembled_haps_renamed.values().cloned().collect::<Vec<_>>();
                        for hap in &two_cp_haps {
                            total_cn += 1;
                            if hap.contains("pseudo") {
                                counter_pseudo += 1;
                            } else {
                                counter_gene += 1;
                            }
                        }
                    }
                }
            }
            if counter_gene == 1 && counter_pseudo == 2 {
                call.total_cn = None;
                gene_cn = None;
                total_cn_is_no_call = true;
            }
        } else if counter_gene == 1 {
            let two_cp_hap_candidate =
                self.compare_depth(&hap_details_renamed, &assembled_haps_renamed, false, false)?;
            for hap in &two_cp_hap_candidate {
                if hap.contains("ncf1_hap1") {
                    counter_gene += 1;
                    total_cn += 1;
                    two_cp_haps = two_cp_hap_candidate.clone();
                }
            }
        }

        // check against genome depth
        if gene_cn.is_some() {
            if let Some(depth) = self.settings.depth.as_ref() {
                let prob = depth_prob(gene_reads, depth.median / 2.0_f64);
                if let Some(prob_value) = prob {
                    if prob_value[0] < 0.9 && counter_gene == 1 {
                        gene_cn = None;
                        call.total_cn = None;
                        total_cn_is_no_call = true;
                    }
                    if prob_value[0] > 0.95 && counter_gene > 1 && !two_cp_haps.is_empty() {
                        gene_cn = None;
                        call.total_cn = None;
                        total_cn_is_no_call = true;
                    }
                }
            }
        }
        if gene_cn.is_some() {
            gene_cn = Some(counter_gene);
        }
        if !total_cn_is_no_call {
            call.total_cn = if total_cn > 0 {
                Some(total_cn as i32)
            } else {
                None
            };
        }

        call.two_copy_haplotypes = two_cp_haps;
        self.fill_in_call(phase_results, &mut call);

        // additional fields to report
        call.region_specific_info
            .insert(String::from("gene_reads"), gene_reads.into());
        call.region_specific_info
            .insert(String::from("pseudo_reads"), pseudo_reads.into());
        call.region_specific_info
            .insert(String::from("gene_cn"), gene_cn.into());
        Ok(call)
    }
}
