// HBA specific caller
use crate::io::json::GeneCall;
use crate::phaser::HapInfoForJson;
use crate::phaser::Phaser;
use crate::toolkit::math::depth_prob;
use crate::toolkit::util::DError;
use itertools::intersperse;
use std::collections::BTreeMap;
use vstr::VString;

impl Phaser {
    fn get_surrounding_depth(&self) -> f32 {
        let depth_region = self.locus_config().depth_region();
        let Ok(mut bam) = self.try_genome_bam() else {
            log::warn!(
                "Failed to open genome BAM while computing surrounding depth for gene {}",
                self.gene_name()
            );
            return f32::NAN;
        };
        let Some(tid) = self.genome_tid().map(|x| x as i32) else {
            log::warn!(
                "Missing genome tid while computing surrounding depth for gene {}",
                self.gene_name()
            );
            return f32::NAN;
        };
        let region_depth = Self::regional_depth(
            &mut bam,
            tid,
            &depth_region,
            /* num_intervals (step) */ None,
            /* exclude_flag */ None,
            /* one_based */ Some(true),
            /* percentile */ Some(crate::phaser::region_depth::PERCENTILE),
        );
        let (surrounding_region_depth, _percentile) = region_depth[0];
        surrounding_region_depth
    }

    /// Run HBA-specific phasing/copy-number workflow with clip-signature-based
    /// haplotype naming and structural-event interpretation.
    pub fn run_hba(&mut self) -> Result<GeneCall, DError> {
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

        let surrounding_region_depth = self.get_surrounding_depth();

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
                /* min_clip_len */ Some(9u32),
            ),
            tid,
            None,
            &hom_sites_to_add,
        )?;
        let (mut phase_results, known_del) =
            self.update_indel_and_phase(init_read_hap_map.clone(), &mut call)?;

        // rename haplotypes
        let mut assembled_haps = BTreeMap::new();
        let main_haps_clone = phase_results.assemblies.main_haps.clone();
        let mod_gene_name =
            intersperse(self.gene_name().split_terminator('-'), ",").collect::<String>();

        let first_clip_3p = self.clip_3p_positions.first().ok_or_else(|| {
            crate::phaser::Exception::new(format!(
                "clip_3p_positions is empty for gene '{}'",
                self.gene_name()
            ))
        })?;
        let second_clip_3p = self.clip_3p_positions.last().ok_or_else(|| {
            crate::phaser::Exception::new(format!(
                "clip_3p_positions is empty for gene '{}'",
                self.gene_name()
            ))
        })?;
        let first_clip_5p = self.clip_5p_positions.first().ok_or_else(|| {
            crate::phaser::Exception::new(format!(
                "clip_5p_positions is empty for gene '{}'",
                self.gene_name()
            ))
        })?;
        let mut count_hba1 = 0;
        let mut count_hba2 = 0;
        let mut count_unknown = 0;
        let mut count_3p7del = 0;
        let mut count_3p7dup = 0;
        let mut count_4p2del = 0;
        let mut count_4p2dup = 0;
        let mut count_homology = 0;
        for hap in main_haps_clone.iter() {
            let clip_5p = self.get_5pclip_from_hap(&hap.vstr())?;
            let clip_3p = self.get_3pclip_from_hap(&hap.vstr())?;
            log::debug!(
                "Classifying HBA haplotype by clip support: hap={}, clip_5p={clip_5p:?}, clip_3p={clip_3p:?}",
                hap
            );
            if clip_5p.is_none() || clip_3p.is_none() {
                count_unknown += 1;
                assembled_haps.insert(
                    hap.vstr(),
                    format!("{mod_gene_name}_unknownhap{}", count_unknown),
                );
            }
            if let Some(clip_5p_value) = clip_5p {
                if let Some(clip_3p_value) = clip_3p {
                    if clip_5p_value == 0 && clip_3p_value == 0 {
                        count_hba2 += 1;
                        assembled_haps
                            .insert(hap.vstr(), format!("{mod_gene_name}_hba2hap{}", count_hba2));
                    } else if clip_5p_value == 0 {
                        if clip_3p_value == *first_clip_3p {
                            count_3p7del += 1;
                            assembled_haps.insert(
                                hap.vstr(),
                                format!("{mod_gene_name}_3p7delhap{}", count_3p7del),
                            );
                        } else if clip_3p_value == *second_clip_3p {
                            count_4p2dup += 1;
                            assembled_haps.insert(
                                hap.vstr(),
                                format!("{mod_gene_name}_4p2duphap{}", count_4p2dup),
                            );
                        }
                    } else if clip_3p_value == 0 {
                        if clip_5p_value == *first_clip_5p {
                            count_3p7dup += 1;
                            assembled_haps.insert(
                                hap.vstr(),
                                format!("{mod_gene_name}_3p7duphap{}", count_3p7dup),
                            );
                        } else {
                            count_4p2del += 1;
                            assembled_haps.insert(
                                hap.vstr(),
                                format!("{mod_gene_name}_4p2delhap{}", count_4p2del),
                            );
                        }
                    } else if clip_5p_value == *first_clip_5p && clip_3p_value == *first_clip_3p {
                        count_hba1 += 1;
                        assembled_haps
                            .insert(hap.vstr(), format!("{mod_gene_name}_hba1hap{}", count_hba1));
                    } else if clip_5p_value != *first_clip_5p && clip_3p_value == *second_clip_3p {
                        count_homology += 1;
                        assembled_haps.insert(
                            hap.vstr(),
                            format!("{mod_gene_name}_homologyhap{}", count_homology),
                        );
                    }
                }
            }
        }
        call.final_haplotypes = assembled_haps
            .clone()
            .into_iter()
            .map(|(k, v)| (k.to_string(), v))
            .collect::<BTreeMap<_, _>>();
        // Phase alleles
        let mut haps_to_exclude: Vec<VString> = Vec::new();
        let mut hap_names_to_consider = BTreeMap::new();
        for (hap_seq, hap_name) in &assembled_haps {
            if hap_name.contains("homology") || hap_name.contains("unknown") {
                haps_to_exclude.push(hap_seq.into());
            } else {
                hap_names_to_consider.insert(*hap_seq, hap_name.clone());
            };
        }
        let mut alleles = Vec::new();
        if self.config.locus.to_phase() {
            let allele_result =
                self.phase_alleles(&mut phase_results, &assembled_haps, Some(haps_to_exclude));
            alleles = allele_result.alleles;
            if self.all_haps_phased_onto_one_allele(&alleles, &hap_names_to_consider)? {
                alleles = Vec::new();
            }
            call.region_specific_info
                .insert(String::from("alleles_final"), alleles.clone().into());
            call.region_specific_info.insert(
                String::from("haplotype_links"),
                serde_json::to_value(&allele_result.haplotype_links)?,
            );
            // use final alleles for YC tag in BAM, to best show 2/1 case
            call.region_specific_info
                .insert(String::from("raw_alleles"), alleles.clone().into());
        }
        // Output variants
        let haps =
            self.output_variants_in_haps(&phase_results, &known_del, assembled_haps.clone())?;
        call.haplotype_details = haps
            .iter()
            .map(|(key, val)| (key.clone(), HapInfoForJson::from(val)))
            .collect::<BTreeMap<_, _>>();

        let mut sv_called = BTreeMap::new();
        let mut sv_coordinate_4p2 = self
            .locus_config()
            .get("4p2_coordinate")
            .and_then(|x| x.as_str())
            .ok_or_else(|| {
                crate::phaser::Exception::new(format!(
                    "Missing 4p2_coordinate in region config for gene '{}'",
                    self.gene_name()
                ))
            })?
            .to_string();
        let mut sv_coordinate_3p7 = self
            .locus_config()
            .get("3p7_coordinate")
            .and_then(|x| x.as_str())
            .ok_or_else(|| {
                crate::phaser::Exception::new(format!(
                    "Missing 3p7_coordinate in region config for gene '{}'",
                    self.gene_name()
                ))
            })?
            .to_string();
        if self.settings.genome == "37" {
            sv_coordinate_4p2 = sv_coordinate_4p2
                .strip_prefix("chr")
                .unwrap_or(&sv_coordinate_4p2)
                .to_string();
            sv_coordinate_3p7 = sv_coordinate_3p7
                .strip_prefix("chr")
                .unwrap_or(&sv_coordinate_3p7)
                .to_string();
        }
        for hap in assembled_haps.values() {
            if hap.contains("del") || hap.contains("dup") {
                let Some(hap_no_prefix) = hap.strip_prefix("hba_") else {
                    continue;
                };
                let Some(sv_name) = hap_no_prefix.split_terminator("hap").next() else {
                    continue;
                };
                let mut sv_type = None;
                if hap.contains("del") {
                    sv_type = sv_name.strip_suffix("del");
                } else if hap.contains("dup") {
                    sv_type = sv_name.strip_suffix("dup");
                }
                if let Some(sv_type) = sv_type {
                    if sv_type == "3p7" {
                        sv_called.insert(
                            hap.to_string(),
                            format!("{sv_name},{}", sv_coordinate_3p7.clone()),
                        );
                    } else if sv_type == "4p2" {
                        sv_called.insert(
                            hap.to_string(),
                            format!("{sv_name},{}", sv_coordinate_4p2.clone()),
                        );
                    }
                }
            }
        }
        let mut two_cp_haps = Vec::new();
        if count_3p7del == 0 && count_4p2del == 1 {
            if count_hba1 == 1 && count_hba2 == 1 {
                for hap in assembled_haps.values() {
                    if hap.contains("hba1") {
                        two_cp_haps.push(hap.to_string());
                        count_hba1 += 1;
                    }
                }
            }
        } else if count_3p7del == 0 && count_4p2del == 0 {
            if count_hba1 == 1 && count_hba2 == 2 {
                two_cp_haps = assembled_haps
                    .values()
                    .filter(|a| a.contains("hba1"))
                    .map(|a| a.to_string())
                    .collect::<Vec<_>>();
                count_hba1 += 1;
            } else if count_hba1 == 2 && count_hba2 == 1 {
                two_cp_haps = assembled_haps
                    .values()
                    .filter(|a| a.contains("hba2"))
                    .map(|a| a.to_string())
                    .collect::<Vec<_>>();
                count_hba2 += 1;
            } else if count_hba1 == 1 && count_hba2 == 1 {
                if let Some(depth) = self.settings.depth.as_ref() {
                    let probs =
                        depth_prob(surrounding_region_depth as i32, (depth.median as f32) / 2.0);
                    log::debug!("HBA surrounding-depth copy-number probabilities: {probs:?}");
                    if let Some(probs_value) = probs {
                        if probs_value[0] < 0.999 {
                            two_cp_haps = assembled_haps
                                .values()
                                .filter(|a| a.contains("hba2") || a.contains("hba1"))
                                .map(|a| a.to_string())
                                .collect::<Vec<_>>();
                            count_hba1 += 1;
                            count_hba2 += 1;
                        }
                    }
                } else {
                    // assume two pairs of identical copies
                    two_cp_haps = assembled_haps
                        .values()
                        .filter(|a| a.contains("hba2") || a.contains("hba1"))
                        .map(|a| a.to_string())
                        .collect::<Vec<_>>();
                    count_hba1 += 1;
                    count_hba2 += 1;
                }
            }
        } else if count_3p7del == 1
            && count_hba1 == 0
            && count_hba2 == 0
            && count_3p7dup == 0
            && count_4p2del == 0
            && count_4p2dup == 0
        {
            count_3p7del += 1;
            two_cp_haps = assembled_haps
                .values()
                .filter(|a| a.contains("del"))
                .map(|a| a.to_string())
                .collect::<Vec<_>>();
        }

        let mut total_cn = assembled_haps.len() + two_cp_haps.len()
            - count_homology
            - count_4p2del
            - count_unknown;
        if self.init_het_sites.is_empty() && total_cn < 2 {
            total_cn = 2;
        }
        // genotype
        let mut genotype: Option<String> = None;
        if count_4p2del == 0 && count_4p2dup == 0 {
            if count_3p7del == 1 && total_cn == 3 {
                genotype = Some(String::from("-a/aa"));
            } else if count_3p7del == 1 && total_cn == 4 && count_3p7dup == 1 {
                genotype = Some(String::from("-a/aaa"));
            } else if count_3p7del == 2 && total_cn == 2 {
                genotype = Some(String::from("-a/-a"));
            } else if count_3p7del == 0 && total_cn == 4 && count_3p7dup == 0 {
                genotype = Some(String::from("aa/aa"));
            } else if count_3p7del == 0 && total_cn == 5 && count_3p7dup == 1 {
                genotype = Some(String::from("aa/aaa"));
            } else if count_3p7del == 0 && total_cn == 6 && count_3p7dup == 2 {
                genotype = Some(String::from("aaa/aaa"));
            } else if count_3p7del == 0 && total_cn == 2 && count_3p7dup == 0 {
                genotype = Some(String::from("--/aa"));
            }
        }
        // 4.2
        else if count_4p2del > 0 && count_4p2dup == 0 {
            if count_3p7del == 1 && count_hba1 == 1 && count_hba2 == 0 {
                // 3.7/4.2
                genotype = Some(String::from("-a/-a"));
            } else if count_hba2 == 1 && count_hba1 == 2 {
                if count_3p7dup == 0 {
                    genotype = Some(String::from("-a/aa"));
                } else {
                    // anti3.7/4.2
                    if check_two_haps_in_cis(&alleles, "4p2del", "3p7dup") {
                        genotype = Some(String::from("aa/aa"));
                    } else {
                        genotype = Some(String::from("-a/aaa"));
                    }
                }
            } else if count_hba2 == 0 && count_hba1 == 2 {
                genotype = Some(String::from("-a/-a"));
            }
        } else if count_4p2del == 0 && count_4p2dup > 0 {
            // we dont consider anti4.2/anti4.2 for now. Should be very rare
            if count_3p7del == 1 && count_hba1 == 1 && count_hba2 == 1 {
                // 3.7/anti4.2
                if check_two_haps_in_cis(&alleles, "4p2dup", "3p7del") {
                    genotype = Some(String::from("aa/aa"));
                } else {
                    genotype = Some(String::from("aaa/-a"));
                }
            } else if count_hba2 == 2 && count_hba1 == 2 {
                if count_3p7dup == 0 {
                    genotype = Some(String::from("aaa/aa"));
                } else {
                    // anti3.7/anti4.2
                    if check_two_haps_in_cis(&alleles, "4p2dup", "3p7dup") {
                        genotype = Some(String::from("aaaa/aa"));
                    } else {
                        genotype = Some(String::from("aaa/aaa"));
                    }
                }
            }
        } else if count_4p2del > 0 && count_4p2dup > 0 {
            // 4.2/anti4.2
            if count_hba1 == 2 && count_hba2 == 1 {
                genotype = Some(String::from("aaa/-a"));
            }
        }

        call.total_cn = if total_cn <= 1 {
            None
        } else {
            Some(total_cn as i32)
        };
        call.two_copy_haplotypes = two_cp_haps;
        // additional fields to report
        call.region_specific_info.insert(
            String::from("surrounding_region_depth"),
            surrounding_region_depth.into(),
        );
        call.region_specific_info
            .insert(String::from("genotype"), genotype.into());
        call.region_specific_info
            .insert(String::from("sv_called"), serde_json::to_value(&sv_called)?);
        // report
        self.fill_in_call(phase_results, &mut call);
        Ok(call)
    }
}

/// check if two types of haplotypes are on the same allele
/// Whether both haplotype names appear on the same phased allele.
#[must_use]
pub fn check_two_haps_in_cis(alleles: &Vec<Vec<String>>, name1: &str, name2: &str) -> bool {
    let mut in_cis = false;
    for allele in alleles {
        let mut found1 = false;
        let mut found2 = false;
        for hap in allele {
            if hap.contains(name1) {
                found1 = true;
            }
            if hap.contains(name2) {
                found2 = true;
            }
        }
        if found1 && found2 {
            in_cis = true;
        }
    }
    in_cis
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_check_two_haps_in_cis() {
        let alleles = vec![
            vec![String::from("hba_hba1hap1"), String::from("hba_hba1hap2")],
            vec![String::from("hba_hba1hap3"), String::from("hba_hba1hap4")],
        ];
        let cis = check_two_haps_in_cis(&alleles, "4p2", "3p7");
        assert!(!cis);

        let alleles = vec![
            vec![
                String::from("hba_3p7delhap1"),
                String::from("hba_3p7delhap2"),
            ],
            vec![
                String::from("hba_4p2delhap1"),
                String::from("hba_4p2delhap2"),
            ],
        ];
        let cis = check_two_haps_in_cis(&alleles, "4p2", "3p7");
        assert!(!cis);

        let alleles = vec![
            vec![
                String::from("hba_3p7delhap1"),
                String::from("hba_4p2delhap2"),
            ],
            vec![
                String::from("hba_4p2delhap1"),
                String::from("hba_4p2delhap2"),
            ],
        ];
        let cis = check_two_haps_in_cis(&alleles, "4p2", "3p7");
        assert!(cis);
    }
}
