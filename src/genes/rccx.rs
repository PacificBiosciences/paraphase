// RCCX specific caller
use crate::io::json::GeneCall;
use crate::phaser::Phaser;
use crate::phaser::{HapInfo, HapInfoForJson};
use crate::toolkit::site_selection::CandidateSite;
use crate::toolkit::util::DError;
use itertools::intersperse;
use std::collections::BTreeMap;
use vstr::VStr;

impl Phaser {
    /// update alleles based on all info available
    #[allow(clippy::type_complexity)]
    fn update_alleles(
        &mut self,
        new_alleles: &Vec<Vec<String>>,
        haplotypes: &BTreeMap<String, HapInfo>,
        assembled_haps: &BTreeMap<VStr<'_>, String>,
        single_copies: &Vec<String>,
        starting_copies: &Vec<String>,
        ending_copies: &Vec<String>,
        hcn: usize,
    ) -> Result<(bool, Vec<Vec<String>>, Vec<String>), DError> {
        let final_haps = &assembled_haps.values().cloned().collect::<Vec<String>>();
        let mut two_cp_haplotypes = Vec::new();
        let mut successful_phasing = false;
        let mut updated_alleles = Vec::new();
        let nhap = final_haps.len();
        let nsingle = single_copies.len();
        // update the case with homozygous deletion
        if nhap == 1 && nsingle == 1 {
            two_cp_haplotypes = final_haps.clone();
            successful_phasing = true;
            if let Some(hap) = single_copies.first() {
                let hap = hap.to_string();
                updated_alleles = vec![vec![hap.clone()], vec![hap.clone()]];
            }
        } else if nsingle == 1 && nhap < 5 {
            // the deletion haplotype will be reported as an allele
            let mut ok_to_phase = false;
            if new_alleles.len() == 1 && !new_alleles.contains(single_copies) {
                if let Some(first_allele) = new_alleles.first() {
                    if first_allele.len() == nhap - 1 {
                        updated_alleles = vec![first_allele.to_vec(), single_copies.to_vec()];
                    } else if first_allele.len() < nhap - 1
                        && starting_copies.len() == 1
                        && ending_copies.len() == 1
                    {
                        ok_to_phase = true;
                    }
                }
            }
            if new_alleles.is_empty() || ok_to_phase {
                let remaining_hap = final_haps
                    .iter()
                    .filter(|x| !single_copies.contains(*x))
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>();
                if remaining_hap.len() == nhap - 1 {
                    updated_alleles = vec![single_copies.clone(), remaining_hap];
                }
            }
        } else if nsingle == 2
            && new_alleles.is_empty()
            && nhap == 2
            && starting_copies.is_empty()
            && ending_copies.is_empty()
        {
            // deletions on each allele
            if let (Some(first_single), Some(second_single)) =
                (single_copies.first(), single_copies.last())
            {
                updated_alleles = vec![
                    vec![first_single.to_string()],
                    vec![second_single.to_string()],
                ];
                successful_phasing = true;
            }
        } else if single_copies.is_empty() {
            // homozygous, each haplotype has cn 2
            if nhap == 2
                && ending_copies.len() == 1
                && starting_copies.len() == 1
                && new_alleles.len() == 1
            {
                two_cp_haplotypes = final_haps.clone();
                successful_phasing = true;
                if let Some(first_allele) = new_alleles.first() {
                    let first_allele = first_allele.to_vec();
                    updated_alleles = vec![first_allele.clone(), first_allele.clone()];
                }
            }
            // depth-based adjustment when found 3 haplotypes or <2 ending haplotypes
            let two_cp_hap_candidate =
                self.compare_depth(haplotypes, assembled_haps, true, false)?;
            if ending_copies.len() == 1 && starting_copies.len() == 2 {
                if let Some(ending_copy) = ending_copies.first() {
                    if two_cp_hap_candidate.len() == 1 && two_cp_hap_candidate.contains(ending_copy)
                    {
                        two_cp_haplotypes = two_cp_hap_candidate.clone();
                        if nhap == 3 {
                            if let (Some(first_starting), Some(second_starting)) =
                                (starting_copies.first(), starting_copies.last())
                            {
                                updated_alleles = vec![
                                    vec![first_starting.to_string(), ending_copy.to_string()],
                                    vec![second_starting.to_string(), ending_copy.to_string()],
                                ];
                                successful_phasing = true;
                            }
                        }
                    }
                }
            } else if ending_copies.len() == 2 && starting_copies.len() == 1 {
                if let Some(starting_copy) = starting_copies.first() {
                    if two_cp_hap_candidate.len() == 1
                        && two_cp_hap_candidate.contains(starting_copy)
                    {
                        two_cp_haplotypes = two_cp_hap_candidate.clone();
                        if nhap == 3 {
                            if let (Some(first_ending), Some(second_ending)) =
                                (ending_copies.first(), ending_copies.last())
                            {
                                updated_alleles = vec![
                                    vec![starting_copy.to_string(), first_ending.to_string()],
                                    vec![starting_copy.to_string(), second_ending.to_string()],
                                ];
                                successful_phasing = true;
                            }
                        }
                    }
                }
            }
            // add missing links when there is no two-cp haplotypes
            if two_cp_haplotypes.is_empty() && ending_copies.len() <= 2 {
                // add the missing link in cn=4
                if (nhap == 3 || nhap == 4) && new_alleles.len() == 1 && hcn == nhap {
                    if let Some(first_allele) = new_alleles.first() {
                        if first_allele.len() == 2 {
                            let remaining_hap = final_haps
                                .iter()
                                .filter(|x| !first_allele.contains(*x))
                                .map(|x| x.to_string())
                                .collect::<Vec<_>>();
                            if remaining_hap.len() == nhap - 2 {
                                updated_alleles = vec![first_allele.to_vec(), remaining_hap];
                            }
                        }
                    }
                }
                // add the missing link in cn=5
                if nhap == 5 {
                    if new_alleles.len() == 1 {
                        if let Some(first_allele) = new_alleles.first() {
                            if first_allele.len() == 2 {
                                if let (
                                    Some(first_allele_first_hap),
                                    Some(first_allele_second_hap),
                                ) = (first_allele.first(), first_allele.last())
                                {
                                    if (starting_copies.contains(first_allele_first_hap)
                                        && ending_copies.contains(first_allele_second_hap))
                                        || (starting_copies.contains(first_allele_second_hap)
                                            && ending_copies.contains(first_allele_first_hap))
                                    {
                                        let remaining_hap = final_haps
                                            .iter()
                                            .filter(|x| !first_allele.contains(*x))
                                            .map(|x| x.to_string())
                                            .collect::<Vec<_>>();
                                        if remaining_hap.len() == 3 {
                                            updated_alleles =
                                                vec![first_allele.to_vec(), remaining_hap];
                                        }
                                    }
                                }
                            }
                        }
                    } else if new_alleles.len() == 2 {
                        if let (Some(first_allele), Some(second_allele)) =
                            (new_alleles.first(), new_alleles.last())
                        {
                            if let (
                                Some(first_allele_first_hap),
                                Some(first_allele_second_hap),
                                Some(second_allele_first_hap),
                                Some(second_allele_second_hap),
                            ) = (
                                first_allele.first(),
                                first_allele.last(),
                                second_allele.first(),
                                second_allele.last(),
                            ) {
                                let allele1 = (starting_copies.contains(first_allele_first_hap)
                                    && ending_copies.contains(first_allele_second_hap))
                                    || starting_copies.contains(first_allele_second_hap)
                                        && ending_copies.contains(first_allele_first_hap);
                                let allele2 = (starting_copies.contains(second_allele_first_hap)
                                    && ending_copies.contains(second_allele_second_hap))
                                    || starting_copies.contains(second_allele_second_hap)
                                        && ending_copies.contains(second_allele_first_hap);
                                if allele1 && !allele2 {
                                    let remaining_hap = final_haps
                                        .iter()
                                        .filter(|x| !first_allele.contains(*x))
                                        .map(|x| x.to_string())
                                        .collect::<Vec<_>>();
                                    if remaining_hap.len() == 3 {
                                        updated_alleles =
                                            vec![first_allele.to_vec(), remaining_hap];
                                    }
                                } else if allele2 && !allele1 {
                                    let remaining_hap = final_haps
                                        .iter()
                                        .filter(|x| !second_allele.contains(*x))
                                        .map(|x| x.to_string())
                                        .collect::<Vec<_>>();
                                    if remaining_hap.len() == 3 {
                                        updated_alleles =
                                            vec![second_allele.to_vec(), remaining_hap];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if updated_alleles.is_empty() {
            updated_alleles = new_alleles.clone();
        }
        log::debug!(
            "RCCX allele update: input_alleles={:?}, output_alleles={:?}",
            new_alleles,
            updated_alleles
        );
        // check wrong phasing
        let mut wrong_allele = false;
        for allele in &updated_alleles {
            // if both copies are starting or ending
            if allele.len() == 2 {
                if let (Some(first_hap), Some(second_hap)) = (allele.first(), allele.last()) {
                    if (starting_copies.contains(first_hap) && starting_copies.contains(second_hap))
                        || (ending_copies.contains(first_hap) && ending_copies.contains(second_hap))
                    {
                        wrong_allele = true;
                    }
                }
            }
            for hap in final_haps {
                if ending_copies.contains(hap) && !two_cp_haplotypes.contains(hap) {
                    if allele.iter().filter(|x| *x == hap).count() > 1 {
                        wrong_allele = true;
                    }
                }
            }
        }
        if wrong_allele {
            log::debug!("RCCX allele phasing is ambiguous; clearing allele calls");
            updated_alleles = vec![];
        }

        if updated_alleles.len() == 2 {
            if let (Some(first_allele), Some(second_allele)) =
                (updated_alleles.first(), updated_alleles.last())
            {
                let a = final_haps
                    .iter()
                    .filter(|x| first_allele.contains(*x))
                    .count();
                let b = final_haps
                    .iter()
                    .filter(|x| second_allele.contains(*x))
                    .count();
                if a + b == nhap {
                    successful_phasing = true;
                }
            }
        }
        Ok((successful_phasing, updated_alleles, two_cp_haplotypes))
    }

    fn annotate_alleles(
        &mut self,
        successful_phasing: bool,
        new_alleles: &Vec<Vec<String>>,
        hap_variants: &BTreeMap<String, Vec<String>>,
        ending_copies: &Vec<String>,
        nhap: usize,
        two_cp_haplotypes: &Vec<String>,
    ) -> Result<Vec<Option<String>>, DError> {
        let mut annotated_alleles = Vec::new();
        //if new_alleles.is_empty() {
        //    return Ok(annotated_alleles);
        //}
        if successful_phasing {
            if let (Some(allele1), Some(allele2)) = (new_alleles.first(), new_alleles.last()) {
                let allele1_var = allele1
                    .iter()
                    .filter_map(|x| hap_variants.get(x).cloned())
                    .collect::<Vec<_>>();
                let allele2_var = allele2
                    .iter()
                    .filter_map(|x| hap_variants.get(x).cloned())
                    .collect::<Vec<_>>();
                let annotated_allele = self.annotate_var(&allele1_var)?;
                annotated_alleles.push(annotated_allele);
                let annotated_allele = self.annotate_var(&allele2_var)?;
                annotated_alleles.push(annotated_allele);
            }
        } else if ending_copies.len() == 2 && nhap == 4 && two_cp_haplotypes.is_empty() {
            for hap in ending_copies {
                let Some(allele_var) = hap_variants.get(hap) else {
                    continue;
                };
                if allele_var.is_empty() {
                    annotated_alleles.push(Some(String::from("WT")));
                } else {
                    annotated_alleles.push(Some(allele_var.join(",")));
                }
            }
        }
        Ok(annotated_alleles)
    }

    /// annotate an allele with variants
    fn annotate_var(&mut self, allele_var: &Vec<Vec<String>>) -> Result<Option<String>, DError> {
        let allele_len = allele_var.len();
        if allele_len == 2 {
            if allele_var.contains(&vec![]) {
                return Ok(Some(String::from("WT")));
            } else {
                let mut tmp = allele_var.clone();
                tmp.sort_by(|a, b| a.len().cmp(&(b.len())));
                if let Some(annotated_allele) = tmp.first() {
                    return Ok(Some(annotated_allele.join(",")));
                }
            }
        } else if allele_len == 1 {
            if allele_var.contains(&vec![]) {
                return Ok(Some(String::from("pseudogene_deletion")));
            } else {
                let mut annotated_allele = String::from("deletion_");
                if let Some(first) = allele_var.first() {
                    annotated_allele += first.join(",").as_str();
                }
                return Ok(Some(annotated_allele));
            }
        } else if allele_len == 3 {
            let mut tmp = allele_var.clone();
            tmp.sort_by(|a, b| a.len().cmp(&(b.len())));
            let first_hap = &tmp[0];
            let second_hap = &tmp[1];
            let third_hap = &tmp[2];
            if first_hap.is_empty() {
                if second_hap.is_empty() {
                    return Ok(Some(String::from("gene_duplication")));
                } else if third_hap.len() >= 6
                    && (second_hap.len() as i64 - third_hap.len() as i64).abs() <= 1
                {
                    return Ok(Some(String::from("pseudogene_duplication")));
                } else {
                    let mut annotated_allele = String::from("duplication_WT_plus_");
                    annotated_allele += second_hap.join(",").as_str();
                    return Ok(Some(annotated_allele));
                }
            } else if third_hap.len() >= 6
                && (second_hap.len() as i64 - third_hap.len() as i64).abs() <= 1
            {
                let mut annotated_allele = first_hap.join(",");
                annotated_allele += "_pseudogene_duplication";
                return Ok(Some(annotated_allele.to_string()));
            } else {
                let mut annotated_allele = String::from("duplication_");
                annotated_allele += first_hap.join(",").as_str();
                annotated_allele += "_plus_";
                annotated_allele += second_hap.join(",").as_str();
                return Ok(Some(annotated_allele));
            }
        }
        Ok(None)
    }

    /// Run RCCX-specific workflow, including CYP21/TNXB-aware hap naming and
    /// region-specific annotation synthesis.
    pub fn run_rccx(&mut self) -> Result<GeneCall, DError> {
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

        let _snp_file = self
            .locus_config()
            .get("data")
            .and_then(|x| x.as_mapping())
            .and_then(|x| x.get("snp_file"))
            .and_then(|x| x.as_str())
            .ok_or_else(|| {
                crate::phaser::Exception::new(format!(
                    "Missing data.snp_file in region config for gene '{}'",
                    self.gene_name()
                ))
            })?;

        const VARIANTDEF_38: &[u8] = std::include_bytes!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/data/38/rccx/cyp21_diff_sites.txt"
        ));
        const VARIANTDEF_19: &[u8] = std::include_bytes!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/data/19/rccx/cyp21_diff_sites.txt"
        ));
        let bed = if self.settings.genome == "19" || self.settings.genome == "37" {
            std::str::from_utf8(VARIANTDEF_19)?
        } else {
            std::str::from_utf8(VARIANTDEF_38)?
        };
        let bed = bed
            .split_terminator('\n')
            .map(std::borrow::ToOwned::to_owned)
            .collect::<Vec<_>>();
        let mut known_variants = BTreeMap::new();
        for line in bed {
            let fields = line
                .split_terminator(' ')
                .map(std::borrow::ToOwned::to_owned)
                .collect::<Vec<_>>();
            if fields.len() < 5 {
                log::warn!("Skipping malformed RCCX variant line: {line}");
                continue;
            }
            let var_pos = match fields[1].parse::<i64>() {
                Ok(v) => v - 1,
                Err(_) => {
                    log::warn!("Skipping RCCX variant line with invalid pos: {line}");
                    continue;
                }
            };
            let ref_base = fields[2].to_string();
            let alt_base = fields[3].to_string();
            let var_name = fields[4].to_string();
            let var_site = CandidateSite::new(var_pos, ref_base, alt_base);
            known_variants.insert(var_site, var_name);
        }
        log::debug!(
            "Loaded RCCX known variants for gene {}: {} entries",
            self.gene_name(),
            known_variants.len()
        );

        let (hom_sites_to_add, add_sites) = self.get_sites(&seq, None, None)?;
        // reverse del_data
        self.del_data.reverse();

        let tid = self.genome_tid().map(|x| x as i32).ok_or_else(|| {
            crate::phaser::Exception::new(format!(
                "Missing chromosome tid while running '{}' workflow",
                self.gene_name()
            ))
        })?;
        let del_reads = self.del_data[1].del_reads_partial.clone();
        let init_read_hap_map = self.haplotypes_from_reads(
            None,
            /* kept_sites */ &hom_sites_to_add,
            Some(&add_sites),
            /* partial_deletion_reads */ Some(&del_reads),
            (
                /* min_mapq= */ 5,
                /* check_clip= */ true,
                /* min_clip_len */ Some(50u32),
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
        // get haps that extend into tnxb
        let mut ending_copies = Vec::new();
        let mut starting_copies = Vec::new();
        let mut single_copies = Vec::new();
        for (idx, hap) in main_haps_clone.iter().enumerate() {
            let hap_name = format!("{mod_gene_name}_hap{}", idx + 1);
            assembled_haps.insert(hap.vstr(), hap_name.clone());
            if hap.len() >= 2 {
                let Some(first_base) = hap.first() else {
                    continue;
                };
                let Some(last_base) = hap.last() else {
                    continue;
                };
                if *first_base != b'x'
                    && *first_base != b'0'
                    && *last_base != b'x'
                    && *last_base != b'0'
                {
                    ending_copies.push(hap_name.clone());
                } else if *first_base == b'0' && *last_base == b'0' {
                    starting_copies.push(hap_name.clone());
                } else if *first_base == b'0' && *last_base != b'x' && *last_base != b'0' {
                    single_copies.push(hap_name.clone());
                }
            }
        }
        call.final_haplotypes = assembled_haps
            .clone()
            .into_iter()
            .map(|(k, v)| (k.to_string(), v))
            .collect::<BTreeMap<_, _>>();
        // Phase alleles
        let allele_result = self.phase_alleles(&mut phase_results, &assembled_haps, None);
        call.region_specific_info.insert(
            String::from("haplotype_links"),
            serde_json::to_value(&allele_result.haplotype_links)?,
        );
        let alleles = allele_result.raw_alleles;

        // Output variants
        let haps =
            self.output_variants_in_haps(&phase_results, &known_del, assembled_haps.clone())?;
        call.haplotype_details = haps
            .iter()
            .map(|(key, val)| (key.clone(), HapInfoForJson::from(val)))
            .collect::<BTreeMap<_, _>>();

        // update alleles
        let hcn = phase_results.assemblies.highest_cn;
        let (successful_phasing, updated_alleles, two_cp_haplotypes) = self.update_alleles(
            &alleles,
            &haps,
            &assembled_haps,
            &single_copies,
            &starting_copies,
            &ending_copies,
            hcn,
        )?;

        // annotate haplotypes by checking the diff sites
        // output variants carried by each haplotype
        let mut hap_variants: BTreeMap<String, Vec<String>> = BTreeMap::new();
        for (hap, hap_info) in &haps {
            hap_variants.entry(hap.to_string()).or_default();
            let variants = &hap_info.variants;
            for var in variants {
                if known_variants.contains_key(var) {
                    if let Some(var_name) = known_variants.get(var) {
                        hap_variants
                            .entry(hap.to_string())
                            .or_default()
                            .push(var_name.to_string());
                    }
                }
            }
        }
        let total_cn = assembled_haps.len() + two_cp_haplotypes.len();
        call.total_cn = Some(total_cn as i32);
        if total_cn < 2 || ending_copies.len() > 2 {
            call.total_cn = None;
        }
        let annotated_alleles = self.annotate_alleles(
            successful_phasing,
            &updated_alleles,
            &hap_variants,
            &ending_copies,
            assembled_haps.len(),
            &two_cp_haplotypes,
        )?;

        call.two_copy_haplotypes = two_cp_haplotypes;
        call.region_specific_info.insert(
            String::from("alleles_final"),
            updated_alleles.clone().into(),
        );
        call.region_specific_info
            .insert(String::from("raw_alleles"), updated_alleles.clone().into());
        self.fill_in_call(phase_results, &mut call);
        call.region_specific_info
            .insert(String::from("phasing_success"), successful_phasing.into());
        call.region_specific_info
            .insert(String::from("starting_hap"), starting_copies.into());
        call.region_specific_info
            .insert(String::from("ending_hap"), ending_copies.into());
        call.region_specific_info
            .insert(String::from("deletion_hap"), single_copies.into());
        call.region_specific_info.insert(
            String::from("hap_variants"),
            serde_json::to_value(&hap_variants)?,
        );
        call.region_specific_info.insert(
            String::from("annotated_alleles"),
            serde_json::to_value(&annotated_alleles)?,
        );
        Ok(call)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config;
    use crate::phaser;
    use crate::toolkit::util;

    fn build_test_phaser(outdir: &std::path::Path) -> Phaser {
        let settings = phaser::Settings::new(
            "TEST",
            (
                util::test_file("ref/smn1_ref.fa"),
                util::test_file("HG00733_smn1_realigned.bam"),
            ),
            outdir,
            "rccx",
            &config::Region::try_load(None).expect("region config should load"),
            None,
            None,
            String::from("38"),
            None,
            0.03,
            false,
        );
        let gene_config = config::Gene::try_load(None).expect("gene config should load");
        Phaser::new(settings, Some(gene_config), None, None).expect("phaser should build")
    }

    #[test]
    fn annotate_var_matches_python_cases() {
        let outdir = tempfile::TempDir::new().expect("tempdir should build");
        let mut phaser = build_test_phaser(outdir.path());

        assert_eq!(
            phaser.annotate_var(&vec![vec![], vec![]]).unwrap(),
            Some(String::from("WT"))
        );
        assert_eq!(
            phaser
                .annotate_var(&vec![vec![], vec![String::from("var1")]])
                .unwrap(),
            Some(String::from("WT"))
        );
        assert_eq!(
            phaser
                .annotate_var(&vec![
                    vec![String::from("var1")],
                    vec![String::from("var1"), String::from("var2")],
                ])
                .unwrap(),
            Some(String::from("var1"))
        );
        assert_eq!(
            phaser
                .annotate_var(&vec![
                    vec![
                        String::from("var1"),
                        String::from("var2"),
                        String::from("var3"),
                    ],
                    vec![String::from("var1"), String::from("var2")],
                ])
                .unwrap(),
            Some(String::from("var1,var2"))
        );
        assert_eq!(
            phaser.annotate_var(&vec![vec![]]).unwrap(),
            Some(String::from("pseudogene_deletion"))
        );
        assert_eq!(
            phaser
                .annotate_var(&vec![vec![String::from("var1")]])
                .unwrap(),
            Some(String::from("deletion_var1"))
        );
        assert_eq!(
            phaser
                .annotate_var(&vec![
                    vec![String::from("var1"), String::from("var2")],
                    vec![],
                    vec![],
                ])
                .unwrap(),
            Some(String::from("gene_duplication"))
        );
        assert_eq!(
            phaser
                .annotate_var(&vec![
                    vec![
                        String::from("var1"),
                        String::from("var2"),
                        String::from("var3"),
                        String::from("var4"),
                        String::from("var5"),
                        String::from("var6"),
                    ],
                    vec![
                        String::from("var1"),
                        String::from("var2"),
                        String::from("var3"),
                        String::from("var4"),
                        String::from("var5"),
                        String::from("var6"),
                    ],
                    vec![],
                ])
                .unwrap(),
            Some(String::from("pseudogene_duplication"))
        );
        assert_eq!(
            phaser
                .annotate_var(&vec![
                    vec![String::from("var1"), String::from("var2")],
                    vec![String::from("var1")],
                    vec![],
                ])
                .unwrap(),
            Some(String::from("duplication_WT_plus_var1"))
        );
        assert_eq!(
            phaser
                .annotate_var(&vec![
                    vec![
                        String::from("var1"),
                        String::from("var2"),
                        String::from("var3"),
                        String::from("var4"),
                        String::from("var5"),
                        String::from("var6"),
                    ],
                    vec![
                        String::from("var1"),
                        String::from("var2"),
                        String::from("var3"),
                        String::from("var4"),
                        String::from("var5"),
                        String::from("var6"),
                    ],
                    vec![String::from("var1"), String::from("var2")],
                ])
                .unwrap(),
            Some(String::from("var1,var2_pseudogene_duplication"))
        );
        assert_eq!(
            phaser
                .annotate_var(&vec![
                    vec![
                        String::from("var1"),
                        String::from("var2"),
                        String::from("var3"),
                        String::from("var4"),
                        String::from("var5"),
                        String::from("var6"),
                    ],
                    vec![
                        String::from("var1"),
                        String::from("var2"),
                        String::from("var3"),
                        String::from("var4"),
                    ],
                    vec![String::from("var1"), String::from("var2")],
                ])
                .unwrap(),
            Some(String::from(
                "duplication_var1,var2_plus_var1,var2,var3,var4"
            ))
        );
        assert_eq!(
            phaser
                .annotate_var(&vec![
                    vec![],
                    vec![],
                    vec![],
                    vec![String::from("var1"), String::from("var2")],
                ])
                .unwrap(),
            None
        );
    }

    #[test]
    fn annotate_alleles_matches_python_cases() {
        let outdir = tempfile::TempDir::new().expect("tempdir should build");
        let mut phaser = build_test_phaser(outdir.path());
        let alleles = vec![
            vec![String::from("hap1"), String::from("hap2")],
            vec![String::from("hap3"), String::from("hap4")],
        ];
        let hap_variants = BTreeMap::from([
            (String::from("hap1"), vec![]),
            (String::from("hap2"), vec![]),
            (String::from("hap3"), vec![]),
            (String::from("hap4"), vec![]),
        ]);
        let empty: Vec<String> = vec![];

        assert_eq!(
            phaser
                .annotate_alleles(true, &alleles, &hap_variants, &empty, 4, &vec![])
                .unwrap(),
            vec![Some(String::from("WT")), Some(String::from("WT"))]
        );

        assert_eq!(
            phaser
                .annotate_alleles(false, &alleles, &hap_variants, &empty, 4, &vec![])
                .unwrap(),
            Vec::<Option<String>>::new()
        );

        let ending = vec![String::from("hap1"), String::from("hap2")];
        assert_eq!(
            phaser
                .annotate_alleles(false, &alleles, &hap_variants, &ending, 4, &vec![])
                .unwrap(),
            vec![Some(String::from("WT")), Some(String::from("WT"))]
        );

        let hap_variants = BTreeMap::from([
            (String::from("hap1"), vec![String::from("var1")]),
            (String::from("hap2"), vec![]),
            (String::from("hap3"), vec![]),
            (String::from("hap4"), vec![]),
        ]);
        assert_eq!(
            phaser
                .annotate_alleles(false, &alleles, &hap_variants, &ending, 4, &vec![])
                .unwrap(),
            vec![Some(String::from("var1")), Some(String::from("WT"))]
        );
    }
}
