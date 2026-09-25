// RCCX specific caller
use crate::io::json::GeneCall;
use crate::phaser::Phaser;
use crate::phaser::{HapInfo, HapInfoForJson};
use crate::toolkit::site_selection::CandidateSite;
use crate::toolkit::util::DError;
use itertools::intersperse;
use std::collections::{BTreeMap, BTreeSet};
use vstr::VStr;

/// Build an RCCX haplotype rename map from the finalized allele order.
///
/// Renumbering only applies when successful phasing produced exactly two alleles
/// whose unique haplotypes fully cover the assembled haplotypes.
fn rccx_haplotype_rename_map(
    successful_phasing: bool,
    updated_alleles: &[Vec<String>],
    assembled_haps: &BTreeMap<VStr<'_>, String>,
    two_cp_haplotypes: &[String],
    gene_name: &str,
) -> Option<BTreeMap<String, String>> {
    if !successful_phasing || updated_alleles.len() != 2 {
        return None;
    }

    let ordered_haps = updated_alleles
        .iter()
        .flatten()
        .cloned()
        .collect::<Vec<_>>();
    if ordered_haps.len() != assembled_haps.len() + two_cp_haplotypes.len() {
        return None;
    }

    let assembled_names = assembled_haps.values().collect::<BTreeSet<_>>();
    let ordered_names = ordered_haps.iter().collect::<BTreeSet<_>>();
    if ordered_names != assembled_names {
        return None;
    }

    let mut rename_map = BTreeMap::new();
    let mut index = 0;
    for hap in &ordered_haps {
        if !rename_map.contains_key(hap) {
            index += 1;
            rename_map.insert(hap.to_string(), format!("{gene_name}_hap{}", index));
        }
    }
    Some(rename_map)
}

/// Reorder an RCCX allele by its configured start and end copy assignments.
///
/// Haplotype order is preserved within each group. A haplotype present in both
/// boundary lists is treated as a starting copy to avoid emitting it twice.
fn reorder_rccx_allele(
    allele: &[String],
    starting_copies: &[String],
    ending_copies: &[String],
) -> Vec<String> {
    let reordered_allele = allele
        .iter()
        .filter(|hap| starting_copies.contains(*hap))
        .chain(
            allele
                .iter()
                .filter(|hap| !starting_copies.contains(*hap) && !ending_copies.contains(*hap)),
        )
        .chain(
            allele
                .iter()
                .filter(|hap| !starting_copies.contains(*hap) && ending_copies.contains(*hap)),
        )
        .cloned()
        .collect::<Vec<_>>();
    reordered_allele
}

/// Clear RCCX allele calls when boundary copies indicate ambiguous phasing.
fn check_wrong_allele(
    updated_alleles: &mut Vec<Vec<String>>,
    starting_copies: &[String],
    ending_copies: &[String],
    single_copies: &[String],
) {
    let mut wrong_allele = false;
    if updated_alleles.len() != 2 {
        wrong_allele = true;
    }
    for allele in updated_alleles.iter() {
        let starting_hap_count = allele
            .iter()
            .filter(|hap| starting_copies.contains(*hap))
            .count();
        let ending_hap_count = allele
            .iter()
            .filter(|hap| ending_copies.contains(*hap))
            .count();
        let single_copy_count = allele
            .iter()
            .filter(|hap| single_copies.contains(*hap))
            .count();
        if starting_hap_count > 1 || ending_hap_count > 1 {
            wrong_allele = true;
        }
        if single_copy_count > 0 && allele.len() > 1 {
            wrong_allele = true;
        }
        if allele.is_empty() {
            wrong_allele = true;
        }
    }
    if wrong_allele {
        log::debug!("RCCX allele phasing is ambiguous; clearing allele calls");
        updated_alleles.clear();
    }
}

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
        haplotype_links: &BTreeMap<String, Vec<String>>,
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
                        updated_alleles = vec![first_allele.clone(), single_copies.clone()];
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
            }
        } else if single_copies.is_empty() {
            // homozygous, each haplotype has cn 2
            if nhap == 2
                && ending_copies.len() == 1
                && starting_copies.len() == 1
                && new_alleles.len() == 1
            {
                two_cp_haplotypes = final_haps.clone();
                if let Some(first_allele) = new_alleles.first() {
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
                            }
                        } else if nhap == 4 {
                            let middle_copies = final_haps
                                .iter()
                                .filter(|x| {
                                    !ending_copies.contains(x) && !starting_copies.contains(x)
                                })
                                .collect::<Vec<_>>();
                            if middle_copies.len() == 1 {
                                let middle_copy = middle_copies.first().unwrap().to_string();
                                if let Some(links) = haplotype_links.get(&middle_copy) {
                                    let middle_copy_linked_starting_copies = links
                                        .iter()
                                        .filter(|x| starting_copies.contains(x))
                                        .collect::<Vec<_>>();
                                    if middle_copy_linked_starting_copies.len() == 1 {
                                        let middle_copy_linked_starting_copy =
                                            middle_copy_linked_starting_copies
                                                .first()
                                                .unwrap()
                                                .to_string();
                                        let other_starting_copy = starting_copies
                                            .iter()
                                            .filter(|x| {
                                                !middle_copy_linked_starting_copies.contains(x)
                                            })
                                            .next()
                                            .unwrap()
                                            .to_string();
                                        updated_alleles = vec![
                                            vec![
                                                middle_copy_linked_starting_copy,
                                                middle_copy,
                                                ending_copy.to_string(),
                                            ],
                                            vec![other_starting_copy, ending_copy.to_string()],
                                        ];
                                    }
                                }
                            }
                        } else if nhap == 5 {
                            let middle_copies = final_haps
                                .iter()
                                .filter(|x| {
                                    !ending_copies.contains(x) && !starting_copies.contains(x)
                                })
                                .collect::<Vec<_>>();
                            if middle_copies.len() == 2 {
                                let first_middle_copy = middle_copies.first().unwrap().to_string();
                                let second_middle_copy = middle_copies.last().unwrap().to_string();
                                if let (Some(links1), Some(links2)) = (
                                    haplotype_links.get(&first_middle_copy),
                                    haplotype_links.get(&second_middle_copy),
                                ) {
                                    let first_middle_copy_linked_starting_copies = links1
                                        .iter()
                                        .filter(|x| starting_copies.contains(x))
                                        .collect::<Vec<_>>();
                                    let second_middle_copy_linked_starting_copies = links2
                                        .iter()
                                        .filter(|x| starting_copies.contains(x))
                                        .collect::<Vec<_>>();
                                    if first_middle_copy_linked_starting_copies.len() == 1
                                        && second_middle_copy_linked_starting_copies.len() == 1
                                    {
                                        let first_middle_copy_linked_starting_copy =
                                            first_middle_copy_linked_starting_copies
                                                .first()
                                                .unwrap()
                                                .to_string();
                                        let second_middle_copy_linked_starting_copy =
                                            second_middle_copy_linked_starting_copies
                                                .first()
                                                .unwrap()
                                                .to_string();
                                        updated_alleles = vec![
                                            vec![
                                                first_middle_copy_linked_starting_copy,
                                                first_middle_copy,
                                                ending_copy.to_string(),
                                            ],
                                            vec![
                                                second_middle_copy_linked_starting_copy,
                                                second_middle_copy,
                                                ending_copy.to_string(),
                                            ],
                                        ];
                                    }
                                }
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
                            }
                        } else if nhap == 4 {
                            let middle_copies = final_haps
                                .iter()
                                .filter(|x| {
                                    !ending_copies.contains(x) && !starting_copies.contains(x)
                                })
                                .collect::<Vec<_>>();
                            if middle_copies.len() == 1 {
                                let middle_copy = middle_copies.first().unwrap().to_string();
                                if let Some(links) = haplotype_links.get(&middle_copy) {
                                    let middle_copy_linked_ending_copies = links
                                        .iter()
                                        .filter(|x| ending_copies.contains(x))
                                        .collect::<Vec<_>>();
                                    if middle_copy_linked_ending_copies.len() == 1 {
                                        let middle_copy_linked_ending_copy =
                                            middle_copy_linked_ending_copies
                                                .first()
                                                .unwrap()
                                                .to_string();
                                        let other_ending_copy = ending_copies
                                            .iter()
                                            .filter(|x| {
                                                !middle_copy_linked_ending_copies.contains(x)
                                            })
                                            .next()
                                            .unwrap()
                                            .to_string();
                                        updated_alleles = vec![
                                            vec![
                                                starting_copy.to_string(),
                                                middle_copy,
                                                middle_copy_linked_ending_copy,
                                            ],
                                            vec![starting_copy.to_string(), other_ending_copy],
                                        ];
                                    }
                                }
                            }
                        } else if nhap == 5 {
                            let middle_copies = final_haps
                                .iter()
                                .filter(|x| {
                                    !ending_copies.contains(x) && !starting_copies.contains(x)
                                })
                                .collect::<Vec<_>>();
                            if middle_copies.len() == 2 {
                                let first_middle_copy = middle_copies.first().unwrap().to_string();
                                let second_middle_copy = middle_copies.last().unwrap().to_string();
                                if let (Some(links1), Some(links2)) = (
                                    haplotype_links.get(&first_middle_copy),
                                    haplotype_links.get(&second_middle_copy),
                                ) {
                                    let first_middle_copy_linked_ending_copies = links1
                                        .iter()
                                        .filter(|x| ending_copies.contains(x))
                                        .collect::<Vec<_>>();
                                    let second_middle_copy_linked_ending_copies = links2
                                        .iter()
                                        .filter(|x| ending_copies.contains(x))
                                        .collect::<Vec<_>>();
                                    if first_middle_copy_linked_ending_copies.len() == 1
                                        && second_middle_copy_linked_ending_copies.len() == 1
                                    {
                                        let first_middle_copy_linked_ending_copy =
                                            first_middle_copy_linked_ending_copies
                                                .first()
                                                .unwrap()
                                                .to_string();
                                        let second_middle_copy_linked_ending_copy =
                                            second_middle_copy_linked_ending_copies
                                                .first()
                                                .unwrap()
                                                .to_string();
                                        updated_alleles = vec![
                                            vec![
                                                starting_copy.to_string(),
                                                first_middle_copy,
                                                first_middle_copy_linked_ending_copy.to_string(),
                                            ],
                                            vec![
                                                starting_copy.to_string(),
                                                second_middle_copy,
                                                second_middle_copy_linked_ending_copy.to_string(),
                                            ],
                                        ];
                                    }
                                }
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
                                updated_alleles = vec![first_allele.clone(), remaining_hap];
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
                                                vec![first_allele.clone(), remaining_hap];
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
                                        updated_alleles = vec![first_allele.clone(), remaining_hap];
                                    }
                                } else if allele2 && !allele1 {
                                    let remaining_hap = final_haps
                                        .iter()
                                        .filter(|x| !second_allele.contains(*x))
                                        .map(|x| x.to_string())
                                        .collect::<Vec<_>>();
                                    if remaining_hap.len() == 3 {
                                        updated_alleles =
                                            vec![second_allele.clone(), remaining_hap];
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
        check_wrong_allele(
            &mut updated_alleles,
            starting_copies,
            ending_copies,
            single_copies,
        );

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
                if a + b == nhap + two_cp_haplotypes.len() {
                    let reordered_first_allele =
                        reorder_rccx_allele(first_allele, starting_copies, ending_copies);
                    log::debug!("reordered_first_allele={:?}", reordered_first_allele);
                    let reordered_second_allele =
                        reorder_rccx_allele(second_allele, starting_copies, ending_copies);
                    log::debug!("reordered_second_allele={:?}", reordered_second_allele);
                    updated_alleles = vec![reordered_first_allele, reordered_second_allele];
                    successful_phasing = true;
                }
            }
        }
        // if phasing is not successful, clear the alleles
        if !successful_phasing {
            updated_alleles = vec![];
        }
        Ok((successful_phasing, updated_alleles, two_cp_haplotypes))
    }

    fn annotate_alleles(
        &mut self,
        successful_phasing: bool,
        new_alleles: &Vec<Vec<String>>,
        hap_variants: &BTreeMap<String, Vec<String>>,
        ending_copies: &Vec<String>,
        single_copies: &Vec<String>,
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
        } else if ending_copies.len() == 2
            && nhap == 4
            && two_cp_haplotypes.is_empty()
            && single_copies.is_empty()
        {
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
        // Phase alleles
        let allele_result = self.phase_alleles(&mut phase_results, &assembled_haps, None);
        let mut haplotype_links = allele_result.haplotype_links;
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
        let (successful_phasing, mut updated_alleles, mut two_cp_haplotypes) = self
            .update_alleles(
                &alleles,
                &haps,
                &assembled_haps,
                &single_copies,
                &starting_copies,
                &ending_copies,
                &haplotype_links,
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
            &single_copies,
            assembled_haps.len(),
            &two_cp_haplotypes,
        )?;

        if let Some(renamed_haps) = rccx_haplotype_rename_map(
            successful_phasing,
            &updated_alleles,
            &assembled_haps,
            &two_cp_haplotypes,
            &mod_gene_name,
        ) {
            let rename = |name: &str| {
                renamed_haps
                    .get(name)
                    .cloned()
                    .expect("RCCX haplotype rename map should cover all assembled haplotypes")
            };
            for hap_name in assembled_haps.values_mut() {
                *hap_name = rename(hap_name);
            }
            for allele in &mut updated_alleles {
                for hap_name in allele {
                    *hap_name = rename(hap_name);
                }
            }
            for hap_name in &mut two_cp_haplotypes {
                *hap_name = rename(hap_name);
            }
            for hap_name in &mut starting_copies {
                *hap_name = rename(hap_name);
            }
            for hap_name in &mut ending_copies {
                *hap_name = rename(hap_name);
            }
            for hap_name in &mut single_copies {
                *hap_name = rename(hap_name);
            }

            hap_variants = hap_variants
                .into_iter()
                .map(|(hap_name, variants)| (rename(&hap_name), variants))
                .collect();
            haplotype_links = haplotype_links
                .into_iter()
                .map(|(hap_name, links)| {
                    (
                        rename(&hap_name),
                        links.into_iter().map(|link| rename(&link)).collect(),
                    )
                })
                .collect();
            call.haplotype_details = std::mem::take(&mut call.haplotype_details)
                .into_iter()
                .map(|(hap_name, hap_info)| (rename(&hap_name), hap_info))
                .collect();
        }

        // Report assembled haplotypes even when allele phasing could not rename them.
        call.final_haplotypes = assembled_haps
            .into_iter()
            .map(|(hap_sequence, hap_name)| (hap_sequence.to_string(), hap_name))
            .collect();
        call.two_copy_haplotypes = two_cp_haplotypes;
        call.region_specific_info.insert(
            String::from("alleles_final"),
            updated_alleles.clone().into(),
        );
        call.region_specific_info
            .insert(String::from("raw_alleles"), updated_alleles.clone().into());
        call.region_specific_info.insert(
            String::from("haplotype_links"),
            serde_json::to_value(&haplotype_links)?,
        );
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

    #[test]
    fn check_wrong_allele_clears_all_calls_for_multiple_boundary_copies() {
        let starting_copies = vec![String::from("start1"), String::from("start2")];
        let ending_copies = vec![String::from("end1"), String::from("end2")];
        for invalid_allele in [
            vec!["start1", "start2"],
            vec!["end1", "end2"],
            vec!["start1", "middle", "start2", "end1"],
            vec!["start1", "end1", "middle", "end2"],
            vec!["start1", "start1", "end1"],
            vec!["start1", "end1", "end1"],
        ] {
            let mut alleles = vec![
                vec![String::from("start1"), String::from("end1")],
                invalid_allele.iter().map(|hap| hap.to_string()).collect(),
            ];
            check_wrong_allele(&mut alleles, &starting_copies, &ending_copies, &[]);
            assert!(alleles.is_empty(), "invalid allele: {invalid_allele:?}");
        }
    }

    #[test]
    fn check_wrong_allele_rejects_single_copy_deletions_mixed_with_other_haplotypes() {
        let single_copies = vec![String::from("deletion1"), String::from("deletion2")];
        for invalid_allele in [
            vec!["deletion1", "middle"],
            vec!["start", "deletion1", "end"],
            vec!["deletion1", "deletion2"],
            vec!["deletion1", "deletion1"],
        ] {
            let mut alleles = vec![
                invalid_allele.iter().map(|hap| hap.to_string()).collect(),
                vec![String::from("deletion2")],
            ];
            check_wrong_allele(
                &mut alleles,
                &[String::from("start")],
                &[String::from("end")],
                &single_copies,
            );
            assert!(alleles.is_empty(), "invalid allele: {invalid_allele:?}");
        }
    }

    #[test]
    fn check_wrong_allele_preserves_two_valid_nonempty_alleles() {
        let starting_copies = vec![String::from("start")];
        let ending_copies = vec![String::from("end")];
        let single_copies = vec![String::from("deletion")];
        for input in [
            vec![vec!["start", "middle", "end"], vec!["deletion"]],
            vec![vec!["start", "end"], vec!["start", "end"]],
            vec![vec!["deletion"], vec!["deletion"]],
            vec![vec!["middle", "middle"], vec!["deletion"]],
            vec![vec!["start", "middle"], vec!["middle", "end"]],
            vec![vec!["middle"], vec!["middle"]],
        ] {
            let expected = input
                .iter()
                .map(|allele| allele.iter().map(|hap| hap.to_string()).collect::<Vec<_>>())
                .collect::<Vec<_>>();
            let mut alleles = expected.clone();
            check_wrong_allele(
                &mut alleles,
                &starting_copies,
                &ending_copies,
                &single_copies,
            );
            assert_eq!(alleles, expected);
        }
    }

    #[test]
    fn check_wrong_allele_requires_exactly_two_alleles() {
        for allele_count in [0, 1, 3, 4] {
            let mut alleles = vec![vec![String::from("deletion")]; allele_count];
            check_wrong_allele(&mut alleles, &[], &[], &[String::from("deletion")]);
            assert!(alleles.is_empty(), "allele count: {allele_count}");
        }
    }

    #[test]
    fn check_wrong_allele_clears_calls_when_either_allele_is_empty() {
        for mut alleles in [
            vec![vec![], vec![String::from("deletion")]],
            vec![vec![String::from("deletion")], vec![]],
            vec![vec![], vec![]],
        ] {
            check_wrong_allele(&mut alleles, &[], &[], &[String::from("deletion")]);
            assert!(alleles.is_empty());
        }
    }

    fn assembled_haps_for_renaming() -> BTreeMap<VStr<'static>, String> {
        BTreeMap::from([
            (VStr::from("1111"), String::from("rccx_hap1")),
            (VStr::from("1112"), String::from("rccx_hap2")),
            (VStr::from("1121"), String::from("rccx_hap3")),
            (VStr::from("1122"), String::from("rccx_hap4")),
        ])
    }

    #[test]
    fn rccx_haplotype_rename_map_follows_updated_allele_order() {
        let updated_alleles = vec![
            vec![String::from("rccx_hap2"), String::from("rccx_hap3")],
            vec![String::from("rccx_hap1"), String::from("rccx_hap4")],
        ];

        assert_eq!(
            rccx_haplotype_rename_map(
                true,
                &updated_alleles,
                &assembled_haps_for_renaming(),
                &[],
                "rccx",
            ),
            Some(BTreeMap::from([
                (String::from("rccx_hap1"), String::from("rccx_hap3")),
                (String::from("rccx_hap2"), String::from("rccx_hap1")),
                (String::from("rccx_hap3"), String::from("rccx_hap2")),
                (String::from("rccx_hap4"), String::from("rccx_hap4")),
            ]))
        );
    }

    #[test]
    fn rccx_haplotype_rename_map_requires_complete_successful_diploid_phasing() {
        let assembled_haps = assembled_haps_for_renaming();
        let complete_alleles = vec![
            vec![String::from("rccx_hap1"), String::from("rccx_hap2")],
            vec![String::from("rccx_hap3"), String::from("rccx_hap4")],
        ];
        assert!(
            rccx_haplotype_rename_map(false, &complete_alleles, &assembled_haps, &[], "rccx")
                .is_none()
        );
        assert!(rccx_haplotype_rename_map(
            true,
            &complete_alleles[..1],
            &assembled_haps,
            &[],
            "rccx",
        )
        .is_none());
        assert!(rccx_haplotype_rename_map(
            true,
            &vec![
                vec![String::from("rccx_hap1")],
                vec![String::from("rccx_hap2")],
            ],
            &assembled_haps,
            &[],
            "rccx",
        )
        .is_none());
    }

    #[test]
    fn rccx_haplotype_rename_map_numbers_shared_haplotypes_only_once() {
        let assembled_haps = assembled_haps_for_renaming();
        let updated_alleles = vec![
            vec![String::from("rccx_hap3"), String::from("rccx_hap2")],
            vec![
                String::from("rccx_hap1"),
                String::from("rccx_hap4"),
                String::from("rccx_hap2"),
            ],
        ];
        let two_cp_haplotypes = vec![String::from("rccx_hap2")];

        assert_eq!(
            rccx_haplotype_rename_map(
                true,
                &updated_alleles,
                &assembled_haps,
                &two_cp_haplotypes,
                "rccx",
            ),
            Some(BTreeMap::from([
                (String::from("rccx_hap1"), String::from("rccx_hap3")),
                (String::from("rccx_hap2"), String::from("rccx_hap2")),
                (String::from("rccx_hap3"), String::from("rccx_hap1")),
                (String::from("rccx_hap4"), String::from("rccx_hap4")),
            ]))
        );
        assert!(
            rccx_haplotype_rename_map(true, &updated_alleles, &assembled_haps, &[], "rccx",)
                .is_none()
        );
    }

    #[test]
    fn reorder_rccx_allele_places_boundary_copies_at_their_respective_ends() {
        let allele = vec![
            String::from("rccx_hap3"),
            String::from("rccx_hap1"),
            String::from("rccx_hap4"),
            String::from("rccx_hap2"),
        ];
        let starting_copies = vec![String::from("rccx_hap2")];
        let ending_copies = vec![String::from("rccx_hap4")];

        assert_eq!(
            reorder_rccx_allele(&allele, &starting_copies, &ending_copies),
            vec![
                String::from("rccx_hap2"),
                String::from("rccx_hap3"),
                String::from("rccx_hap1"),
                String::from("rccx_hap4"),
            ]
        );
    }

    #[test]
    fn reorder_rccx_allele_allows_missing_or_multiple_boundary_copies() {
        let allele = vec![
            String::from("hap3"),
            String::from("hap1"),
            String::from("hap2"),
        ];
        for (starting_copies, ending_copies, expected) in [
            (vec![], vec![], vec!["hap3", "hap1", "hap2"]),
            (
                vec![],
                vec![String::from("hap3")],
                vec!["hap1", "hap2", "hap3"],
            ),
            (
                vec![String::from("hap1")],
                vec![],
                vec!["hap1", "hap3", "hap2"],
            ),
            (
                vec![String::from("hap2"), String::from("hap1")],
                vec![String::from("hap3")],
                vec!["hap1", "hap2", "hap3"],
            ),
            (
                vec![String::from("hap1")],
                vec![String::from("hap2"), String::from("hap3")],
                vec!["hap1", "hap3", "hap2"],
            ),
        ] {
            assert_eq!(
                reorder_rccx_allele(&allele, &starting_copies, &ending_copies),
                expected
            );
        }
    }

    #[test]
    fn reorder_rccx_allele_preserves_alleles_without_boundary_copies() {
        let allele = vec![String::from("hap1")];
        assert_eq!(reorder_rccx_allele(&allele, &[], &[]), allele);
        assert_eq!(
            reorder_rccx_allele(&[String::from("hap1"), String::from("hap2")], &[], &[],),
            vec![String::from("hap1"), String::from("hap2")]
        );
    }

    #[test]
    fn reorder_rccx_allele_emits_overlapping_boundary_copies_only_once() {
        let allele = vec![
            String::from("hap1"),
            String::from("hap2"),
            String::from("hap3"),
        ];
        let starting_copies = vec![String::from("hap2")];
        let ending_copies = vec![String::from("hap2"), String::from("hap1")];
        assert_eq!(
            reorder_rccx_allele(&allele, &starting_copies, &ending_copies),
            vec![
                String::from("hap2"),
                String::from("hap3"),
                String::from("hap1")
            ]
        );
        assert!(reorder_rccx_allele(&[], &starting_copies, &ending_copies).is_empty());
    }

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
                .annotate_alleles(true, &alleles, &hap_variants, &empty, &empty, 4, &vec![])
                .unwrap(),
            vec![Some(String::from("WT")), Some(String::from("WT"))]
        );

        assert_eq!(
            phaser
                .annotate_alleles(false, &alleles, &hap_variants, &empty, &empty, 4, &vec![])
                .unwrap(),
            Vec::<Option<String>>::new()
        );

        let ending = vec![String::from("hap1"), String::from("hap2")];
        assert_eq!(
            phaser
                .annotate_alleles(false, &alleles, &hap_variants, &ending, &empty, 4, &vec![])
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
                .annotate_alleles(false, &alleles, &hap_variants, &ending, &empty, 4, &vec![])
                .unwrap(),
            vec![Some(String::from("var1")), Some(String::from("WT"))]
        );
    }

    #[test]
    fn annotate_alleles_requires_successful_phasing_when_single_copies_are_present() {
        let outdir = tempfile::TempDir::new().expect("tempdir should build");
        let mut phaser = build_test_phaser(outdir.path());
        let alleles = vec![
            vec![String::from("deletion")],
            vec![
                String::from("start"),
                String::from("middle"),
                String::from("end"),
            ],
        ];
        let hap_variants = BTreeMap::from([
            (String::from("deletion"), vec![String::from("var1")]),
            (String::from("start"), vec![]),
            (String::from("middle"), vec![]),
            (String::from("end"), vec![]),
        ]);
        // Satisfy every other fallback condition to isolate the single-copy guard.
        let ending_copies = vec![String::from("middle"), String::from("end")];
        let single_copies = vec![String::from("deletion")];
        assert_eq!(
            phaser
                .annotate_alleles(
                    false,
                    &vec![],
                    &hap_variants,
                    &ending_copies,
                    &vec![],
                    4,
                    &vec![]
                )
                .unwrap(),
            vec![Some(String::from("WT")), Some(String::from("WT"))]
        );
        assert!(phaser
            .annotate_alleles(
                false,
                &vec![],
                &hap_variants,
                &ending_copies,
                &single_copies,
                4,
                &vec![]
            )
            .unwrap()
            .is_empty());
        assert_eq!(
            phaser
                .annotate_alleles(
                    true,
                    &alleles,
                    &hap_variants,
                    &ending_copies,
                    &single_copies,
                    4,
                    &vec![]
                )
                .unwrap(),
            vec![
                Some(String::from("deletion_var1")),
                Some(String::from("gene_duplication")),
            ]
        );
    }
}
