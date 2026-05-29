// OPN1LW specific caller
use crate::depth::Sex;
use crate::io::json::GeneCall;
use crate::phaser::HapInfoForJson;
use crate::phaser::Phaser;
use crate::toolkit::site_selection::CandidateSite;
use crate::toolkit::util::DError;
use itertools::{intersperse, Itertools};
use std::collections::BTreeMap;
use std::str::FromStr;

impl Phaser {
    /// Call known variants in Exon 3
    fn call_exon3(
        &mut self,
        hap_vars: &Vec<CandidateSite>,
        //exon3_vars: BTreeMap<(String, String), Vec<CandidateSite>>,
    ) -> Result<String, DError> {
        // let mut exon3_vars = BTreeMap::new();
        let mut annotated_vars = Vec::new();
        let exon3_variants = self
            .locus_config()
            .get("exon3_variants")
            .and_then(|x| x.as_sequence())
            .map(|seq| seq.iter().collect::<Vec<_>>())
            .unwrap_or_default();
        log::debug!(
            "Loaded {} configured OPN1LW exon3 variant groups",
            exon3_variants.len()
        );
        for sites in exon3_variants {
            let Some(sites_vec) = sites.as_sequence() else {
                continue;
            };
            if sites_vec.len() < 3 {
                continue;
            }
            let variants = sites_vec
                .first()
                .and_then(|v| v.as_sequence())
                .map(|x| {
                    x.iter()
                        .filter_map(|x| x.as_str().map(CandidateSite::from_str))
                        .flatten()
                        .collect::<Vec<_>>()
                })
                .unwrap_or(vec![]);
            let alt_aa = sites_vec[1].as_str().unwrap_or("NA").to_string();
            let ref_aa = sites_vec[2].as_str().unwrap_or("NA").to_string();
            log::debug!(
                "Evaluating exon3 group alt_aa={alt_aa}, ref_aa={ref_aa}, variant_count={}",
                variants.len()
            );
            //exon3_vars.insert((alt_aa, ref_aa), variants);
            //}
            //for ((alt_aa, ref_aa), variants) in exon3_vars.iter() {
            let num_var_overlap = variants
                .iter()
                .filter(|x| hap_vars.contains(*x))
                .collect::<Vec<_>>()
                .len();
            if num_var_overlap == variants.len() {
                annotated_vars.push(alt_aa.to_string());
            } else {
                annotated_vars.push(ref_aa.to_string());
            }
        }
        Ok(annotated_vars.iter().join(""))
    }

    /// Find the second copy on one allele where the first copy is known
    fn phase_single_allele(
        &mut self,
        hap: String,
        all_haps_on_this_allele: Vec<String>,
        last_copies: &Vec<String>,
        middle_copies: &Vec<String>,
        hap_links: &BTreeMap<String, Vec<String>>,
        alleles: &Vec<Vec<String>>,
        dir_links: &BTreeMap<(String, String), usize>,
        dir_links_loose: &BTreeMap<String, BTreeMap<String, i32>>,
    ) -> Result<Option<Vec<String>>, DError> {
        let remaining_haps = all_haps_on_this_allele
            .clone()
            .iter()
            .filter(|x| **x != hap)
            .map(|x| x.to_string())
            .collect::<Vec<String>>();
        // two copies
        if remaining_haps.len() == 1 {
            if let Some(remaining) = remaining_haps.first() {
                return Ok(Some(vec![hap, remaining.to_string()]));
            }
        }
        let middle_copies_this_allele = all_haps_on_this_allele
            .clone()
            .iter()
            .filter(|x| middle_copies.contains(*x))
            .map(|x| x.to_string())
            .collect::<Vec<String>>();
        // three copies
        let remaining_noending_haps = remaining_haps
            .clone()
            .iter()
            .filter(|x| !last_copies.contains(*x))
            .map(|x| x.to_string())
            .collect::<Vec<String>>();
        // three copies, last one known
        if remaining_noending_haps.len() == 1 {
            if let Some(next_hap) = remaining_noending_haps.first() {
                return Ok(Some(vec![hap, next_hap.to_string()]));
            }
        }
        // three total copies. we don't know the last copy but we know the middle copy
        if remaining_haps.len() == 2 && middle_copies_this_allele.len() == 1 {
            if let Some(next_hap) = middle_copies_this_allele.first() {
                return Ok(Some(vec![hap, next_hap.to_string()]));
            }
        }

        // four copies and we know the one before the last copy
        if last_copies.len() == 1 {
            if let Some(last_copy) = last_copies.first() {
                if let Some(last_copy_before) = hap_links.get(last_copy) {
                    if last_copy_before.len() == 1 {
                        let remaining_noending_not_before_ending_haps = remaining_noending_haps
                            .iter()
                            .filter(|x| !last_copy_before.contains(*x))
                            .collect::<Vec<_>>();
                        if remaining_noending_not_before_ending_haps.len() == 1 {
                            if let Some(next_hap) =
                                remaining_noending_not_before_ending_haps.first()
                            {
                                return Ok(Some(vec![hap, next_hap.to_string()]));
                            }
                        }
                    }
                }
            }
        }
        // four total copies. we know the directional link between the middle two copies
        if remaining_noending_haps.len() == 2 {
            let Some(hap1) = remaining_noending_haps.first().map(ToString::to_string) else {
                return Ok(None);
            };
            let Some(hap2) = remaining_noending_haps.last().map(ToString::to_string) else {
                return Ok(None);
            };
            let link1 = (hap1.clone(), hap2.clone());
            let link2 = (hap2.clone(), hap1.clone());
            if dir_links.contains_key(&link1) && !dir_links.contains_key(&link2) {
                return Ok(Some(vec![hap, hap1.clone()]));
            } else if dir_links.contains_key(&link2) && !dir_links.contains_key(&link1) {
                return Ok(Some(vec![hap, hap2.clone()]));
            }
        }
        // if other haps are phased on this allele
        log::debug!(
            "OPN1LW phase_single_allele context: allele_count={}, last_copy_count={}, remaining_noending_count={}",
            alleles.len(),
            last_copies.len(),
            remaining_noending_haps.len()
        );
        if alleles.len() == 1 && !last_copies.is_empty() {
            let Some(last_copy) = last_copies.first() else {
                return Ok(None);
            };
            let Some(allele) = alleles.first() else {
                return Ok(None);
            };
            if allele.contains(last_copy) {
                let remaining_haps_not_phased = remaining_noending_haps
                    .iter()
                    .filter(|x| !allele.contains(*x))
                    .collect::<Vec<_>>();
                log::debug!(
                    "OPN1LW unphased remaining hap count on allele containing last-copy hap: {}",
                    remaining_haps_not_phased.len()
                );
                if remaining_haps_not_phased.len() == 1 {
                    if let Some(next_hap) = remaining_haps_not_phased.first() {
                        return Ok(Some(vec![hap, next_hap.to_string()]));
                    }
                }
                if remaining_haps_not_phased.len() == 2 {
                    let Some(hap1) = remaining_noending_haps.first().map(ToString::to_string)
                    else {
                        return Ok(None);
                    };
                    let Some(hap2) = remaining_noending_haps.last().map(ToString::to_string) else {
                        return Ok(None);
                    };
                    let link1 = (hap1.clone(), hap2.clone());
                    let link2 = (hap2.clone(), hap1.clone());
                    if dir_links.contains_key(&link1) && !dir_links.contains_key(&link2) {
                        return Ok(Some(vec![hap, hap1.clone()]));
                    } else if dir_links.contains_key(&link2) && !dir_links.contains_key(&link1) {
                        return Ok(Some(vec![hap, hap2.clone()]));
                    }
                }
            }
        }
        //loose links of the first copy
        let hap_next_loose = dir_links_loose.get(&hap);
        if let Some(hap_next_loose) = hap_next_loose {
            let mut hap_next_loose_noending = Vec::new();
            for (hap_candidate, read_count) in hap_next_loose {
                if !last_copies.contains(hap_candidate) && *read_count >= 3 {
                    hap_next_loose_noending.push(hap_candidate.to_string());
                }
            }
            if hap_next_loose_noending.len() == 1 {
                if let Some(hap_next_loose_noending_hap) = hap_next_loose_noending.first() {
                    if all_haps_on_this_allele.contains(hap_next_loose_noending_hap) {
                        return Ok(Some(vec![hap, hap_next_loose_noending_hap.to_string()]));
                    }
                }
            }
        }
        Ok(None)
    }

    pub fn run_opn1lw(&mut self) -> Result<GeneCall, DError> {
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

        let pivot_vars = self
            .locus_config()
            .get("pivot_vars")
            .and_then(|x| x.as_sequence())
            .map(|x| {
                x.iter()
                    .filter_map(|x| x.as_str().map(CandidateSite::from_str))
                    .flatten()
                    .collect::<Vec<_>>()
            })
            .unwrap_or(vec![]);
        let last_copy_vars = self
            .locus_config()
            .get("last_copy_vars")
            .and_then(|x| x.as_sequence())
            .map(|x| {
                x.iter()
                    .filter_map(|x| x.as_str().map(CandidateSite::from_str))
                    .flatten()
                    .collect::<Vec<_>>()
            })
            .unwrap_or(vec![]);

        let (hom_sites_to_add, add_sites) = self.get_sites(&seq, None, Some(0.08))?;
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
        let (mut phase_results, known_del) =
            self.update_indel_and_phase(init_read_hap_map.clone(), &mut call)?;

        // rename haplotypes
        let mut assembled_haps = BTreeMap::new();
        let main_haps_clone = phase_results.assemblies.main_haps.clone();
        let mod_gene_name =
            intersperse(self.gene_name().split_terminator('-'), ",").collect::<String>();
        for (idx, hap) in main_haps_clone.iter().enumerate() {
            assembled_haps.insert(hap.vstr(), format!("{mod_gene_name}_hap{}", idx + 1));
        }
        let haps =
            self.output_variants_in_haps(&phase_results, &known_del, assembled_haps.clone())?;

        // rename haplotypes
        let mut counter_lw = 0;
        let mut counter_mw = 0;
        let mut counter_unknown = 0;
        let mut first_copies = Vec::new();
        let mut last_copies = Vec::new();
        let mut middle_copies = Vec::new();
        let mut annotated_haps = BTreeMap::new();
        let mut renamed_haps = BTreeMap::new();
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
            let mut num_var_overlap = 0;
            let renamed_hap: String;
            let mut gene_annotated: String;
            for var in &pivot_vars {
                if hap_vars.contains(var) {
                    num_var_overlap += 1;
                }
            }
            if num_var_overlap == 0 {
                counter_lw += 1;
                renamed_hap = format!("{mod_gene_name}_opn1lwhap{}", counter_lw);
                gene_annotated = String::from("opn1lw");
            } else if num_var_overlap == 2 {
                counter_mw += 1;
                renamed_hap = format!("{mod_gene_name}_opn1mwhap{}", counter_mw);
                gene_annotated = String::from("opn1mw");
            } else {
                counter_unknown += 1;
                renamed_hap = format!("{mod_gene_name}_opnunknownhap{}", counter_unknown);
                gene_annotated = String::from("opnunknown");
            }
            renamed_haps.insert(hap_name, renamed_hap.clone());
            let hap_seq_first_base = hap_seq.first().ok_or_else(|| {
                crate::phaser::Exception::new(format!(
                    "Haplotype sequence for '{}' is unexpectedly empty",
                    hap_name
                ))
            })?;
            if *hap_seq_first_base != b'x' && *hap_seq_first_base != b'0' {
                first_copies.push(renamed_hap.clone());
            }
            let num_var_last_copy = hap_vars
                .iter()
                .filter(|x| last_copy_vars.contains(*x))
                .collect::<Vec<_>>()
                .len();
            if num_var_last_copy > 0 {
                last_copies.push(renamed_hap.clone());
            } else if !first_copies.contains(&renamed_hap) && !hap_seq.contains(&b'x') {
                middle_copies.push(renamed_hap.clone());
            }
            gene_annotated += "_";
            gene_annotated += self.call_exon3(hap_vars)?.as_str(); //exon3_vars.clone()
            annotated_haps.insert(renamed_hap, gene_annotated);
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

        // Phase alleles
        let allele_result = self.phase_alleles(&mut phase_results, &assembled_haps_renamed, None);
        let hap_links = allele_result.haplotype_links;
        let dir_links = allele_result.directed_links;
        let dir_links_loose = allele_result.directed_links_loose;
        let mut alleles = allele_result.raw_alleles;
        // do not look for two-copy haplotypes for now
        call.total_cn = Some(assembled_haps_renamed.len() as i32);
        let sample_sex = self.settings.sample_sex;
        // incorrect phasing suggests haplotypes with cn > 1
        if self.all_haps_phased_onto_one_allele(&alleles, &assembled_haps_renamed)? {
            alleles = vec![];
            call.total_cn = None;
        }
        // focus on first and second copies
        let mut alleles_1st_2nd = Vec::new();
        let mut phasing_success = false;
        if (sample_sex == Sex::Male && first_copies.len() == 1)
            || (sample_sex == Sex::Female && first_copies.len() == 2)
        {
            for hap in &first_copies {
                let hap_next = hap_links.get(hap);
                if hap_next.is_none() {
                    if last_copies.contains(hap) {
                        alleles_1st_2nd.push(vec![hap.to_string()]);
                    } else {
                        let hap_next_loose = dir_links_loose.get(hap);
                        if let Some(hap_next_loose_value) = hap_next_loose {
                            if hap_next_loose_value.len() == 1 {
                                if let Some((next_hap_name, read_count)) =
                                    hap_next_loose_value.iter().next()
                                {
                                    if *read_count >= 3 {
                                        alleles_1st_2nd
                                            .push(vec![hap.to_string(), next_hap_name.to_string()]);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    let Some(hap_next_value) = hap_next else {
                        continue;
                    };
                    if hap_next_value.len() == 1 {
                        if let Some(hap_next_hap) = hap_next_value.first() {
                            alleles_1st_2nd.push(vec![hap.to_string(), hap_next_hap.to_string()]);
                        }
                    }
                }
            }
            // easier phasing for males as there is one allele
            if sample_sex == Sex::Male && alleles_1st_2nd.is_empty() {
                let Some(first_copy) = first_copies.first().map(ToString::to_string) else {
                    return Ok(call);
                };
                let all_haps_on_this_allele = assembled_haps_renamed
                    .values()
                    .cloned()
                    .collect::<Vec<String>>();
                let allele_found = self.phase_single_allele(
                    first_copy,
                    all_haps_on_this_allele,
                    &last_copies,
                    &middle_copies,
                    &hap_links,
                    &alleles,
                    &dir_links,
                    &dir_links_loose,
                )?;
                if let Some(allele_found) = allele_found {
                    alleles_1st_2nd.push(allele_found);
                }
            }
            // for females, if one allele is completely phased and the other is not phased yet
            if sample_sex == Sex::Female && alleles_1st_2nd.len() == 1 {
                let mut complete_alleles = Vec::new();
                let mut incomplete_alleles = Vec::new();
                for each_allele in &alleles {
                    let this_allele_first_copy = each_allele
                        .iter()
                        .filter(|x| first_copies.contains(*x))
                        .collect::<Vec<_>>();
                    let this_allele_last_copy = each_allele
                        .iter()
                        .filter(|x| last_copies.contains(*x))
                        .collect::<Vec<_>>();
                    if !this_allele_first_copy.is_empty() && !this_allele_last_copy.is_empty() {
                        complete_alleles.push(each_allele.clone());
                    } else {
                        incomplete_alleles.push(each_allele.clone());
                    }
                }
                if complete_alleles.len() == 1 {
                    let Some(complete_allele) = complete_alleles.first() else {
                        return Ok(call);
                    };
                    let other_allele_copies = assembled_haps_renamed
                        .clone()
                        .values()
                        .filter(|x| !complete_allele.contains(*x))
                        .map(|x| x.to_string())
                        .collect::<Vec<_>>();
                    let other_first_copy = first_copies
                        .iter()
                        .filter(|x| !complete_allele.contains(*x))
                        .map(|x| x.to_string())
                        .collect::<Vec<_>>();
                    let other_last_copy = last_copies
                        .iter()
                        .filter(|x| !complete_allele.contains(*x))
                        .map(|x| x.to_string())
                        .collect::<Vec<_>>();
                    log::debug!(
                        "OPN1LW other-allele copy candidates: first={other_first_copy:?}, last={other_last_copy:?}, all={other_allele_copies:?}"
                    );
                    if other_first_copy.len() == 1 && other_last_copy.len() == 1 {
                        let allele_found = self.phase_single_allele(
                            other_first_copy
                                .first()
                                .map(ToString::to_string)
                                .unwrap_or_default(),
                            other_allele_copies,
                            &other_last_copy,
                            &middle_copies,
                            &hap_links,
                            &incomplete_alleles,
                            &dir_links,
                            &dir_links_loose,
                        )?;
                        if let Some(allele_found) = allele_found {
                            alleles_1st_2nd.push(allele_found);
                        }
                    }
                }
            }
            if (sample_sex == Sex::Male && alleles_1st_2nd.len() == 1)
                || (sample_sex == Sex::Female && alleles_1st_2nd.len() == 2)
            {
                // phasing success is defined as the first two copies being phased on each allele
                phasing_success = true;
            } else {
                for hap in &first_copies {
                    let hap_in_allele = alleles_1st_2nd
                        .iter()
                        .map(|x| x.contains(hap))
                        .collect::<Vec<_>>();
                    if !hap_in_allele.contains(&true) {
                        alleles_1st_2nd.push(vec![hap.clone(), String::from("Unknown")]);
                    }
                }
            }
        }

        let mut annotated_alleles = Vec::new();
        for each_allele in &alleles_1st_2nd {
            let mut each_allele_annotated = Vec::new();
            for each_hap in each_allele {
                if each_hap != &String::from("Unknown") {
                    if let Some(annotated_hap) = annotated_haps.get(each_hap) {
                        each_allele_annotated.push(annotated_hap.to_string());
                    } else {
                        each_allele_annotated.push(each_hap.to_string());
                    }
                } else {
                    each_allele_annotated.push(each_hap.to_string());
                }
            }
            annotated_alleles.push(each_allele_annotated.clone());
        }
        let unknown_to_null = |alleles: &[Vec<String>]| -> serde_json::Value {
            serde_json::Value::Array(
                alleles
                    .iter()
                    .map(|allele| {
                        serde_json::Value::Array(
                            allele
                                .iter()
                                .map(|hap| {
                                    if hap == "Unknown" {
                                        serde_json::Value::Null
                                    } else {
                                        serde_json::Value::String(hap.clone())
                                    }
                                })
                                .collect::<Vec<_>>(),
                        )
                    })
                    .collect::<Vec<_>>(),
            )
        };

        // homozygous cases should be rare. Not considered for now.
        call.region_specific_info.insert(
            String::from("alleles_final"),
            unknown_to_null(&alleles_1st_2nd),
        );
        call.region_specific_info.insert(
            String::from("haplotype_links"),
            serde_json::to_value(&hap_links)?,
        );
        call.region_specific_info.insert(
            String::from("raw_alleles"),
            unknown_to_null(&alleles_1st_2nd),
        );

        self.fill_in_call(phase_results, &mut call);
        call.region_specific_info
            .insert(String::from("opn1lw_cn"), counter_lw.into());
        call.region_specific_info
            .insert(String::from("opn1mw_cn"), counter_mw.into());
        call.region_specific_info
            .insert(String::from("first_copies"), first_copies.into());
        call.region_specific_info
            .insert(String::from("middle_copies"), middle_copies.into());
        call.region_specific_info
            .insert(String::from("last_copies"), last_copies.into());
        call.region_specific_info.insert(
            String::from("annotated_haplotypes"),
            serde_json::to_value(&annotated_haps)?,
        );
        call.region_specific_info
            .insert(String::from("phasing_success"), phasing_success.into());
        call.region_specific_info.insert(
            String::from("annotated_alleles"),
            unknown_to_null(&annotated_alleles),
        );
        call.region_specific_info
            .insert(String::from("alleles_all_haplotypes"), alleles.into());

        let mut dir_links_reformat = BTreeMap::new();
        for ((a, b), c) in dir_links.iter() {
            let new_key = format!("{a}-{b}");
            dir_links_reformat.insert(new_key, c);
        }
        call.region_specific_info.insert(
            String::from("directional_links"),
            serde_json::to_value(&dir_links_reformat)?,
        );
        call.region_specific_info.insert(
            String::from("links_loose"),
            serde_json::to_value(&dir_links_loose)?,
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
            "opn1lw",
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
    fn phase_single_allele_matches_python_cases() {
        let outdir = tempfile::TempDir::new().expect("tempdir should build");
        let mut phaser = build_test_phaser(outdir.path());
        let first_copy = String::from("h1");

        let phased = phaser
            .phase_single_allele(
                first_copy.clone(),
                vec![String::from("h1"), String::from("h2")],
                &vec![String::from("h2")],
                &vec![],
                &BTreeMap::new(),
                &vec![],
                &BTreeMap::new(),
                &BTreeMap::new(),
            )
            .expect("phase_single_allele should succeed");
        assert_eq!(phased, Some(vec![String::from("h1"), String::from("h2")]));

        let phased = phaser
            .phase_single_allele(
                first_copy.clone(),
                vec![String::from("h1"), String::from("h2"), String::from("h3")],
                &vec![String::from("h3")],
                &vec![],
                &BTreeMap::new(),
                &vec![],
                &BTreeMap::new(),
                &BTreeMap::new(),
            )
            .expect("phase_single_allele should succeed");
        assert_eq!(phased, Some(vec![String::from("h1"), String::from("h2")]));

        let phased = phaser
            .phase_single_allele(
                first_copy.clone(),
                vec![String::from("h1"), String::from("h2"), String::from("h3")],
                &vec![],
                &vec![String::from("h3")],
                &BTreeMap::new(),
                &vec![],
                &BTreeMap::new(),
                &BTreeMap::new(),
            )
            .expect("phase_single_allele should succeed");
        assert_eq!(phased, Some(vec![String::from("h1"), String::from("h3")]));

        let mut hap_links = BTreeMap::new();
        hap_links.insert(String::from("h4"), vec![String::from("h2")]);
        let phased = phaser
            .phase_single_allele(
                first_copy.clone(),
                vec![
                    String::from("h1"),
                    String::from("h2"),
                    String::from("h3"),
                    String::from("h4"),
                ],
                &vec![String::from("h4")],
                &vec![],
                &hap_links,
                &vec![],
                &BTreeMap::new(),
                &BTreeMap::new(),
            )
            .expect("phase_single_allele should succeed");
        assert_eq!(phased, Some(vec![String::from("h1"), String::from("h3")]));

        let mut dir_links = BTreeMap::new();
        dir_links.insert((String::from("h2"), String::from("h3")), 1usize);
        let phased = phaser
            .phase_single_allele(
                first_copy.clone(),
                vec![
                    String::from("h1"),
                    String::from("h2"),
                    String::from("h3"),
                    String::from("h4"),
                ],
                &vec![String::from("h4")],
                &vec![],
                &BTreeMap::new(),
                &vec![],
                &dir_links,
                &BTreeMap::new(),
            )
            .expect("phase_single_allele should succeed");
        assert_eq!(phased, Some(vec![String::from("h1"), String::from("h2")]));

        let phased = phaser
            .phase_single_allele(
                first_copy.clone(),
                vec![
                    String::from("h1"),
                    String::from("h2"),
                    String::from("h3"),
                    String::from("h4"),
                ],
                &vec![String::from("h4")],
                &vec![],
                &BTreeMap::new(),
                &vec![vec![String::from("h2"), String::from("h4")]],
                &BTreeMap::new(),
                &BTreeMap::new(),
            )
            .expect("phase_single_allele should succeed");
        assert_eq!(phased, Some(vec![String::from("h1"), String::from("h3")]));

        let mut dir_links_loose = BTreeMap::new();
        dir_links_loose.insert(
            String::from("h1"),
            BTreeMap::from([(String::from("h2"), 3i32)]),
        );
        let phased = phaser
            .phase_single_allele(
                first_copy,
                vec![
                    String::from("h1"),
                    String::from("h2"),
                    String::from("h3"),
                    String::from("h4"),
                ],
                &vec![String::from("h4")],
                &vec![],
                &BTreeMap::new(),
                &vec![],
                &BTreeMap::new(),
                &dir_links_loose,
            )
            .expect("phase_single_allele should succeed");
        assert_eq!(phased, Some(vec![String::from("h1"), String::from("h2")]));
    }
}
