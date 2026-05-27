// IKBKG specific caller
use crate::depth::Sex;
use crate::io::json::GeneCall;
use crate::phaser::HapInfoForJson;
use crate::phaser::Phaser;
use crate::toolkit::util::DError;
use itertools::intersperse;
use std::collections::BTreeMap;

impl Phaser {
    /// Run IKBKG-specific phasing/copy-number workflow with pseudogene and
    /// deletion-aware haplotype interpretation.
    pub fn run_ikbkg(&mut self) -> Result<GeneCall, DError> {
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

        let (hom_sites_to_add, add_sites) = self.get_sites(&seq, Some(6000), Some(0.095))?;
        let tid = self.genome_tid().map(|x| x as i32).ok_or_else(|| {
            crate::phaser::Exception::new(format!(
                "Missing chromosome tid while running '{}' workflow",
                self.gene_name()
            ))
        })?;
        let del_reads = self.del_data[0].del_reads_partial.clone();
        let del_reads_count = del_reads.len();
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
        let mut counter_gene = 0;
        let mut counter_pseudo = 0;
        let mut counter_unknown = 0;
        let mut counter_dup = 0;
        let mut deletion_haplotypes = Vec::new();
        let first_clip_5p = self.clip_5p_positions.first().ok_or_else(|| {
            crate::phaser::Exception::new(format!(
                "clip_5p_positions is empty for gene '{}'",
                self.gene_name()
            ))
        })?;
        let second_clip_5p = self.clip_5p_positions.last().ok_or_else(|| {
            crate::phaser::Exception::new(format!(
                "clip_5p_positions is empty for gene '{}'",
                self.gene_name()
            ))
        })?;
        for hap in main_haps_clone.iter() {
            let clip_5p = self.get_5pclip_from_hap(&hap.vstr())?;
            log::debug!(
                "Classifying IKBKG haplotype by 5' clip support: hap={}, clip_5p={clip_5p:?}",
                hap
            );
            let hap_name: String;
            if clip_5p.is_none() {
                counter_unknown += 1;
                hap_name = format!("{mod_gene_name}_unknownhap{}", counter_unknown);
            } else if let Some(clip_5p_value) = clip_5p {
                if clip_5p_value == *first_clip_5p {
                    counter_pseudo += 1;
                    hap_name = format!("{mod_gene_name}_pseudohap{}", counter_pseudo);
                } else if clip_5p_value == *second_clip_5p {
                    counter_dup += 1;
                    hap_name = format!("{mod_gene_name}_duphap{}", counter_dup);
                } else {
                    assert_eq!(clip_5p_value, 0);
                    counter_gene += 1;
                    hap_name = format!("{mod_gene_name}_ikbkghap{}", counter_gene);
                }
            } else {
                counter_unknown += 1;
                hap_name = format!("{mod_gene_name}_unknownhap{}", counter_unknown);
            }
            assembled_haps.insert(hap.vstr(), hap_name.clone());
            if hap.contains(&b'3') {
                deletion_haplotypes.push(hap_name.to_string());
            }
        }
        call.final_haplotypes = assembled_haps
            .clone()
            .into_iter()
            .map(|(k, v)| (k.to_string(), v))
            .collect::<BTreeMap<_, _>>();

        // Output variants
        let haps =
            self.output_variants_in_haps(&phase_results, &known_del, assembled_haps.clone())?;
        call.haplotype_details = haps
            .iter()
            .map(|(key, val)| (key.clone(), HapInfoForJson::from(val)))
            .collect::<BTreeMap<_, _>>();

        // phase alleles
        let allele_result = self.phase_alleles(&mut phase_results, &assembled_haps, None);
        call.region_specific_info.insert(
            String::from("haplotype_links"),
            serde_json::to_value(&allele_result.haplotype_links)?,
        );

        let mut total_cn = assembled_haps.len();
        let mut two_cp_haps = Vec::new();
        let sample_sex = self.settings.sample_sex;
        if counter_unknown == 0 && sample_sex == Sex::Female {
            if counter_gene == 1 && counter_pseudo == 1 {
                for hap in assembled_haps.values() {
                    if !hap.contains("dup") {
                        two_cp_haps.push(hap.to_string());
                    }
                }
            } else if (counter_gene > 1 && counter_pseudo == 1)
                || (counter_gene == 1 && counter_pseudo > 1)
            {
                two_cp_haps = self.compare_depth(&haps, &assembled_haps, true, false)?;
                if two_cp_haps.is_empty() && !phase_results.read_counts.0.is_empty() {
                    two_cp_haps = self.compare_depth_by_read_count(
                        &assembled_haps,
                        &phase_results,
                        0.15,
                        &[],
                    );
                }
            }
        }
        for hap in &two_cp_haps {
            total_cn += 1;
            if hap.contains("ikbkghap") {
                counter_gene += 1;
            }
        }
        let mut gene_cn = Some(counter_gene);
        call.total_cn = Some(total_cn as i32);
        call.two_copy_haplotypes = two_cp_haps;
        if counter_gene == 1
            && counter_pseudo != 1
            && counter_unknown == 0
            && sample_sex == Sex::Female
        {
            gene_cn = None;
            call.total_cn = None;
        }
        if counter_gene == 0 || counter_pseudo == 0 || total_cn == 0 {
            gene_cn = None;
            call.total_cn = None;
        }

        // this is on chrX, males have one copy of gene and one copy of pseudogene
        if sample_sex == Sex::Male {
            if counter_unknown == 0 && (counter_gene > 1 || counter_pseudo > 1) {
                gene_cn = None;
                call.total_cn = None;
            }
        }
        // all haplotypes are phased into one allele
        let alleles = allele_result.alleles;
        let mut raw_alleles = allele_result.raw_alleles;
        if self.all_haps_phased_onto_one_allele(&alleles, &assembled_haps)? {
            raw_alleles = vec![];
        }
        call.region_specific_info
            .insert(String::from("raw_alleles"), raw_alleles.into());

        self.fill_in_call(phase_results, &mut call);
        // additional fields to report
        call.region_specific_info.insert(
            String::from("deletion_haplotypes"),
            deletion_haplotypes.into(),
        );
        call.region_specific_info
            .insert(String::from("del_read_number"), del_reads_count.into());
        call.region_specific_info
            .insert(String::from("gene_cn"), gene_cn.into());
        Ok(call)
    }
}
