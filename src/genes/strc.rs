// STRC specific caller
use crate::io::json::GeneCall;
use crate::phaser::HapInfoForJson;
use crate::phaser::Phaser;
use crate::toolkit::math::depth_prob;
use crate::toolkit::util::DError;
use itertools::intersperse;
use std::collections::BTreeMap;

impl Phaser {
    /// Compute median depth over STRC intergenic normalization intervals.
    pub fn get_intergenic_depth(&self) -> f32 {
        let depth_region = self.locus_config().depth_region();
        let Ok(mut bam) = self.try_genome_bam() else {
            log::warn!(
                "Failed to open genome BAM while computing intergenic depth for gene {}",
                self.gene_name()
            );
            return f32::NAN;
        };
        let Some(tid) = self.genome_tid().map(|x| x as i32) else {
            log::warn!(
                "Missing genome tid while computing intergenic depth for gene {}",
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
        let (intergenic_region_depth, _percentile) = region_depth[0];
        intergenic_region_depth
    }

    /// Run STRC-specific workflow with pseudogene-aware hap naming and depth-driven
    /// copy-number adjustment.
    pub fn run_strc(&mut self) -> Result<GeneCall, DError> {
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

        // check intergenic depth
        let intergenic_depth = self.get_intergenic_depth();

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

        // rename haplotypes
        let mut assembled_haps = BTreeMap::new();
        let main_haps_clone = phase_results.assemblies.main_haps.clone();
        let mod_gene_name =
            intersperse(self.gene_name().split_terminator('-'), ",").collect::<String>();
        let mut counter_gene = 0;
        let mut counter_pseudo = 0;
        for hap in main_haps_clone.iter() {
            if hap.contains(&b'3') {
                counter_pseudo += 1;
                assembled_haps.insert(
                    hap.vstr(),
                    format!("{mod_gene_name}_strcp1hap{}", counter_pseudo),
                );
            } else {
                counter_gene += 1;
                assembled_haps.insert(
                    hap.vstr(),
                    format!("{mod_gene_name}_strchap{}", counter_gene),
                );
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

        let mut two_cp_haps = Vec::new();
        // one haplotype, identical on both alleles
        if assembled_haps.len() == 1 && self.init_het_sites.is_empty() {
            two_cp_haps = assembled_haps.values().cloned().collect::<Vec<_>>();
        } else if counter_gene == 1 || counter_pseudo == 1 {
            two_cp_haps = self.compare_depth(&haps, &assembled_haps, false, false)?;
        }
        if intergenic_depth > 5.0
            && counter_gene == 1
            && counter_pseudo == 1
            && two_cp_haps.is_empty()
        {
            two_cp_haps = assembled_haps.values().cloned().collect::<Vec<_>>();
        } else if two_cp_haps.is_empty() && counter_gene == 1 && counter_pseudo > 1 {
            two_cp_haps =
                self.compare_depth_by_read_count(&assembled_haps, &phase_results, 0.15, &[]);
            two_cp_haps = two_cp_haps
                .iter()
                .filter(|x| !x.contains("strcp1"))
                .cloned()
                .collect::<Vec<_>>();
        }
        for hap in &two_cp_haps {
            if hap.contains("strcp1") {
                counter_pseudo += 1;
            } else {
                counter_gene += 1;
            }
        }

        let total_cn = assembled_haps.len() + two_cp_haps.len();
        call.total_cn = Some(total_cn as i32);
        let mut gene_cn = Some(counter_gene);
        call.two_copy_haplotypes = two_cp_haps;
        // check depth between STRC and pseudogene
        if let Some(depth) = self.settings.depth.as_ref() {
            let genome_depth = depth.median;
            let prob = depth_prob(intergenic_depth as i32, genome_depth / 2.0_f64);
            if let Some(prob_value) = prob {
                log::debug!(
                    "STRC depth comparison against genome coverage: prob_value={prob_value:?}, intergenic_depth={intergenic_depth}, genome_depth={genome_depth}"
                );
                if prob_value[0] < 0.9 && counter_gene == 1 && counter_pseudo == 2 {
                    gene_cn = None;
                    call.total_cn = None;
                }
                if prob_value[0] > 0.95 && counter_gene > 1 && counter_pseudo > 1 {
                    gene_cn = None;
                    call.total_cn = None;
                }
            }
        }

        self.fill_in_call(phase_results, &mut call);
        // additional fields to report
        call.region_specific_info
            .insert(String::from("intergenic_depth"), intergenic_depth.into());
        call.region_specific_info
            .insert(String::from("gene_cn"), gene_cn.into());
        Ok(call)
    }
}
