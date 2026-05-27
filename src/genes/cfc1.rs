// CFC1 specific caller
use crate::io::json::GeneCall;
use crate::phaser::HapInfoForJson;
use crate::phaser::Phaser;
use crate::toolkit::util::DError;
use itertools::intersperse;
use std::collections::BTreeMap;

impl Phaser {
    pub fn run_cfc1(&mut self) -> Result<GeneCall, DError> {
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

        // rename haplotypes
        let mut assembled_haps = BTreeMap::new();
        let main_haps_clone = phase_results.assemblies.main_haps.clone();
        let mod_gene_name =
            intersperse(self.gene_name().split_terminator('-'), ",").collect::<String>();
        for (idx, hap) in main_haps_clone.iter().enumerate() {
            assembled_haps.insert(hap.vstr(), format!("{mod_gene_name}_hap{}", idx + 1));
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

        let (two_cp_haps, mut total_cn) = self.adjust_depth(
            assembled_haps.clone(),
            haps,
            phase_results.clone(),
            true,
            false,
            true,
            0.15,
        )?;
        if self.init_het_sites.is_empty() {
            total_cn = 4;
        }
        call.total_cn = if total_cn < 4 {
            None
        } else {
            Some(total_cn as i32)
        };
        call.two_copy_haplotypes = two_cp_haps;

        self.fill_in_call(phase_results, &mut call);
        Ok(call)
    }
}
