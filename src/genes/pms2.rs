// PMS2 specific caller
use crate::io::json::GeneCall;
use crate::phaser::HapInfoForJson;
use crate::phaser::Phaser;
use crate::toolkit::util::DError;
use itertools::intersperse;
use std::collections::BTreeMap;

impl Phaser {
    pub fn run_pms2(&mut self) -> Result<GeneCall, DError> {
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

        self.settings.big_deletion_settings.min_size = 1000;
        let (hom_sites_to_add, add_sites) = self.get_sites(&seq, Some(6000), None)?;
        let pms2cl_clip_site = *self.clip_3p_positions.first().ok_or_else(|| {
            crate::phaser::Exception::new(format!(
                "clip_3p_positions is empty for gene '{}'",
                self.gene_name()
            ))
        })?;
        self.find_clip_site(None, None, None)?;
        self.clip_3p_positions = vec![pms2cl_clip_site];
        log::debug!(
            "Using PMS2CL clip anchor for {}: clip_3p_positions={:?}, clip_5p_positions={:?}",
            self.gene_name(),
            self.clip_3p_positions,
            self.clip_5p_positions
        );

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
                /* min_clip_len */ Some(100u32),
            ),
            tid,
            Some(50),
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
        let mut counter_unknown = 0;
        for hap in main_haps_clone.iter() {
            let clip_3p = self.get_3pclip_from_hap(&hap.vstr())?;
            if clip_3p.is_none() {
                let mut num_pms2_bases = 0;
                for (i, base) in hap.iter().enumerate() {
                    let var_pos = self.het_sites[i].pos;
                    if var_pos > pms2cl_clip_site + 50 && base != &b'0' && base != &b'x' {
                        num_pms2_bases += 1;
                    }
                }
                if num_pms2_bases >= 3 {
                    counter_gene += 1;
                    assembled_haps.insert(
                        hap.vstr(),
                        format!("{mod_gene_name}_pms2hap{}", counter_gene),
                    );
                } else {
                    counter_unknown += 1;
                    assembled_haps.insert(
                        hap.vstr(),
                        format!("{mod_gene_name}_unknownhap{}", counter_unknown),
                    );
                }
            }
            if let Some(clip_3p_value) = clip_3p {
                if self.clip_3p_positions.contains(&clip_3p_value) {
                    counter_pseudo += 1;
                    assembled_haps.insert(
                        hap.vstr(),
                        format!("{mod_gene_name}_pms2clhap{}", counter_pseudo),
                    );
                } else {
                    //assert_eq!(clip_3p_value, 0);
                    counter_gene += 1;
                    assembled_haps.insert(
                        hap.vstr(),
                        format!("{mod_gene_name}_pms2hap{}", counter_gene),
                    );
                }
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

        let mut total_cn = assembled_haps.len();
        let mut two_cp_haps = Vec::new();
        if total_cn < 4 && counter_unknown == 0 {
            if counter_gene == 1 && counter_pseudo == 1 {
                two_cp_haps = assembled_haps.values().cloned().collect::<Vec<_>>();
            } else if total_cn == 3 {
                if counter_gene == 2 && counter_pseudo == 1 {
                    for hap in assembled_haps.values() {
                        if hap.contains("pms2cl") {
                            two_cp_haps.push(hap.to_string());
                        }
                    }
                } else if counter_gene == 1 && counter_pseudo == 2 {
                    for hap in assembled_haps.values() {
                        if hap.contains("pms2hap") {
                            two_cp_haps.push(hap.to_string());
                        }
                    }
                }
            }
        }

        total_cn = assembled_haps.len() + two_cp_haps.len();
        let pms2_cn = assembled_haps
            .values()
            .filter(|x| !x.contains("cl") && !x.contains("unknown"))
            .collect::<Vec<_>>()
            .len()
            + two_cp_haps
                .iter()
                .filter(|x| !x.contains("cl") && !x.contains("unknown"))
                .collect::<Vec<_>>()
                .len();

        call.total_cn = Some(total_cn as i32);
        let mut gene_cn = Some(pms2_cn as i32);
        call.two_copy_haplotypes = two_cp_haps;

        if pms2_cn != 2 || counter_unknown > 0 {
            gene_cn = None;
        }
        // homozygous case
        if total_cn == 0 {
            call.total_cn = None;
        }

        self.fill_in_call(phase_results, &mut call);
        call.region_specific_info
            .insert(String::from("gene_cn"), gene_cn.into());

        Ok(call)
    }
}
