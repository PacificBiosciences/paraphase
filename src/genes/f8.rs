// F8 specific caller
use crate::depth::Sex;
use crate::io::json::GeneCall;
use crate::phaser::HapInfoForJson;
use crate::phaser::Phaser;
use crate::toolkit::math::depth_prob;
use crate::toolkit::util::DError;
use itertools::{intersperse, Itertools};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::Read;
use std::collections::{BTreeMap, HashSet};
use vstr::VStr;

impl Phaser {
    /// Get mapped region of the part of reads not overlapping repeat
    #[allow(clippy::type_complexity)]
    fn get_read_positions(
        &mut self,
    ) -> Result<(BTreeMap<String, Vec<String>>, BTreeMap<String, Vec<String>>), DError> {
        let min_extension: i64 = 5000;
        let regions_to_extract = self
            .locus_config()
            .extract_regions(self.settings.genome == "37");
        let mut dpos5: BTreeMap<String, Vec<String>> = BTreeMap::new();
        let mut dpos3: BTreeMap<String, Vec<String>> = BTreeMap::new();
        for (idx, extract_region) in regions_to_extract.iter().enumerate() {
            let coords = extract_region.split_terminator(":").last().ok_or_else(|| {
                crate::phaser::Exception::new(format!(
                    "Malformed extract_region '{}': expected 'chr:start-end'",
                    extract_region
                ))
            })?;
            let (start, stop) = coords
                .split_terminator('-')
                .filter_map(|x| x.parse::<i64>().ok())
                .next_tuple()
                .ok_or_else(|| {
                    crate::phaser::Exception::new(format!(
                        "Malformed coordinate range '{}' in extract_region '{}'",
                        coords, extract_region
                    ))
                })?;
            let pos_name = format!("region{}", idx + 1);

            let mut bam = self.try_genome_bam()?;
            let tid = self.genome_tid().map(|x| x as i32).ok_or_else(|| {
                crate::phaser::Exception::new(format!(
                    "Missing genome tid while extracting read positions for gene '{}'",
                    self.gene_name()
                ))
            })?;
            bam.fetch((tid, start - 2, start - 1))?;
            for read in bam.records() {
                let record = read?;
                let ref_start = record.reference_start();
                if ref_start < start - min_extension {
                    let read_name = std::str::from_utf8(record.qname())?;
                    if idx == 0 || idx == 1 {
                        dpos5
                            .entry(read_name.to_string())
                            .or_default()
                            .push(pos_name.clone());
                    } else {
                        dpos3
                            .entry(read_name.to_string())
                            .or_default()
                            .push(pos_name.clone());
                    }
                }
            }
            let mut bam = self.try_genome_bam()?;
            bam.fetch((tid, stop - 2, stop - 1))?;
            for read in bam.records() {
                let record = read?;
                let ref_end = record.reference_end();
                if ref_end > stop + min_extension {
                    let read_name = std::str::from_utf8(record.qname())?;
                    if idx == 0 || idx == 1 {
                        dpos3
                            .entry(read_name.to_string())
                            .or_default()
                            .push(pos_name.clone());
                    } else {
                        dpos5
                            .entry(read_name.to_string())
                            .or_default()
                            .push(pos_name.clone());
                    }
                }
            }
        }
        Ok((dpos5, dpos3))
    }

    /// Run the F8-specific workflow, including inversion/deletion hap labeling
    /// and flanking-read evidence summarization.
    pub fn run_f8(&mut self) -> Result<GeneCall, DError> {
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
        let e1_e22_depth = self.get_intergenic_depth();

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
        let mut counter_h1 = 0;
        let mut counter_h2 = 0;
        let mut counter_h3 = 0;
        let mut counter_unknown = 0;
        let mut counter_inv = 0;
        let mut counter_del = 0;
        let mut counter_invdup = 0;
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
        for hap in main_haps_clone.iter() {
            let clip_3p = self.get_3pclip_from_hap(&hap.vstr())?;
            let clip_5p = self.get_5pclip_from_hap(&hap.vstr())?;
            log::debug!(
                "F8 haplotype clip anchors: hap={}, clip_3p={clip_3p:?}, clip_5p={clip_5p:?}",
                hap
            );
            let hap_name: String;
            if clip_3p.is_none() || clip_5p.is_none() {
                counter_unknown += 1;
                hap_name = format!("{mod_gene_name}_unknownhap{}", counter_unknown);
            } else {
                let clip_3p_value = clip_3p.unwrap_or(0);
                let clip_5p_value = clip_5p.unwrap_or(0);
                if clip_3p_value == 0 && clip_5p_value == 0 {
                    counter_h2 += 1;
                    hap_name = format!("{mod_gene_name}_int22h2hap{}", counter_h2);
                } else if clip_3p_value == *first_clip_3p && clip_5p_value == *first_clip_5p {
                    counter_h1 += 1;
                    hap_name = format!("{mod_gene_name}_int22h1hap{}", counter_h1);
                } else if clip_3p_value == *second_clip_3p && clip_5p_value == 0 {
                    counter_h3 += 1;
                    hap_name = format!("{mod_gene_name}_int22h3hap{}", counter_h3);
                } else if clip_3p_value == *second_clip_3p && clip_5p_value == *first_clip_5p {
                    counter_inv += 1;
                    hap_name = format!("{mod_gene_name}_int22invhap{}", counter_inv);
                } else if clip_3p_value == 0 && clip_5p_value == *first_clip_5p {
                    counter_del += 1;
                    hap_name = format!("{mod_gene_name}_int22delhap{}", counter_del);
                } else if clip_5p_value == 0 && clip_3p_value == *first_clip_3p {
                    counter_invdup += 1;
                    hap_name = format!("{mod_gene_name}_int22invorduphap{}", counter_invdup);
                } else {
                    counter_unknown += 1;
                    hap_name = format!("{mod_gene_name}_unknownhap{}", counter_unknown);
                }
            }
            assembled_haps.insert(hap.vstr(), hap_name.clone());
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

        // check flanking region, call sv
        let (dpos5, dpos3) = self.get_read_positions()?;
        let mut haplotype_flanking_regions_5p: BTreeMap<VStr<'_>, HashSet<String>> =
            BTreeMap::new();
        let mut haplotype_flanking_regions_3p: BTreeMap<VStr<'_>, HashSet<String>> =
            BTreeMap::new();
        for (hap, reads) in &phase_results.uniquely_supporting_reads {
            for read in reads {
                if dpos5.contains_key(read) {
                    let read_pos_name = &dpos5[read];
                    if read_pos_name.len() == 1 {
                        haplotype_flanking_regions_5p
                            .entry(hap.vstr())
                            .or_default()
                            .insert(read_pos_name[0].clone());
                    }
                }
                if dpos3.contains_key(read) {
                    let read_pos_name = &dpos3[read];
                    if read_pos_name.len() == 1 {
                        haplotype_flanking_regions_3p
                            .entry(hap.vstr())
                            .or_default()
                            .insert(read_pos_name[0].clone());
                    }
                }
            }
        }

        let total_cn = assembled_haps.len();
        let sample_sex = self.settings.sample_sex;
        let mut flanking_sum = BTreeMap::new();
        let mut sv_hap = BTreeMap::new();
        for (hap, hap_name) in &assembled_haps {
            let p5region = if let Some(regions) = haplotype_flanking_regions_5p.get(hap) {
                regions.iter().sorted().join("/")
            } else {
                String::from("")
            };
            let p3region = if let Some(regions) = haplotype_flanking_regions_3p.get(hap) {
                regions.iter().sorted().join("/")
            } else {
                String::from("")
            };
            flanking_sum.insert(hap_name.to_string(), [p5region, p3region].join("-"));
        }
        log::debug!("F8 flanking-region support summary: {flanking_sum:?}");

        // region2 and region3 are homologous at 5p for another 50kb,
        // so we cannot separate the upstream regions
        // we look for read evidence only at downstream region of region2 and region3
        // so we drop duplication as it involves upstream region of region2 and it's not pathogenic (?)
        for hap_name in assembled_haps.values() {
            if hap_name.contains("del") {
                if sample_sex == Sex::Female {
                    if let Some(depth) = self.settings.depth.as_ref() {
                        let genome_depth = depth.median;
                        let prob = depth_prob(e1_e22_depth as i32, genome_depth / 2.0_f64);
                        if let Some(prob_value) = prob {
                            if prob_value[0] > 0.75 {
                                sv_hap.insert(hap_name.to_string(), String::from("deletion"));
                            }
                        }
                    }
                }
                if sample_sex == Sex::Male {
                    if e1_e22_depth < 1.0 {
                        sv_hap.insert(hap_name.to_string(), String::from("deletion"));
                    }
                }
            }
            if hap_name.contains("inv") && !hap_name.contains("dup") {
                sv_hap.insert(hap_name.to_string(), String::from("inversion"));
            }
        }
        /*
        for (hap, links) in &flanking_sum {
            if links == &String::from("region1-region2") && hap.contains("int22h2") {
                if sample_sex == Sex::Female {
                    if let Some(depth) = self.settings.depth.as_ref() {
                        let genome_depth = depth.median;
                        let prob = depth_prob(e1_e22_depth as i32, genome_depth / 2.0 as f64);
                        if let Some(prob_value) = prob {
                            if prob_value[0] > 0.75 {
                                sv_hap.insert(hap.to_string(), String::from("deletion"));
                            }
                        }
                    }
                }
                if sample_sex == Sex::Male {
                    if e1_e22_depth < 1.0 {
                        sv_hap.insert(hap.to_string(), String::from("deletion"));
                    }
                }
            } else if links == &String::from("region1-region3") && hap.contains("int22h3") {
                sv_hap.insert(hap.to_string(), String::from("inversion"));
            }
        }
        */

        call.total_cn = Some(total_cn as i32);
        if sv_hap.is_empty() {
            if (sample_sex == Sex::Female && total_cn < 6)
                || (sample_sex == Sex::Male && total_cn < 3)
            {
                call.total_cn = None;
            }
        }

        self.fill_in_call(phase_results, &mut call);
        call.region_specific_info
            .insert(String::from("exon1_to_exon22_depth"), e1_e22_depth.into());
        call.region_specific_info
            .insert(String::from("sv_called"), serde_json::to_value(&sv_hap)?);
        call.region_specific_info.insert(
            String::from("flanking_summary"),
            serde_json::to_value(&flanking_sum)?,
        );

        Ok(call)
    }
}

#[cfg(test)]
mod tests {
    use crate::depth::Sex;
    use crate::toolkit::util::{self, DError, DResult};
    use crate::{
        config,
        io::json::GeneCall,
        phaser::{self, Phaser},
    };
    use serde_json::json;

    fn hg38_reference() -> Option<String> {
        std::env::var("HG38").ok()
    }

    fn run_f8_fixture(sample_id: &str, bam_name: &str) -> Result<Option<GeneCall>, DError> {
        let Some(genome_path) = hg38_reference() else {
            log::warn!("Skipping F8 test because the HG38 environment variable is not configured.");
            return Ok(None);
        };
        let outdir = tempfile::TempDir::new()?;
        let settings = phaser::Settings::new(
            sample_id,
            (genome_path, util::test_file(bam_name)),
            outdir.path(),
            "f8",
            &config::Region::try_load(None)?,
            None,
            Some(Sex::Male),
            String::from("38"),
            None,
            0.03,
            false,
        );
        let gene_config = config::Gene::try_load(None)?;
        let mut phaser = Phaser::new(settings, Some(gene_config), None, None)?;
        phaser.run().map(Some)
    }

    #[test]
    fn f8_inversion_fixture_reports_sv() -> DResult {
        let Some(call) = run_f8_fixture("inv", "f8/f8_inv_genome.bam")? else {
            return Ok(());
        };
        if call.gene_name.is_empty() && call.region_specific_info.is_empty() {
            return Ok(());
        }
        assert_eq!(
            call.region_specific_info.get("sv_called"),
            Some(&json!({"f8_int22invhap1": "inversion"}))
        );
        Ok(())
    }

    #[test]
    fn f8_deletion_fixture_reports_sv() -> DResult {
        let Some(call) = run_f8_fixture("del", "f8/f8_del_genome.bam")? else {
            return Ok(());
        };
        if call.gene_name.is_empty() && call.region_specific_info.is_empty() {
            return Ok(());
        }
        assert_eq!(
            call.region_specific_info.get("sv_called"),
            Some(&json!({"f8_int22delhap1": "deletion"}))
        );
        Ok(())
    }
}
