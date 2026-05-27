// NEB specific caller
use crate::io::json::GeneCall;
use crate::phaser::HapInfoForJson;
use crate::phaser::Phaser;
use crate::toolkit::util::DError;
use itertools::intersperse;
use std::collections::BTreeMap;

impl Phaser {
    pub fn run_neb(&mut self) -> Result<GeneCall, DError> {
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

        // assign tri-1,2,3
        let mut tri1 = Vec::new();
        let mut tri2 = Vec::new();
        let mut tri3 = Vec::new();
        let hap_len = self.het_sites.len();
        let nsite = std::cmp::min(hap_len / 2, 10);
        for (hap_seq, hap_name) in &assembled_haps {
            let start_seq = &hap_seq[..nsite];
            let end_seq = &hap_seq[(hap_len - nsite)..];
            let start_seq_count1 = start_seq
                .iter()
                .filter(|x| **x == b'1')
                .collect::<Vec<_>>()
                .len();
            let start_seq_count2 = start_seq
                .iter()
                .filter(|x| **x == b'2')
                .collect::<Vec<_>>()
                .len();
            let end_seq_count1 = end_seq
                .iter()
                .filter(|x| **x == b'1')
                .collect::<Vec<_>>()
                .len();
            let end_seq_count2 = end_seq
                .iter()
                .filter(|x| **x == b'2')
                .collect::<Vec<_>>()
                .len();
            if start_seq_count1 >= start_seq_count2 && end_seq_count1 >= end_seq_count2 {
                tri1.push(hap_name.to_string());
            } else if start_seq_count1 < start_seq_count2 && end_seq_count1 < end_seq_count2 {
                tri3.push(hap_name.to_string());
            } else {
                tri2.push(hap_name.to_string());
            }
        }

        let mut two_cp_haps = Vec::new();
        let nhap = assembled_haps.len();
        if nhap == 3 && tri1.len() == 1 && tri2.len() == 1 && tri3.len() == 1 {
            for hap in assembled_haps.values() {
                two_cp_haps.push(hap.to_string());
            }
        } else if nhap < 6 && nhap > 1 {
            two_cp_haps = self.compare_depth(&haps, &assembled_haps, true, false)?;
            if two_cp_haps.is_empty() && !phase_results.read_counts.0.is_empty() {
                two_cp_haps =
                    self.compare_depth_by_read_count(&assembled_haps, &phase_results, 0.15, &[]);
            }
        }
        for hap in &two_cp_haps {
            if tri1.contains(hap) {
                tri1.push(hap.to_string());
            }
            if tri2.contains(hap) {
                tri2.push(hap.to_string());
            }
            if tri3.contains(hap) {
                tri3.push(hap.to_string());
            }
        }
        if tri1.len() == 1 {
            let tri1_first = tri1[0].clone();
            two_cp_haps.push(tri1_first.clone());
            tri1.push(tri1_first);
        }
        if tri3.len() == 1 {
            let tri3_first = tri3[0].clone();
            two_cp_haps.push(tri3_first.clone());
            tri3.push(tri3_first);
        }

        // Phase alleles
        let mut alleles: Vec<Vec<String>> = Vec::new();
        let mut raw_alleles: Vec<Vec<String>> = Vec::new();
        let mut hap_links: BTreeMap<String, Vec<String>> = BTreeMap::new();
        if two_cp_haps.is_empty() {
            let allele_result = self.phase_alleles(&mut phase_results, &assembled_haps, None);
            hap_links = allele_result.haplotype_links;
            alleles = allele_result.alleles;
            raw_alleles = allele_result.raw_alleles;
        }

        let total_cn = nhap + two_cp_haps.len();
        call.total_cn = if total_cn == 0 {
            None
        } else {
            Some(total_cn as i32)
        };
        // incorrect phasing suggests haplotypes with cn > 1
        if self.all_haps_phased_onto_one_allele(&alleles, &assembled_haps)? {
            alleles = vec![];
            raw_alleles = vec![];
            call.total_cn = None;
        }
        if tri1.len() > 2 || tri3.len() > 2 {
            alleles = vec![];
            call.total_cn = None;
        }
        //if tri2 has links to both haps in tri1 or tri3
        if two_cp_haps.is_empty() && tri1.len() == 2 && tri2.len() == 1 && tri3.len() == 2 {
            let tri2_hap = &tri2[0];
            if hap_links.contains_key(tri2_hap) {
                let tri2_links = &hap_links[tri2_hap];
                let mut link_to_tri1 = Vec::new();
                let mut link_to_tri3 = Vec::new();
                for hap in &tri1 {
                    if tri2_links.contains(hap) {
                        link_to_tri1.push(hap);
                    }
                }
                for hap in &tri3 {
                    if tri2_links.contains(hap) {
                        link_to_tri3.push(hap);
                    }
                }
                if link_to_tri1.len() > 1 || link_to_tri3.len() > 1 {
                    call.total_cn = Some(6);
                    two_cp_haps = vec![tri2_hap.to_string()];
                    alleles = vec![];
                    raw_alleles = vec![];
                }
            }
        }

        call.region_specific_info
            .insert(String::from("alleles_final"), alleles.into());
        call.region_specific_info.insert(
            String::from("haplotype_links"),
            serde_json::to_value(&hap_links)?,
        );
        call.region_specific_info
            .insert(String::from("raw_alleles"), raw_alleles.into());
        call.two_copy_haplotypes = two_cp_haps;

        self.fill_in_call(phase_results, &mut call);
        let mut repeat_name = BTreeMap::new();
        repeat_name.insert(String::from("tri1"), tri1);
        repeat_name.insert(String::from("tri2"), tri2);
        repeat_name.insert(String::from("tri3"), tri3);
        call.region_specific_info.insert(
            String::from("repeat_name"),
            serde_json::to_value(&repeat_name)?,
        );

        Ok(call)
    }
}
