// STRC specific caller
use crate::io::json::GeneCall;
use crate::phaser::HapInfoForJson;
use crate::phaser::Phaser;
use crate::toolkit::math::depth_prob;
use crate::toolkit::util::DError;
use itertools::intersperse;
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum StrcHaplotypeIdentity {
    Strc,
    Strcp1,
    Discordant,
    Unresolved,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum StrcDeletionEvidence {
    Absent,
    Present,
    Unresolved,
}

fn strc_deletion_evidence(
    haplotype: &[u8],
    synthetic_index: Option<usize>,
    internal_indices: &[usize],
) -> StrcDeletionEvidence {
    if let Some(index) = synthetic_index {
        return match haplotype.get(index) {
            Some(b'1') => StrcDeletionEvidence::Absent,
            Some(b'3') => StrcDeletionEvidence::Present,
            _ => StrcDeletionEvidence::Unresolved,
        };
    }

    let mut deletion_present = false;
    let mut deletion_absent = false;
    for base in internal_indices
        .iter()
        .filter_map(|index| haplotype.get(*index))
    {
        match base {
            b'3' => deletion_present = true,
            b'1' | b'2' => deletion_absent = true,
            _ => {}
        }
    }
    match (deletion_present, deletion_absent) {
        (true, false) => StrcDeletionEvidence::Present,
        (false, true) => StrcDeletionEvidence::Absent,
        _ => StrcDeletionEvidence::Unresolved,
    }
}

fn classify_strc_haplotype(
    haplotype: &[u8],
    pivot_index: Option<usize>,
    deletion_evidence: StrcDeletionEvidence,
) -> StrcHaplotypeIdentity {
    match (
        pivot_index.and_then(|index| haplotype.get(index)),
        deletion_evidence,
    ) {
        (Some(b'1'), StrcDeletionEvidence::Absent) => StrcHaplotypeIdentity::Strc,
        (Some(b'2'), StrcDeletionEvidence::Present) => StrcHaplotypeIdentity::Strcp1,
        (Some(b'1'), StrcDeletionEvidence::Present)
        | (Some(b'2'), StrcDeletionEvidence::Absent) => StrcHaplotypeIdentity::Discordant,
        _ => StrcHaplotypeIdentity::Unresolved,
    }
}

fn confident_identity_from_haplotype_name(name: &str) -> Option<StrcHaplotypeIdentity> {
    if name.contains("_strcp1hap") {
        Some(StrcHaplotypeIdentity::Strcp1)
    } else if name.contains("_strchap") {
        Some(StrcHaplotypeIdentity::Strc)
    } else {
        None
    }
}

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
        let pivot_index = self
            .get_pivot_index()
            .and_then(|index| usize::try_from(index).ok());
        let known_strc_deletion = self.del_data.first().filter(|deletion| {
            known_del
                .get(&'3')
                .is_some_and(|known_name| known_name == &deletion.name())
        });
        let synthetic_deletion_index = known_strc_deletion.and_then(|deletion| {
            let deletion_name = deletion.name();
            self.het_sites
                .iter()
                .position(|site| site.to_string() == deletion_name)
        });
        let internal_deletion_indices = known_strc_deletion
            .map(|deletion| {
                self.het_sites
                    .iter()
                    .enumerate()
                    .filter_map(|(index, site)| {
                        (site.pos > deletion.threep().start && site.pos <= deletion.fivep().end)
                            .then_some(index)
                    })
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();
        let mut counter_gene = 0;
        let mut counter_pseudo = 0;
        let mut counter_unknown = 0;
        for hap in main_haps_clone.iter() {
            let deletion_evidence =
                strc_deletion_evidence(hap, synthetic_deletion_index, &internal_deletion_indices);
            let identity = classify_strc_haplotype(hap, pivot_index, deletion_evidence);
            log::debug!(
                "Classifying STRC haplotype from configured markers: hap={hap}, pivot_index={pivot_index:?}, synthetic_deletion_index={synthetic_deletion_index:?}, internal_deletion_indices={internal_deletion_indices:?}, deletion_evidence={deletion_evidence:?}, identity={identity:?}"
            );
            let hap_name = match identity {
                StrcHaplotypeIdentity::Strc => {
                    counter_gene += 1;
                    format!("{mod_gene_name}_strchap{}", counter_gene)
                }
                StrcHaplotypeIdentity::Strcp1 => {
                    counter_pseudo += 1;
                    format!("{mod_gene_name}_strcp1hap{}", counter_pseudo)
                }
                StrcHaplotypeIdentity::Discordant | StrcHaplotypeIdentity::Unresolved => {
                    counter_unknown += 1;
                    format!("{mod_gene_name}_unknownhap{}", counter_unknown)
                }
            };
            assembled_haps.insert(hap.vstr(), hap_name);
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
        } else if counter_gene == 1 || counter_pseudo == 1 || counter_unknown > 0 {
            two_cp_haps = self.compare_depth(&haps, &assembled_haps, false, false)?;
        }
        if intergenic_depth > 5.0
            && counter_gene == 1
            && counter_pseudo == 1
            && counter_unknown == 0
            && two_cp_haps.is_empty()
        {
            two_cp_haps = assembled_haps.values().cloned().collect::<Vec<_>>();
        } else if two_cp_haps.is_empty() && counter_gene == 1 && counter_pseudo > 1 {
            two_cp_haps =
                self.compare_depth_by_read_count(&assembled_haps, &phase_results, 0.15, &[]);
            two_cp_haps = two_cp_haps
                .iter()
                .filter(|name| name.contains("_strchap"))
                .cloned()
                .collect::<Vec<_>>();
        }
        for hap in &two_cp_haps {
            match confident_identity_from_haplotype_name(hap) {
                Some(StrcHaplotypeIdentity::Strc) => counter_gene += 1,
                Some(StrcHaplotypeIdentity::Strcp1) => counter_pseudo += 1,
                Some(StrcHaplotypeIdentity::Discordant | StrcHaplotypeIdentity::Unresolved)
                | None => {}
            }
        }

        let total_cn = assembled_haps.len() + two_cp_haps.len();
        call.total_cn = Some(total_cn as i32);
        let mut gene_cn = (counter_unknown == 0).then_some(counter_gene);
        call.two_copy_haplotypes = two_cp_haps;
        // check depth between STRC and pseudogene
        if counter_unknown == 0 {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classifies_matching_strc_markers() {
        assert_eq!(
            classify_strc_haplotype(b"1", Some(0), StrcDeletionEvidence::Absent),
            StrcHaplotypeIdentity::Strc
        );
    }

    #[test]
    fn classifies_matching_strcp1_markers() {
        assert_eq!(
            classify_strc_haplotype(b"2", Some(0), StrcDeletionEvidence::Present),
            StrcHaplotypeIdentity::Strcp1
        );
    }

    #[test]
    fn flags_discordant_marker_pairs() {
        for (haplotype, deletion_evidence) in [
            (b"1".as_slice(), StrcDeletionEvidence::Present),
            (b"2".as_slice(), StrcDeletionEvidence::Absent),
        ] {
            assert_eq!(
                classify_strc_haplotype(haplotype, Some(0), deletion_evidence),
                StrcHaplotypeIdentity::Discordant,
                "haplotype {} should be discordant",
                String::from_utf8_lossy(haplotype)
            );
        }
    }

    #[test]
    fn reads_synthetic_deletion_marker_strictly() {
        assert_eq!(
            strc_deletion_evidence(b"11", Some(1), &[]),
            StrcDeletionEvidence::Absent
        );
        assert_eq!(
            strc_deletion_evidence(b"13", Some(1), &[]),
            StrcDeletionEvidence::Present
        );
        assert_eq!(
            strc_deletion_evidence(b"12", Some(1), &[]),
            StrcDeletionEvidence::Unresolved
        );
    }

    #[test]
    fn reads_deletion_encoded_across_internal_sites() {
        assert_eq!(
            strc_deletion_evidence(b"133", None, &[1, 2]),
            StrcDeletionEvidence::Present
        );
        assert_eq!(
            strc_deletion_evidence(b"112", None, &[1, 2]),
            StrcDeletionEvidence::Absent
        );
        assert_eq!(
            strc_deletion_evidence(b"131", None, &[1, 2]),
            StrcDeletionEvidence::Unresolved
        );
    }

    #[test]
    fn leaves_incomplete_or_unrecognized_markers_unresolved() {
        let cases = [
            (b"x".as_slice(), Some(0), StrcDeletionEvidence::Absent),
            (b"0".as_slice(), Some(0), StrcDeletionEvidence::Absent),
            (b"3".as_slice(), Some(0), StrcDeletionEvidence::Absent),
            (b"1".as_slice(), None, StrcDeletionEvidence::Absent),
            (b"1".as_slice(), Some(0), StrcDeletionEvidence::Unresolved),
        ];

        for (haplotype, pivot_index, deletion_evidence) in cases {
            assert_eq!(
                classify_strc_haplotype(haplotype, pivot_index, deletion_evidence),
                StrcHaplotypeIdentity::Unresolved,
                "haplotype {} with pivot {pivot_index:?} and deletion evidence {deletion_evidence:?}",
                String::from_utf8_lossy(haplotype)
            );
        }
    }

    #[test]
    fn duplicate_copy_counting_ignores_unknown_names() {
        assert_eq!(
            confident_identity_from_haplotype_name("strc_strchap1"),
            Some(StrcHaplotypeIdentity::Strc)
        );
        assert_eq!(
            confident_identity_from_haplotype_name("strc_strcp1hap1"),
            Some(StrcHaplotypeIdentity::Strcp1)
        );
        assert_eq!(
            confident_identity_from_haplotype_name("strc_unknownhap1"),
            None
        );
    }
}
