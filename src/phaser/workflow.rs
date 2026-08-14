use crate::assembly::assembly_result::{AssembledPaths, AssemblyResult};
use crate::assembly::variant_graph;
use crate::io::json::{GeneCall, ReadFingerprintMap};
use crate::phaser::{self, Exception, Phaser};
use crate::phaser::{HapInfo, HapInfoForJson, PhasedResult};
use crate::realign::align_mm2_intrinsic;
use crate::toolkit::low_complexity::LowConfidenceSites;
use crate::toolkit::range;
use crate::toolkit::site_selection::CandidateSite;
use crate::toolkit::util::{DError, DResult};

use itertools::Itertools;
use rust_htslib::bam;
use vstr::{VStr, VString};

use std::collections::BTreeMap;

/// Convert a read fingerprint map into a form that can be displayed easily with `serde`.
#[must_use]
pub fn to_string_map(map: &ReadFingerprintMap) -> BTreeMap<String, String> {
    map.clone()
        .into_iter()
        .map(|(id, path)| (id.unique_name(), path.to_string()))
        .collect::<_>()
}

impl Phaser {
    /// Align reads against local reference.
    pub fn align(&mut self) -> DResult {
        let local_bam = self.realigned_bam_path();
        log::debug!("Preparing locus-specific reference for alignment");
        let (local_ref, _was_new) = self.local_reference()?;
        let regions_to_extract = self
            .locus_config()
            .extract_regions(self.settings.genome == "37");
        let regions_to_extract_view = regions_to_extract
            .iter()
            .map(|x| &x[..])
            .collect::<Vec<_>>();
        log::debug!(
            "Aligning reads to the locus-specific reference: bam={local_bam:?}, ref={local_ref:?}"
        );

        let opts = (
            1,
            self.locus_config().chain_bandwidth(),
            self.realign_settings,
            self.try_offset()?,
        );
        align_mm2_intrinsic(
            &self.genome_bam_path(),
            &local_bam,
            &self.settings.genome_reference,
            &local_ref,
            &regions_to_extract_view,
            opts,
        )?;

        bam::index::build(&local_bam, None, bam::index::Type::Bai, 1)?;
        if self
            .locus_config()
            .gene2_region(self.settings.genome == "37")
            .is_some()
        {
            let (secondary_ref, _was_new) = self.secondary_reference()?;
            let opts = (
                1,
                self.locus_config().chain_bandwidth(),
                self.realign_settings,
                self.secondary_offset().ok_or_else(|| {
                    Exception::new(
                        "Required: secondary offset for secondary reference region alignment",
                    )
                })?,
            );
            let secondary_bam = self.realigned_gene2_bam_path();
            align_mm2_intrinsic(
                &self.genome_bam_path(),
                &secondary_bam,
                &self.settings.genome_reference,
                &secondary_ref,
                &regions_to_extract_view,
                opts,
            )?;
            bam::index::build(&secondary_bam, None, bam::index::Type::Bai, 1)?;
        }
        Ok(())
    }

    /// Build/fetch the local reference, run intrinsic realignment, and return
    /// the uppercase local reference sequence used for downstream pileups.
    pub fn realign(&mut self) -> Result<Vec<u8>, DError> {
        log::debug!("Loading local reference sequence for downstream pileups");
        let local_chr = self.local_chr().ok_or(Exception::new(format!(
            "Missing target name for region {}",
            self.realign_region
        )))?;
        let faidx = self.make_local_faidx()?;
        log::debug!(
            "Made local faidx from file at {:?}. Querying with {local_chr}:0-{}",
            self.local_reference()?,
            i32::MAX
        );
        let seq = faidx
            .fetch_seq(&local_chr, 0, i32::MAX as usize)
            .map_err(|e| {
                Exception::new(format!(
                    "Failed to query region {local_chr} for faidx at {:?}. Error: {e:?}",
                    self.local_reference_path()
                ))
            })?
            .to_ascii_uppercase();
        log::debug!("Loaded local reference sequence");
        log::debug!(
            "Reference seq of size {}. Sampled first 100bp: {}",
            seq.len(),
            VStr::from(&seq[..std::cmp::min(seq.len(), 100)])
        );

        log::debug!("Starting read alignment on the locus-specific reference");
        self.align()?;
        assert!(
            self.realigned_bam_path().exists(),
            "{:?} was not created as expected.",
            self.realigned_bam_path()
        );
        Ok(seq)
    }

    /// Discover candidate phasing/homozygous sites after read alignment to the local reference.
    ///
    /// Returns `(hom_sites_to_add, add_sites)` used by downstream read abstraction.
    pub fn get_sites(
        &mut self,
        seq: &[u8],
        min_no_var_region_size: Option<i64>,
        min_vaf: Option<f64>,
    ) -> Result<(Vec<CandidateSite>, Vec<CandidateSite>), DError> {
        let offset = self.try_offset()?;
        log::debug!("Coverage gate passed; computing low-complexity masks");
        self.low_complexity_sites =
            LowConfidenceSites::new(seq, offset, self.settings.homopolymer_window_size);
        self.parse_deletions_from_config()?;
        log::debug!("Discovering candidate large deletions");
        if self.del_data.is_empty() {
            self.discover_big_dels()?;
        }
        self.label_big_dels()?;

        let regions_to_check = self
            .del_data
            .iter()
            .flat_map(|x| {
                let mut ret = vec![];
                if !x.del_reads_partial.is_empty() {
                    ret.push(x.threep());
                    ret.push(x.fivep());
                }
                ret.into_iter()
            })
            .collect::<Vec<range::I64>>();
        log::debug!("Regions to check for deletions: {regions_to_check:?}");
        let (_filtered_variants, _raw_variant_counts) =
            self.get_candidate_pos(&regions_to_check, seq, min_vaf)?;

        self.remove_noisy_sites();
        self.init_het_sites = self.het_sites.clone();
        let hom_sites_to_add = self.add_hom_sites(min_no_var_region_size, None, seq);
        self.remove_noisy_sites();

        log::debug!(
            "Homozygous-site augmentation candidates identified: count={}",
            hom_sites_to_add.len()
        );
        let mut add_sites = self.add_sites.clone();
        if let Some(pivot_site) = self.pivot_site_0based() {
            let het_sites_all_pos = self.het_sites.iter().map(|x| x.pos).collect::<Vec<_>>();
            let add_sites_all_pos = add_sites.iter().map(|x| x.pos).collect::<Vec<_>>();
            if !het_sites_all_pos.contains(&pivot_site) && !add_sites_all_pos.contains(&pivot_site)
            {
                let pos_on_ref = pivot_site - offset;
                let ref_base_u8 = seq[pos_on_ref as usize];
                let non_ref_bases = [b'A', b'C', b'G', b'T']
                    .iter()
                    .filter(|x| **x != ref_base_u8)
                    .copied()
                    .collect::<Vec<_>>();
                let var_base_u8 = non_ref_bases.first().copied().unwrap_or(b'N');
                let ref_base = char::from(ref_base_u8).to_string();
                let var_base = char::from(var_base_u8).to_string();
                let new_variant = CandidateSite::new(pivot_site, ref_base, var_base);
                add_sites.push(new_variant);
            }
        }

        log::debug!(
            "Adding configured/derived sites to fingerprint set: count={}",
            add_sites.len()
        );

        log::debug!(
            "Assigning read fingerprints across {} heterozygous sites",
            self.het_sites.len()
        );
        Ok((hom_sites_to_add, add_sites))
    }

    /// Apply deletion-aware read-fingerprint updates and phase resulting haplotypes.
    ///
    /// Returns `(phase_results, known_deletions_by_code)`.
    pub fn update_indel_and_phase(
        &mut self,
        init_read_hap_map: BTreeMap<crate::io::json::ReadAlignmentId, VString>,
        call: &mut GeneCall,
    ) -> Result<(PhasedResult, BTreeMap<char, String>), DError> {
        log::debug!(
            "Update read abstractions with deletions and phase haplotypes: initial_read_hap_count={}",
            init_read_hap_map.len()
        );
        let mut read_hap_map = init_read_hap_map;
        let known_del = self.update_for_deletions(&mut read_hap_map)?;
        call.read_details = to_string_map(&read_hap_map);
        let phase_results = self.phase_haps(&read_hap_map)?;
        Ok((phase_results, known_del))
    }

    /// Take read fingerprint maps and generate phased haplotypes from them.
    pub(crate) fn phase_haps(&self, reads: &ReadFingerprintMap) -> Result<PhasedResult, DError> {
        let mut min_hap_support = self.settings.min_hap_support as f32;
        if self.settings.targeted {
            let total_depth = self.region_avg_depth[0].0;
            min_hap_support =
                min_hap_support.max(total_depth * self.settings.min_haplotype_frequency as f32);
        }

        let het_sites = self.het_sites.clone();
        let (haps_to_reads, _raw_read_haps) = Self::simplify_read_haps(reads);
        assert!(_raw_read_haps
            .iter()
            .all(|(k, v)| reads.get(k) == Some(&VString::from(v))));
        let mut ret = PhasedResult {
            raw_read_haps: reads.clone(),
            ..Default::default()
        };

        let nvar = het_sites.len();
        if nvar == 0 {
            return Ok(ret);
        }
        ret.assemblies = if nvar == 1 {
            let final_haps = AssembledPaths::from_seqs(["1", "2"].into_iter().map(VString::from));
            let main_haps = final_haps.clone();
            AssemblyResult {
                final_haps,
                main_haps,
                highest_cn: 2,
            }
        } else {
            let pivot_index = self.get_pivot_index();
            let mut settings = variant_graph::Settings::from_pivot(pivot_index);
            settings.min_hap_support = min_hap_support.ceil() as i32;
            let mut graph = variant_graph::Graph::new(reads.clone(), Some(settings));
            graph.construct()?
        };
        if ret.assemblies.main_haps.is_empty() {
            return Ok(ret);
        }

        let mut flat_haps = ret.assemblies.main_haps.iter().cloned().collect::<Vec<_>>();
        let mut read_support = Self::get_read_support(reads, &haps_to_reads, &flat_haps[..])?;
        let mut assembled_haps = self.adjust_spurious_haplotypes(&read_support.0, None, None)?;

        flat_haps = assembled_haps.iter().cloned().collect::<Vec<_>>();
        read_support = Self::get_read_support(reads, &haps_to_reads, &flat_haps[..])?;
        let uniquely_supporting_reads = &read_support.0;
        let mut flat_read_counts = uniquely_supporting_reads
            .values()
            .map(std::vec::Vec::len)
            .sorted();
        let min_hap_support = if min_hap_support == 4.0
            && flat_read_counts.len() > 2
            && flat_read_counts.next().is_some_and(|x| x <= 4)
            && flat_read_counts.next().is_some_and(|x| x >= 12)
            && !assembled_haps.iter().any(|hap| hap.contains(&b'x'))
        {
            5.0
        } else {
            min_hap_support
        };
        assembled_haps =
            AssembledPaths::from_seqs(uniquely_supporting_reads.iter().filter_map(|x| {
                if x.1.len() as f32 >= min_hap_support {
                    Some(x.0)
                } else {
                    None
                }
            }));
        flat_haps = assembled_haps.iter().cloned().collect::<Vec<_>>();
        let read_support = Self::get_read_support(reads, &haps_to_reads, &flat_haps[..])?;
        let (ref uniquely_supporting_reads, ref nonuniquely_supporting_reads, ref read_counts) =
            read_support;

        ret.assemblies.main_haps = assembled_haps;
        ret.uniquely_supporting_reads = uniquely_supporting_reads
            .iter()
            .map(|(k, v)| {
                (
                    VString::from(k),
                    v.iter()
                        .map(std::string::ToString::to_string)
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<BTreeMap<VString, Vec<String>>>();
        ret.nonuniquely_supporting_reads = nonuniquely_supporting_reads
            .iter()
            .map(|(name, ids)| {
                (
                    name.to_string(),
                    ids.iter()
                        .map(|id| VString::from(&flat_haps[*id as usize]))
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<BTreeMap<String, Vec<VString>>>();
        ret.read_counts = (
            read_counts
                .0
                .iter()
                .map(|(k, v)| (VString::from(k), *v))
                .collect::<BTreeMap<VString, i32>>(),
            read_counts
                .1
                .iter()
                .map(|(k, v)| (VString::from(k), *v))
                .collect::<BTreeMap<VString, f64>>(),
        );
        Ok(ret)
    }

    fn try_update_indel_and_phase_with_fallback(
        &mut self,
        init_read_hap_map: ReadFingerprintMap,
        call: &mut GeneCall,
    ) -> Option<(PhasedResult, BTreeMap<char, String>)> {
        match self.update_indel_and_phase(init_read_hap_map, call) {
            Ok(results) => Some(results),
            Err(e) => {
                log::debug!(
                    "Failed to update indel and phase with the following error: {:?}",
                    e
                );
                self.populate_call_site_fields(call);
                None
            }
        }
    }

    fn build_assembled_haplotypes<'a, I>(&self, main_haps_clone: I) -> BTreeMap<VStr<'a>, String>
    where
        I: IntoIterator<Item = &'a VString>,
    {
        let mut assembled_haps = BTreeMap::new();
        let mod_gene_name = self.gene_name().split_terminator('-').join(",");
        for (idx, hap) in main_haps_clone.into_iter().enumerate() {
            assembled_haps.insert(hap.vstr(), format!("{mod_gene_name}_hap{}", idx + 1));
        }
        assembled_haps
    }

    fn apply_fusion_step_if_enabled<'a>(
        &mut self,
        assembled_haps: BTreeMap<VStr<'a>, String>,
        call: &mut GeneCall,
    ) -> Result<(BTreeMap<VStr<'a>, String>, bool), DError> {
        if let Some(fusion_direction) = self.config.locus.call_fusion() {
            let (assembled_haps_renamed, two_cp_haps, fusions_called) =
                self.find_fusion(&assembled_haps, fusion_direction.to_string())?;
            call.region_specific_info.insert(
                String::from("fusions_called"),
                serde_json::to_value(&fusions_called)?,
            );
            let total_cn = assembled_haps_renamed.len() + two_cp_haps.len();
            call.total_cn = Some(total_cn as i32);
            call.two_copy_haplotypes = two_cp_haps;
            Ok((assembled_haps_renamed, true))
        } else {
            Ok((assembled_haps, false))
        }
    }

    fn assign_final_haplotypes_from_assembled<'a>(
        &self,
        assembled_haps: &BTreeMap<VStr<'a>, String>,
        call: &mut GeneCall,
    ) {
        call.final_haplotypes = assembled_haps
            .clone()
            .into_iter()
            .map(|(k, v)| (k.to_string(), v))
            .collect::<BTreeMap<_, _>>();
    }

    fn assign_haplotype_details_from_phase<'a>(
        &mut self,
        phase_results: &PhasedResult,
        known_del: &BTreeMap<char, String>,
        assembled_haps: &BTreeMap<VStr<'a>, String>,
        call: &mut GeneCall,
    ) -> Result<BTreeMap<String, HapInfo>, DError> {
        let haps =
            self.output_variants_in_haps(phase_results, known_del, assembled_haps.clone())?;
        call.haplotype_details = haps
            .iter()
            .map(|(key, val)| (key.clone(), HapInfoForJson::from(val)))
            .collect::<BTreeMap<_, _>>();
        Ok(haps)
    }

    fn maybe_adjust_copy_number_from_depth<'a>(
        &mut self,
        assembled_haps: &BTreeMap<VStr<'a>, String>,
        haps: &BTreeMap<String, HapInfo>,
        phase_results: &PhasedResult,
        fusion_call_step: bool,
        call: &mut GeneCall,
    ) -> Result<(), DError> {
        if fusion_call_step {
            return Ok(());
        }
        let (two_cp_haps, total_cn) = self.adjust_depth(
            assembled_haps.clone(),
            haps.clone(),
            phase_results.clone(),
            false,
            true,
            false,
            0.05,
        )?;
        call.two_copy_haplotypes = two_cp_haps;
        call.total_cn = if total_cn <= 1 {
            None
        } else {
            Some(total_cn as i32)
        };
        Ok(())
    }

    fn maybe_phase_alleles_for_call<'a>(
        &mut self,
        phase_results: &mut PhasedResult,
        assembled_haps: &BTreeMap<VStr<'a>, String>,
        call: &mut GeneCall,
    ) -> Result<(), DError> {
        if !self.config.locus.to_phase() {
            return Ok(());
        }
        let allele_result = self.phase_alleles(phase_results, assembled_haps, None);
        let mut alleles = allele_result.alleles;
        let mut raw_alleles = allele_result.raw_alleles;
        if self.all_haps_phased_onto_one_allele(&alleles, assembled_haps)? {
            alleles = vec![];
            raw_alleles = vec![];
        }
        call.region_specific_info
            .insert(String::from("alleles_final"), alleles.into());
        call.region_specific_info.insert(
            String::from("haplotype_links"),
            serde_json::to_value(&allele_result.haplotype_links)?,
        );
        call.region_specific_info
            .insert(String::from("raw_alleles"), raw_alleles.into());
        Ok(())
    }

    /// Run the default phasing pipeline for genes that do not have a dedicated gene-specific runner.
    pub(crate) fn run_default_gene_pipeline(&mut self) -> Result<GeneCall, DError> {
        let seq = self.realign()?;
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
        let tid = self
            .genome_tid()
            .map(|x| x as i32)
            .ok_or_else(|| phaser::Exception::new("run: missing chr tid in BAM header"))?;
        let init_read_hap_map = self.haplotypes_from_reads(
            None,
            &hom_sites_to_add,
            Some(&add_sites),
            None,
            (5, true, Some(50u32)),
            tid,
            None,
            &hom_sites_to_add,
        )?;

        let Some((mut phase_results, known_del)) =
            self.try_update_indel_and_phase_with_fallback(init_read_hap_map.clone(), &mut call)
        else {
            return Ok(call);
        };

        let main_haps_clone = phase_results.assemblies.main_haps.clone();
        let assembled_haps = self.build_assembled_haplotypes(main_haps_clone.iter());
        let (assembled_haps, fusion_call_step) =
            self.apply_fusion_step_if_enabled(assembled_haps, &mut call)?;
        self.assign_final_haplotypes_from_assembled(&assembled_haps, &mut call);
        let haps = self.assign_haplotype_details_from_phase(
            &phase_results,
            &known_del,
            &assembled_haps,
            &mut call,
        )?;
        self.maybe_adjust_copy_number_from_depth(
            &assembled_haps,
            &haps,
            &phase_results,
            fusion_call_step,
            &mut call,
        )?;
        self.maybe_phase_alleles_for_call(&mut phase_results, &assembled_haps, &mut call)?;

        self.fill_in_call(phase_results, &mut call);
        Ok(call)
    }
}
