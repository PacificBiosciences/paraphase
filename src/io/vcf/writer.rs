//! Stateful VCF writing orchestration.
//!
//! This file contains `VcfWriter`, which ties together:
//! - region/gene context from `Phaser`
//! - per-haplotype pileup-derived variant collection
//! - merged record emission (small variants + symbolic SV records)
//!
//! Lower-level parsing/formatting and small pure helpers live in `helpers.rs`.

use crate::io::json::GeneCall;
use crate::phaser::Exception;
use crate::phaser::Phaser;
use crate::toolkit::util::{DError, DResult, CLI_COMMAND, FULL_VERSION};
use rust_htslib::bcf::{self, record::GenotypeAllele, Format};
use rust_htslib::{bam, bam::Read};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::path::PathBuf;

use super::helpers::*;
use super::types::*;

pub struct VcfWriter<'a> {
    phaser: &'a Phaser,
    call: &'a GeneCall,
    write_no_calls: bool,
    gene1only: bool,
}

impl<'a> VcfWriter<'a> {
    #[must_use]
    /// Create a VCF writer facade over one phaser/call result pair.
    pub fn new(
        phaser: &'a Phaser,
        call: &'a GeneCall,
        write_no_calls: bool,
        gene1only: bool,
    ) -> Self {
        Self {
            phaser,
            call,
            write_no_calls,
            gene1only,
        }
    }

    /// Get header for VCF
    /// # Arguments
    /// * `bam_header` - header of the realigned tagged bam
    /// * `chr` - chromosome of the archetype gene
    pub fn get_vcf_header(
        &self,
        bam_reader: &bam::IndexedReader,
        chr: &str,
        sample_ids: &Vec<String>,
        has_sv: bool,
    ) -> bcf::header::Header {
        let mut vcf_header = bcf::header::Header::new();
        // add header
        for line in VCF_LINES.iter() {
            vcf_header.push_record(line.as_bytes());
        }
        if has_sv {
            for line in VCF_LINES_SV.iter() {
                vcf_header.push_record(line.as_bytes());
            }
        }
        vcf_header.push_record(format!("##paraphase_version={}", &**FULL_VERSION).as_bytes());
        vcf_header.push_record(
            format!(
                "##paraphase_command={}",
                serde_json::to_string(&*CLI_COMMAND).unwrap_or_else(|_| "\"<unavailable>\"".into())
            )
            .as_bytes(),
        );
        let alleles = self.call.region_specific_info.get("alleles_final");
        if alleles
            .and_then(|value| value.as_array())
            .is_some_and(|items| !items.is_empty())
        {
            for line in VCF_LINE_ALLELE.iter() {
                vcf_header.push_record(line.as_bytes());
            }
        }
        if let Some(records) = bam::Header::from_template(bam_reader.header())
            .to_hashmap()
            .get("SQ")
        {
            for record in records {
                if record["SN"] == chr {
                    let contig_line =
                        format!(r#"##contig=<ID={},length={}>"#, record["SN"], record["LN"]);
                    vcf_header.push_record(contig_line.as_bytes());
                }
            }
        }
        // write haplotypes as "sample" columns
        for hap in sample_ids {
            vcf_header.push_sample(hap.to_string().as_bytes());
        }
        vcf_header
    }

    /// Get sequence of the reference region
    pub fn get_ref_seq(&self) -> Result<Vec<u8>, DError> {
        log::debug!("Loading primary locus reference sequence for VCF emission");
        let local_chr = self.phaser.local_chr().ok_or(Exception::new(format!(
            "Missing target name for region {}",
            self.phaser.realign_region
        )))?;
        let faidx = self.phaser.make_local_faidx()?;
        let ref_seq = faidx
            .fetch_seq(&local_chr, 0, i32::MAX as usize)
            .map_err(|e| {
                Exception::new(format!(
                    "Failed to fetch reference sequence for region {local_chr} using faidx at {:?}: {e:?}",
                    self.phaser.local_reference_path()
                ))
            })?
            .to_ascii_uppercase();
        Ok(ref_seq)
    }

    /// Get sequence of the gene2 reference region
    pub fn get_ref_seq_gene2(&self) -> Result<Vec<u8>, DError> {
        log::debug!("Loading secondary (gene2) locus reference sequence for VCF emission");
        let (nchr, start, end) = self.phaser.parsed_nchr_secondary_0based().ok_or_else(|| {
            Exception::new(format!(
                "Failed to parse secondary region for gene '{}' from config",
                self.phaser.gene_name()
            ))
        })?;
        let local_chr = format!("{}_{}_{}", nchr, start + 1, end + 1);
        let faidx = self.phaser.make_local_faidx_gene2()?;
        let ref_seq = faidx
            .fetch_seq(&local_chr, 0, i32::MAX as usize)?
            .to_ascii_uppercase();
        Ok(ref_seq)
    }

    /// Resolve haplotype output boundaries for VCF annotation fields.
    ///
    /// When `match_range` is true, boundaries are projected to the secondary region.
    pub fn get_hap_bound(
        &self,
        hap_name: &str,
        match_range: bool,
    ) -> Result<HapBoundForVcf, DError> {
        let is_truncated = &self
            .call
            .haplotype_details
            .get(hap_name)
            .ok_or_else(|| {
                Exception::new(format!(
                    "Haplotype '{hap_name}' missing from haplotype_details"
                ))
            })?
            .is_truncated;
        // note this is 1-based
        let hap_bound = self
            .call
            .haplotype_details
            .get(hap_name)
            .ok_or_else(|| {
                Exception::new(format!(
                    "Haplotype '{hap_name}' missing from haplotype_details"
                ))
            })?
            .boundary;
        let sites = &self.call.sites_for_phasing;
        let start = hap_bound[0] - 1;
        let mut start_strict = start;
        if start > self.phaser.left_boundary_0based() && !is_truncated.contains(&String::from("5p"))
        {
            for var in sites {
                let pos = var
                    .split("_")
                    .map(std::borrow::ToOwned::to_owned)
                    .collect::<Vec<_>>()
                    .first()
                    .ok_or_else(|| {
                        Exception::new(format!(
                            "Malformed variant '{}' in sites_for_phasing; expected 'pos_ref_alt'",
                            var
                        ))
                    })?
                    .parse::<i64>()?;
                if pos > start {
                    start_strict = pos;
                    break;
                }
            }
        }
        let end: i64 = hap_bound[1] - 1;
        let mut end_strict = end;
        if end < self.phaser.right_boundary_0based() && !is_truncated.contains(&String::from("3p"))
        {
            for var in sites.iter().rev() {
                let pos = var
                    .split("_")
                    .map(std::borrow::ToOwned::to_owned)
                    .collect::<Vec<_>>()
                    .first()
                    .ok_or_else(|| {
                        Exception::new(format!(
                            "Malformed variant '{}' in sites_for_phasing; expected 'pos_ref_alt'",
                            var
                        ))
                    })?
                    .parse::<i64>()?;
                if pos < end {
                    end_strict = pos;
                    break;
                }
            }
        }
        if match_range {
            // 0-based
            let start_gene2_tmp = self.phaser.get_range_in_other_gene(start, None);
            let start_strict_gene2_tmp = self.phaser.get_range_in_other_gene(start_strict, None);
            let end_gene2_tmp = self.phaser.get_range_in_other_gene(end, None);
            let end_strict_gene2_tmp = self.phaser.get_range_in_other_gene(end_strict, None);

            let (_nchr, gene2_region_start, gene2_region_end) =
                self.phaser.parsed_nchr_secondary_0based().ok_or_else(|| {
                    Exception::new(format!(
                        "Failed to parse secondary region for gene '{}' from config",
                        self.phaser.gene_name()
                    ))
                })?;
            let start_gene2: i64;
            let end_gene2: i64;
            let start_strict_gene2: i64;
            let end_strict_gene2: i64;
            if let (Some(start_gene2_tmp), Some(end_gene2_tmp)) = (start_gene2_tmp, end_gene2_tmp) {
                start_gene2 = start_gene2_tmp.min(end_gene2_tmp);
                end_gene2 = start_gene2_tmp.max(end_gene2_tmp);
            } else {
                start_gene2 = gene2_region_start;
                end_gene2 = gene2_region_end;
            }
            if let (Some(start_strict_gene2_tmp), Some(end_strict_gene2_tmp)) =
                (start_strict_gene2_tmp, end_strict_gene2_tmp)
            {
                start_strict_gene2 = start_strict_gene2_tmp.min(end_strict_gene2_tmp);
                end_strict_gene2 = start_strict_gene2_tmp.max(end_strict_gene2_tmp);
            } else {
                start_strict_gene2 = gene2_region_start;
                end_strict_gene2 = gene2_region_end;
            }
            return Ok(HapBoundForVcf {
                hap_name: hap_name.to_string(),
                start: start_gene2,
                start_strict: start_strict_gene2,
                end: end_gene2,
                end_strict: end_strict_gene2,
                is_truncated: vec![],
            });
        }
        Ok(HapBoundForVcf {
            hap_name: hap_name.to_string(),
            start,
            start_strict,
            end,
            end_strict,
            is_truncated: is_truncated.to_vec(),
        })
    }

    /// Collect Python-compatible symbolic SV calls that should be emitted in the VCF.
    ///
    /// These records are synthesized from gene-specific call metadata instead of pileup
    /// consensus and are emitted for the same genes/conditions as the Python implementation.
    pub(super) fn get_special_variants(
        &self,
        gene2: bool,
    ) -> Result<BTreeMap<String, String>, DError> {
        let mut special_variants = BTreeMap::new();
        let gene_name = self.phaser.gene_name();

        if gene_name == "ikbkg" {
            if let Some(del_haps) = self
                .call
                .region_specific_info
                .get("deletion_haplotypes")
                .and_then(|value| value.as_array())
            {
                let deletion_name_gene1 = self
                    .phaser
                    .locus_config()
                    .get("deletion1_in_gene1")
                    .and_then(|value| value.as_str())
                    .ok_or_else(|| Exception::new("Missing deletion1_in_gene1 in IKBKG config"))?;
                let deletion_name_gene2 = self
                    .phaser
                    .locus_config()
                    .get("deletion1_in_gene2")
                    .and_then(|value| value.as_str())
                    .ok_or_else(|| Exception::new("Missing deletion1_in_gene2 in IKBKG config"))?;
                for del_hap in del_haps.iter().filter_map(|value| value.as_str()) {
                    let del_name = if del_hap.contains("pseudo") && gene2 {
                        deletion_name_gene2
                    } else {
                        deletion_name_gene1
                    };
                    special_variants.insert(del_hap.to_string(), del_name.to_string());
                }
            }
        } else if gene_name == "f8" {
            if let Some(sv_called) = self
                .call
                .region_specific_info
                .get("sv_called")
                .and_then(|value| value.as_object())
            {
                if !sv_called.is_empty() {
                    let extract_regions = self
                        .phaser
                        .locus_config()
                        .extract_regions(self.phaser.settings.genome == "37");
                    let extract_region1 = extract_regions
                        .first()
                        .ok_or_else(|| Exception::new("F8 config is missing extract region 1"))?;
                    let extract_region2 = extract_regions
                        .get(1)
                        .ok_or_else(|| Exception::new("F8 config is missing extract region 2"))?;
                    let extract_region3 = extract_regions
                        .get(2)
                        .ok_or_else(|| Exception::new("F8 config is missing extract region 3"))?;
                    let extract_region1_end = extract_region1
                        .split('-')
                        .next_back()
                        .ok_or_else(|| Exception::new("Malformed F8 extract region 1"))?;
                    let extract_region2_start = extract_region2
                        .split(':')
                        .nth(1)
                        .and_then(|part| part.split('-').next())
                        .ok_or_else(|| Exception::new("Malformed F8 extract region 2"))?;
                    let extract_region3_start = extract_region3
                        .split(':')
                        .nth(1)
                        .and_then(|part| part.split('-').next())
                        .ok_or_else(|| Exception::new("Malformed F8 extract region 3"))?;
                    for (sv_hap, sv_value) in sv_called {
                        if let Some(sv_type) = sv_value.as_str() {
                            let sv_name = match sv_type {
                                "inversion" => {
                                    format!("{extract_region1_end}_INV_{extract_region3_start}")
                                }
                                "deletion" => {
                                    format!("{extract_region1_end}_DEL_{extract_region2_start}")
                                }
                                _ => continue,
                            };
                            special_variants.insert(sv_hap.to_string(), sv_name);
                        }
                    }
                }
            }
        }

        Ok(special_variants)
    }

    /// Emit one symbolic SV VCF record with already-prepared sample fields.
    fn write_sv_record(
        &self,
        writer: &mut bcf::Writer,
        start_1based: i64,
        sv_type: &str,
        end_1based: i64,
        gts: &[GenotypeAllele],
        dps: &[Vec<u8>],
        ads: &[Vec<u8>],
        hap_bounds: &str,
        allele_info: &str,
    ) -> DResult {
        let mut record = writer.empty_record();
        let contig = self
            .phaser
            .chr()
            .ok_or_else(|| {
                Exception::new(format!(
                    "Chromosome not found for gene '{}' while writing merged VCF records",
                    self.phaser.gene_name()
                ))
            })?
            .as_bytes();
        let rid = writer.header().name2rid(contig)?;
        record.set_rid(Some(rid));
        record.set_pos(start_1based - 1);
        record.set_qual(f32::from_bits(0x7F800001));
        record.set_alleles(&[b"N".as_slice(), format!("<{sv_type}>").as_bytes()])?;
        if gts.contains(&GenotypeAllele::Phased(1)) {
            record.set_filters(&["PASS".as_bytes()])?;
        } else {
            record.set_filters(&["LowQual".as_bytes()])?;
        }
        record.push_genotypes(gts)?;
        record.push_format_string(b"DP", dps)?;
        record.push_format_string(b"AD", ads)?;
        record.push_info_string(b"SVTYPE", &[sv_type.as_bytes()])?;
        record.push_info_integer(b"END", &[end_1based as i32])?;
        record.push_info_integer(b"SVLEN", &[(end_1based - start_1based) as i32])?;
        record.push_info_string(b"HPBOUND", &[hap_bounds.as_bytes()])?;
        if !allele_info.is_empty() {
            record.push_info_string(b"ALLELE", &[allele_info.as_bytes()])?;
        }
        writer.write(&record)?;
        Ok(())
    }

    /// Emit one small-variant VCF record with already-prepared sample fields.
    fn write_small_variant_record(
        &self,
        writer: &mut bcf::Writer,
        pos: i64,
        ref_base: &str,
        alt_base: &str,
        gts: &[GenotypeAllele],
        dps: &[Vec<u8>],
        ads: &[Vec<u8>],
        hap_bounds: &str,
        allele_info: &str,
    ) -> DResult {
        let mut record = writer.empty_record();
        let contig = self
            .phaser
            .chr()
            .ok_or_else(|| {
                Exception::new(format!(
                    "Chromosome not found for gene '{}' while writing merged VCF records",
                    self.phaser.gene_name()
                ))
            })?
            .as_bytes();
        let rid = writer.header().name2rid(contig)?;
        record.set_rid(Some(rid));
        record.set_pos(pos);
        record.set_qual(f32::from_bits(0x7F800001));
        let allele1 = &[ref_base.as_bytes(), alt_base.as_bytes()];
        let allele2 = &[ref_base.as_bytes(), ".".as_bytes()];
        let alleles: &[&[u8]] = if alt_base != ref_base {
            allele1
        } else {
            allele2
        };
        record.set_alleles(alleles)?;
        if gts.contains(&GenotypeAllele::Phased(1)) {
            record.set_filters(&["PASS".as_bytes()])?;
        } else {
            record.set_filters(&["LowQual".as_bytes()])?;
        }
        record.push_genotypes(gts)?;
        record.push_format_string(b"DP", dps)?;
        record.push_format_string(b"AD", ads)?;
        record.push_info_string(b"HPBOUND", &[hap_bounds.as_bytes()])?;
        if !allele_info.is_empty() {
            record.push_info_string(b"ALLELE", &[allele_info.as_bytes()])?;
        }
        writer.write(&record)?;
        Ok(())
    }

    /// Merge variant calls from all haplotypes and write VCF entries.
    ///
    /// Handles both small variants and symbolic SV records with Python-compatible
    /// sample-field padding across gene1/gene2 column groups.
    pub fn merge_vcf(
        &self,
        mut writer: bcf::Writer,
        hap_variant_info: Vec<HapVariantInfo>,
    ) -> DResult {
        let mut haps_ids = Vec::new();
        let mut haps_ids1 = Vec::new();
        let mut haps_ids2 = Vec::new();
        // get hap bound info field
        let mut hap_bounds = String::new();
        for (counter, this_hap_variant_info) in hap_variant_info.iter().enumerate() {
            for hap in &this_hap_variant_info.hap_info {
                haps_ids.push(hap.hap_name.clone());
                if counter == 0 {
                    haps_ids1.push(hap.hap_name.clone());
                } else {
                    haps_ids2.push(hap.hap_name.clone());
                }
                let hap_bound = format_hap_bound(hap.start + 1, hap.end + 1, &hap.is_truncated);
                if !hap_bounds.is_empty() {
                    hap_bounds += ",";
                }
                hap_bounds += &hap_bound;
            }
        }

        // get allele info field
        let mut allele_info = String::new();
        let alleles = self.call.region_specific_info.get("alleles_final");
        if let Some(raw_alleles) = alleles.and_then(|x| x.as_array()) {
            let mut parsed_alleles = Vec::new();
            for allele in raw_alleles {
                if let Some(allele_parsed) = allele.as_array() {
                    let allele_string = allele_parsed
                        .iter()
                        .filter_map(|x| x.as_str().map(str::to_string))
                        .collect::<Vec<_>>();
                    if !allele_string.is_empty() {
                        parsed_alleles.push(allele_string);
                    }
                }
            }
            for allele in parsed_alleles {
                let allele_joined = allele.join("+");
                if !allele_info.is_empty() {
                    allele_info += ",";
                }
                allele_info += &allele_joined;
            }
        }

        for (counter, this_hap_variant_info) in hap_variant_info.iter().enumerate() {
            let hap_info = &this_hap_variant_info.hap_info;
            let variants_info = &this_hap_variant_info.variants_info;
            let sv_variants = &this_hap_variant_info.sv_variants;
            let all_positions: BTreeSet<i64> = sv_variants
                .keys()
                .copied()
                .chain(variants_info.keys().copied())
                .collect();
            for pos in all_positions {
                if let Some(pos_calls) = sv_variants.get(&pos) {
                    let mut variant_observed = HashSet::new();
                    for var_name in pos_calls.iter().flatten() {
                        variant_observed.insert(var_name.clone());
                    }
                    for variant_name in variant_observed {
                        let (start_1based, sv_type, end_1based) = parse_symbolic_sv(&variant_name)?;
                        assert_eq!(pos_calls.len(), hap_info.len());
                        let (valid_gts, mut gts, mut dps, mut ads) = sv_sample_fields(
                            &variant_name,
                            start_1based,
                            pos_calls,
                            hap_info,
                            variants_info,
                        );
                        let mut write_variant = false;
                        if self.write_no_calls {
                            if gts.contains(&GenotypeAllele::Phased(1))
                                || valid_gts.contains(&GenotypeAllele::UnphasedMissing)
                            {
                                write_variant = true;
                            }
                        } else if gts.contains(&GenotypeAllele::Phased(1)) {
                            write_variant = true;
                        }
                        if write_variant {
                            pad_sample_formats_for_counter(
                                counter, &haps_ids, &haps_ids1, &haps_ids2, &mut gts, &mut dps,
                                &mut ads,
                            );
                            self.write_sv_record(
                                &mut writer,
                                start_1based,
                                &sv_type,
                                end_1based,
                                &gts,
                                &dps,
                                &ads,
                                &hap_bounds,
                                &allele_info,
                            )?;
                        }
                    }
                }
                if let Some(pos_calls) = variants_info.get(&pos) {
                    let mut variant_observed = HashSet::new();
                    for pos_call_value in pos_calls.iter().flatten() {
                        if let Some(ref variant_call) = pos_call_value.base {
                            let ref_base = pos_call_value.ref_base.clone();
                            let alt_base = variant_call.to_string();
                            variant_observed.insert((ref_base, alt_base));
                        } else {
                            let ref_base = pos_call_value.ref_base.clone();
                            variant_observed.insert((ref_base.clone(), ref_base));
                        }
                    }
                    let ref_only = test_ref_only(&variant_observed)?;
                    //log::debug!("pos {pos}, variant_observed {variant_observed:?}, ref_only {ref_only}");
                    for (ref_base, alt_base) in variant_observed {
                        assert_eq!(pos_calls.len(), hap_info.len());
                        let (valid_gts, mut gts, mut dps, mut ads) = small_variant_sample_fields(
                            pos, &ref_base, &alt_base, pos_calls, hap_info,
                        );
                        let mut write_variant = false;
                        if self.write_no_calls {
                            if (ref_base != alt_base || ref_only)
                                && alt_base != "*"
                                && (gts.contains(&GenotypeAllele::Phased(1))
                                    || valid_gts.contains(&GenotypeAllele::UnphasedMissing))
                            {
                                write_variant = true;
                            }
                        } else if ref_base != alt_base
                            && alt_base != "*"
                            && gts.contains(&GenotypeAllele::Phased(1))
                        {
                            write_variant = true;
                        }
                        if write_variant {
                            pad_sample_formats_for_counter(
                                counter, &haps_ids, &haps_ids1, &haps_ids2, &mut gts, &mut dps,
                                &mut ads,
                            );
                            self.write_small_variant_record(
                                &mut writer,
                                pos,
                                &ref_base,
                                &alt_base,
                                &gts,
                                &dps,
                                &ads,
                                &hap_bounds,
                                &allele_info,
                            )?;
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// Write one per-gene VCF file.
    ///
    /// For two-reference genes, emits merged columns across gene1/gene2 hap sets.
    pub fn write_vcf(&self, vcf_dir: &PathBuf) -> DResult {
        if self.call.failed_for_coverage {
            return Ok(());
        }

        let suffix = format!(
            "{}_{}.vcf",
            self.phaser.settings.sample_id,
            self.phaser.gene_name()
        );
        let chr = self.phaser.chr().ok_or_else(|| {
            Exception::new(format!(
                "Chromosome not found for gene '{}' while preparing VCF output",
                self.phaser.gene_name()
            ))
        })?;
        let vcf_file = vcf_dir.join(suffix);
        let gene1_output_bam = self.phaser.realigned_tagged_bam_path();
        let bam_reader = bam::IndexedReader::from_path(gene1_output_bam)?;

        let hap_variant_info: Vec<HapVariantInfo> = if !self
            .phaser
            .gene_config()
            .two_reference_regions_genes
            .contains(self.phaser.gene_name())
            || self.gene1only
        {
            vec![self.get_variants_for_vcf(&self.call.final_haplotypes, false, false)?]
        } else {
            let (gene1_haps, gene2_haps) = self.separate_two_genes();
            let hap_variant_info_gene1 = self.get_variants_for_vcf(&gene1_haps, true, false)?;
            let hap_variant_info_gene2 = self.get_variants_for_vcf(&gene2_haps, true, true)?;
            let hap_variant_info_gene1_pos = hap_variant_info_gene1
                .variants_info
                .first_key_value()
                .map(|x| *x.0);
            let hap_variant_info_gene2_pos = hap_variant_info_gene2
                .variants_info
                .first_key_value()
                .map(|x| *x.0);
            match (hap_variant_info_gene1_pos, hap_variant_info_gene2_pos) {
                (Some(gene1_pos), Some(gene2_pos)) => {
                    if gene1_pos < gene2_pos {
                        vec![hap_variant_info_gene1, hap_variant_info_gene2]
                    } else {
                        vec![hap_variant_info_gene2, hap_variant_info_gene1]
                    }
                }
                _ => vec![hap_variant_info_gene1, hap_variant_info_gene2],
            }
        };
        let has_sv = hap_variant_info
            .iter()
            .any(|entry| !entry.sv_variants.is_empty());

        let mut haplotype_ids = Vec::new();
        for this_entry in &hap_variant_info {
            for hap in &this_entry.hap_info {
                haplotype_ids.push(hap.hap_name.clone());
            }
        }
        if !haplotype_ids.is_empty() {
            let vcf_header = self.get_vcf_header(&bam_reader, chr, &haplotype_ids, has_sv);
            let writer = bcf::Writer::from_path(&vcf_file, &vcf_header, true, Format::Vcf)
                .map_err(|_| format!("Invalid VCF output path: {}", &vcf_file.display()))?;
            self.merge_vcf(writer, hap_variant_info)?;
        }
        Ok(())
    }

    /// Build ordered haplotype boundary metadata for VCF columns, including
    /// two-copy expansion and the homozygous fallback shape.
    fn build_hap_info(
        &self,
        final_haplotypes: &BTreeMap<String, String>,
        two_cp_haplotypes: &[String],
        is_gene2: bool,
        match_range: bool,
    ) -> Result<Vec<HapBoundForVcf>, DError> {
        let mut hap_info = Vec::new();
        for hap_name in final_haplotypes.values() {
            let hap_boundaries = self.get_hap_bound(hap_name, match_range)?;
            hap_info.push(hap_boundaries.clone());
            if two_cp_haplotypes.contains(hap_name) {
                let cp2_boundaries = HapBoundForVcf {
                    hap_name: format!("{hap_name}_cp2"),
                    start: hap_boundaries.start,
                    start_strict: hap_boundaries.start_strict,
                    end: hap_boundaries.end,
                    end_strict: hap_boundaries.end_strict,
                    is_truncated: hap_boundaries.is_truncated,
                };
                hap_info.push(cp2_boundaries);
            }
        }
        if !is_gene2 && self.call.heterozygous_sites.is_empty() && final_haplotypes.is_empty() {
            let homozygous_boundaries = HapBoundForVcf {
                hap_name: format!("{}_homozygous_hap1", self.phaser.gene_name()),
                start: self.phaser.left_boundary(),
                start_strict: self.phaser.left_boundary(),
                end: self.phaser.right_boundary(),
                end_strict: self.phaser.right_boundary(),
                is_truncated: vec![],
            };
            hap_info.push(homozygous_boundaries);
            let homozygous_boundaries = HapBoundForVcf {
                hap_name: format!("{}_homozygous_hap1_cp2", self.phaser.gene_name()),
                start: self.phaser.left_boundary(),
                start_strict: self.phaser.left_boundary(),
                end: self.phaser.right_boundary(),
                end_strict: self.phaser.right_boundary(),
                is_truncated: vec![],
            };
            hap_info.push(homozygous_boundaries);
        }
        Ok(hap_info)
    }

    /// Prepare variants for each haplotype. Consider gene1/gene2 scenarios.
    pub(super) fn get_variants_for_vcf(
        &self,
        final_haplotypes: &BTreeMap<String, String>,
        is_gene2: bool,
        match_range: bool,
    ) -> Result<HapVariantInfo, DError> {
        // Six variables below vary based on gene2 or gene1.
        // read tagged bam
        let mut bam_reader = if !is_gene2 || !match_range {
            bam::IndexedReader::from_path(self.phaser.realigned_tagged_bam_path())?
        } else {
            bam::IndexedReader::from_path(self.phaser.realigned_tagged_gene2_bam_path())?
        };
        let ref_seq = if !is_gene2 || !match_range {
            self.get_ref_seq()?
        } else {
            self.get_ref_seq_gene2()?
        };
        let ref_seq_len = ref_seq.len();
        let offset = if !is_gene2 || !match_range {
            self.phaser.offset()
        } else {
            self.phaser.secondary_offset().ok_or_else(|| {
                Exception::new(format!(
                    "Secondary offset not found for gene '{}' in gene2 VCF mode",
                    self.phaser.gene_name()
                ))
            })?
        };
        let chr = if !is_gene2 || !match_range {
            self.phaser.chr().ok_or_else(|| {
                Exception::new(format!(
                    "Chromosome not found for gene '{}' while extracting variant pileups",
                    self.phaser.gene_name()
                ))
            })?
        } else {
            self.phaser
                .parsed_nchr_secondary_0based()
                .ok_or_else(|| {
                    Exception::new(format!(
                        "Gene2 region not found for gene '{}' while extracting variant pileups",
                        self.phaser.gene_name()
                    ))
                })?
                .0
        };
        let region_start = if !is_gene2 || !match_range {
            self.phaser.left_boundary_0based()
        } else {
            self.phaser
                .parsed_nchr_secondary_0based()
                .ok_or_else(|| {
                    Exception::new(format!(
                        "Gene2 region not found for gene '{}' while extracting region start",
                        self.phaser.gene_name()
                    ))
                })?
                .1
        };
        let region_end = if !is_gene2 || !match_range {
            self.phaser.right_boundary_0based()
        } else {
            self.phaser
                .parsed_nchr_secondary_0based()
                .ok_or_else(|| {
                    Exception::new(format!(
                        "Gene2 region not found for gene '{}' while extracting region end",
                        self.phaser.gene_name()
                    ))
                })?
                .2
        };

        bam_reader.fetch((chr, region_start, region_end))?;

        // boundary information for each haplotype
        let two_cp_haplotypes = final_haplotypes
            .values()
            .filter(|x| self.call.two_copy_haplotypes.contains(*x))
            .map(|x| x.to_string())
            .collect::<Vec<_>>();
        log::debug!("Two-copy haplotypes for VCF sample columns: {two_cp_haplotypes:?}");

        let nhap = if !is_gene2
            && self.call.heterozygous_sites.is_empty()
            && final_haplotypes.is_empty()
        {
            2
        } else {
            final_haplotypes.len() + two_cp_haplotypes.len()
        };
        let hap_info =
            self.build_hap_info(final_haplotypes, &two_cp_haplotypes, is_gene2, match_range)?;
        let special_variants = self.get_special_variants(is_gene2)?;

        let use_supplementary = self.phaser.use_supplementary();
        // variant info per hp per pos
        let mut variants_info: BTreeMap<i64, Vec<Option<VariantInfoByHP>>> = BTreeMap::new();
        let mut sv_variants: BTreeMap<i64, Vec<Option<String>>> = BTreeMap::new();
        // hp -> pos -> vec of bases
        let mut raw_piles: BTreeMap<String, BTreeMap<i64, BTreeMap<String, Vec<u8>>>> =
            BTreeMap::new();
        // store read seq
        let mut aln2seq = BTreeMap::new();
        for x in bam_reader.pileup() {
            let x = x?;
            let pos = i64::from(x.pos());
            if pos >= region_start && pos < region_end {
                let pileup_per_pos = query_seq_pileup(
                    &x,
                    &ref_seq,
                    &mut aln2seq,
                    25,
                    offset as usize,
                    use_supplementary,
                )?;
                for (hp, hp_bases) in pileup_per_pos.iter() {
                    for (read_name, read_base) in hp_bases.iter() {
                        raw_piles
                            .entry(hp.to_string())
                            .or_default()
                            .entry(pos)
                            .or_default()
                            .entry(read_name.to_string())
                            .or_insert(read_base.to_vec());
                    }
                }
            }
        }
        // unique reads
        let mut uniq_reads = Vec::new();
        for read_set in self.call.unique_supporting_reads.values() {
            for read_name in read_set {
                uniq_reads.push(read_name.to_string());
            }
        }
        log::debug!("Unique supporting reads collected for VCF consensus: {uniq_reads:?}");

        // follow order in final_haplotypes
        let empty_vec = vec![None; nhap];
        for pos in region_start..region_end {
            variants_info.entry(pos).or_insert(empty_vec.clone());
        }
        let mut hap_index = 0;
        for hap_name in final_haplotypes.values() {
            hap_index += 1;
            let hap_bound = hap_info.get(hap_index - 1).ok_or_else(|| {
                Exception::new(format!(
                    "Missing hap boundary for hap_index={} (hap='{}')",
                    hap_index, hap_name
                ))
            })?;

            let hap_start = hap_bound.start;
            let hap_end = hap_bound.end;
            let this_special_variant = special_variants.get(hap_name).cloned();
            let sv_bounds = this_special_variant
                .as_deref()
                .map(parse_symbolic_sv)
                .transpose()?;
            if let Some((start_1based, _sv_type, _end_1based)) = sv_bounds.as_ref() {
                let start_0based = *start_1based - 1;
                let entry = sv_variants.entry(start_0based).or_insert(vec![None; nhap]);
                entry[hap_index - 1] = this_special_variant.clone();
                if two_cp_haplotypes.contains(hap_name) {
                    entry[hap_index] = this_special_variant.clone();
                }
            }
            for pos in (hap_start + 1)..hap_end {
                if let Some((start_1based, _sv_type, end_1based)) = sv_bounds.as_ref() {
                    let true_pos_1based = pos + 1;
                    let skip_for_sv = if self.phaser.gene_name() == "ikbkg" {
                        true_pos_1based >= *start_1based && true_pos_1based <= *end_1based
                    } else {
                        true_pos_1based == *start_1based
                    };
                    if skip_for_sv {
                        continue;
                    }
                }
                if raw_piles.contains_key(hap_name) {
                    let hp_bases = raw_piles.get(hap_name).ok_or_else(|| {
                        Exception::new(format!(
                            "Haplotype '{}' missing from raw_piles while processing position {}",
                            hap_name, pos
                        ))
                    })?;
                    if let Some(pos_bases_by_read) = hp_bases.get(&pos) {
                        let mut pos_bases = Vec::new();
                        let mut pos_bases_uniq_reads = Vec::new();
                        for (read_name, read_base) in pos_bases_by_read.iter() {
                            pos_bases.push(read_base.clone());
                            if uniq_reads.contains(read_name) {
                                pos_bases_uniq_reads.push(read_base.clone());
                            }
                        }
                        let mut base_consensus: VariantInfoByHP =
                            get_consensus_var(pos_bases.clone(), pos - offset, &ref_seq)?;
                        // use uniq reads
                        if base_consensus.depth >= MIN_DEPTH {
                            if pos < hap_bound.start_strict
                                || pos > hap_bound.end_strict
                                || (base_consensus.nread as f32)
                                    < (base_consensus.depth as f32) * 0.7
                            {
                                let base_consensus_uniq_reads: VariantInfoByHP = get_consensus_var(
                                    pos_bases_uniq_reads.clone(),
                                    pos - offset,
                                    &ref_seq,
                                )?;
                                if base_consensus_uniq_reads.depth >= MIN_DEPTH
                                    && (base_consensus_uniq_reads.nread as f32)
                                        >= (base_consensus_uniq_reads.depth as f32) * 0.7
                                {
                                    base_consensus = base_consensus_uniq_reads;
                                }
                            }
                        }
                        /*
                        log::debug!(
                            "hap_name {hap_name} pos {pos} bases {:?} consensus {:?}",
                            pos_bases.clone(),
                            base_consensus.clone(),
                        );
                        */
                        let a: &mut Vec<std::option::Option<VariantInfoByHP>> =
                            variants_info.entry(pos).or_insert(empty_vec.clone());
                        let b = &mut a[hap_index - 1];
                        *b = Some(base_consensus.clone());
                        if two_cp_haplotypes.contains(hap_name) {
                            let b = &mut a[hap_index];
                            *b = Some(base_consensus.clone());
                        }
                    } else if pos - offset >= 0 && pos - offset < ref_seq_len as i64 {
                        let ref_base = vec![ref_base_at(&ref_seq, pos - offset)?];
                        let ref_base_string = std::str::from_utf8(&ref_base)?.to_string();
                        let none_var = VariantInfoByHP {
                            base: None,
                            ref_base: ref_base_string.clone(),
                            depth: 0,
                            nread: 0,
                            original_base: None,
                            bases_count: HashMap::new(),
                        };
                        let a: &mut Vec<std::option::Option<VariantInfoByHP>> =
                            variants_info.entry(pos).or_insert(empty_vec.clone());
                        let b = &mut a[hap_index - 1];
                        *b = Some(none_var.clone());
                        if two_cp_haplotypes.contains(hap_name) {
                            let b = &mut a[hap_index];
                            *b = Some(none_var.clone());
                        }
                    }
                } else if pos - offset >= 0 && pos - offset < ref_seq_len as i64 {
                    log::warn!(
                        "Haplotype missing from raw pileups while writing VCF: hap_name={hap_name}"
                    );
                    let ref_base = vec![ref_base_at(&ref_seq, pos - offset)?];
                    let ref_base_string = std::str::from_utf8(&ref_base)?.to_string();
                    let none_var = VariantInfoByHP {
                        base: None,
                        ref_base: ref_base_string.clone(),
                        depth: 0,
                        nread: 0,
                        original_base: None,
                        bases_count: HashMap::new(),
                    };
                    let a: &mut Vec<std::option::Option<VariantInfoByHP>> =
                        variants_info.entry(pos).or_insert(empty_vec.clone());
                    let b = &mut a[hap_index - 1];
                    *b = Some(none_var.clone());
                    if two_cp_haplotypes.contains(hap_name) {
                        let b = &mut a[hap_index];
                        *b = Some(none_var.clone());
                    }
                }
            }
            if two_cp_haplotypes.contains(hap_name) {
                hap_index += 1;
            }
        }

        // homozygous case
        if !is_gene2 && self.call.heterozygous_sites.is_empty() && final_haplotypes.is_empty() {
            let hap_name = &String::from("Unassigned");
            if raw_piles.contains_key(hap_name) {
                let hp_bases = raw_piles
                    .get(hap_name)
                    .ok_or_else(|| {
                        Exception::new(
                            "Homozygous mode expected 'Unassigned' entry in raw_piles but none was found",
                        )
                    })?;
                for pos in (region_start + 1)..region_end {
                    if let Some(pos_bases_by_read) = hp_bases.get(&pos) {
                        let mut pos_bases = Vec::new();
                        for (_read_name, read_base) in pos_bases_by_read.iter() {
                            pos_bases.push(read_base.clone());
                        }
                        let base_consensus: VariantInfoByHP =
                            get_consensus_var(pos_bases.clone(), pos - offset, &ref_seq)?;
                        /*
                        log::debug!(
                            "hap_name {hap_name} pos {pos} bases {:?} consensus {:?}",
                            pos_bases.clone(),
                            base_consensus.clone(),
                        );
                        */
                        let a: &mut Vec<std::option::Option<VariantInfoByHP>> =
                            variants_info.entry(pos).or_insert(empty_vec.clone());
                        let b = &mut a[0];
                        *b = Some(base_consensus.clone());
                        let b = &mut a[1];
                        *b = Some(base_consensus.clone());
                    } else if pos - offset >= 0 && pos - offset < ref_seq_len as i64 {
                        let ref_base = vec![ref_base_at(&ref_seq, pos - offset)?];
                        let ref_base_string = std::str::from_utf8(&ref_base)?.to_string();
                        let none_var = VariantInfoByHP {
                            base: None,
                            ref_base: ref_base_string.clone(),
                            depth: 0,
                            nread: 0,
                            original_base: None,
                            bases_count: HashMap::new(),
                        };
                        let a: &mut Vec<std::option::Option<VariantInfoByHP>> =
                            variants_info.entry(pos).or_insert(empty_vec.clone());
                        let b = &mut a[0];
                        *b = Some(none_var.clone());
                        let b = &mut a[1];
                        *b = Some(none_var.clone());
                    }
                }
            }
        }

        let hap_variant_info = HapVariantInfo {
            variants_info,
            sv_variants,
            hap_info,
        };
        Ok(hap_variant_info)
    }

    /// Get haplotypes for gene1 and gene2
    /// Split final haplotypes into gene1/gene2 groups for two-region genes.
    fn separate_two_genes(&self) -> (BTreeMap<String, String>, BTreeMap<String, String>) {
        let all_haplotypes = &self.call.final_haplotypes.clone();
        let mut gene1_haps = BTreeMap::new();
        let mut gene2_haps = BTreeMap::new();
        let gene_name = self.phaser.gene_name();
        if gene_name == "smn1" {
            for (hap, hap_name) in all_haplotypes {
                if hap_name.contains("smn1hap") {
                    gene1_haps.insert(hap.to_string(), hap_name.to_string());
                } else {
                    gene2_haps.insert(hap.to_string(), hap_name.to_string());
                }
            }
        } else if gene_name == "pms2" {
            for (hap, hap_name) in all_haplotypes {
                if hap_name.contains("cl") {
                    gene2_haps.insert(hap.to_string(), hap_name.to_string());
                } else {
                    gene1_haps.insert(hap.to_string(), hap_name.to_string());
                }
            }
        } else if gene_name == "ncf1" || gene_name == "ikbkg" {
            for (hap, hap_name) in all_haplotypes {
                if hap_name.contains("pseudo") {
                    gene2_haps.insert(hap.to_string(), hap_name.to_string());
                } else {
                    gene1_haps.insert(hap.to_string(), hap_name.to_string());
                }
            }
        } else if gene_name == "strc" {
            for (hap, hap_name) in all_haplotypes {
                if hap_name.contains("strcp1") {
                    gene2_haps.insert(hap.to_string(), hap_name.to_string());
                } else {
                    gene1_haps.insert(hap.to_string(), hap_name.to_string());
                }
            }
        }
        (gene1_haps, gene2_haps)
    }
}
