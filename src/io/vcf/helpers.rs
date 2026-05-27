//! Pure helpers and constants for VCF generation.
//!
//! Includes:
//! - VCF header line constants
//! - consensus/GT helper logic
//! - symbolic SV parsing and per-sample field helpers
//! - pileup-base extraction utilities

use crate::phaser::{Exception, Phaser};
use crate::toolkit::util::DError;
use anyhow::anyhow;
use regex::Regex;
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::{bam, bam::record::Aux, htslib};
use std::collections::{BTreeMap, HashMap, HashSet};

use super::types::{HapBoundForVcf, VariantInfoByHP};

type SampleFormatFields = (
    Vec<GenotypeAllele>,
    Vec<GenotypeAllele>,
    Vec<Vec<u8>>,
    Vec<Vec<u8>>,
);

/// Header lines defining the INFO and FORMAT fields for the VCF file.
pub(crate) const VCF_LINES: [&str; 6] = [
    r#"##FILTER=<ID=PASS,Description=\"All filters passed\">"#,
    r#"##FILTER=<ID=LowQual,Description=\"Nonpassing variant\">"#,
    r#"##INFO=<ID=HPBOUND,Number=.,Type=String,Description=\"Boundary coordinates of the phased haplotype\">"#,
    r#"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype per haplotype\">"#,
    r#"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth per haplotype\">"#,
    r#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read depth for REF and ALT per haplotype\">"#,
];
pub(crate) const VCF_LINE_ALLELE: [&str; 1] =
    [r#"##INFO=<ID=ALLELE,Number=.,Type=String,Description=\"Haplotypes phased into alleles\">"#];
pub(crate) const VCF_LINES_SV: [&str; 6] = [
    r#"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV\">"#,
    r#"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">"#,
    r#"##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">"#,
    r#"##ALT=<ID=DEL,Description=\"Deletion\">"#,
    r#"##ALT=<ID=DUP,Description=\"Duplication\">"#,
    r#"##ALT=<ID=INV,Description=\"Inversion\">"#,
];
/// minimum depth for variant calling
pub(crate) const MIN_DEPTH: usize = 4;

/// Safely get reference base at index. Returns a descriptive error instead of panicking on OOB.
pub(crate) fn ref_base_at(ref_seq: &[u8], idx: i64) -> Result<u8, DError> {
    let idx = idx as usize;
    ref_seq.get(idx).copied().ok_or_else(|| {
        anyhow!(
            "Reference sequence index out of bounds: index={idx} len={} (position may extend beyond reference region)",
            ref_seq.len()
        )
        .into()
    })
}

/// Parse a symbolic structural-variant name in Paraphase's internal format.
///
/// Expected format is `"{start}_{SVTYPE}_{end}"`, using 1-based coordinates.
pub(crate) fn parse_symbolic_sv(var_name: &str) -> Result<(i64, String, i64), DError> {
    let parts = var_name.split('_').collect::<Vec<_>>();
    if parts.len() != 3 {
        return Err(Exception::new(format!(
            "Malformed symbolic SV name '{var_name}'; expected 'start_type_end'"
        ))
        .into());
    }
    let start = parts[0].parse::<i64>()?;
    let sv_type = parts[1].to_string();
    let end = parts[2].parse::<i64>()?;
    Ok((start, sv_type, end))
}

/// Convert indel representations to compatible with pileups
/// A, ACT -> A, A+2CT
/// ACT, A -> A, A-2CT
pub(crate) fn convert_alt_record(ref_base: String, alt_base: String) -> String {
    let alt_len = alt_base.len();
    if alt_len > 1 {
        let ins_len = alt_len - 1;
        let ins_seq = &alt_base[1..alt_len];
        return format!("{ref_base}+{ins_len}{ins_seq}");
    }
    let ref_len = ref_base.len();
    if ref_len > 1 {
        let del_len = ref_base.len() - 1;
        let del_seq = &ref_base[1..ref_len];
        let ref_base_first = &ref_base[0..1];
        return format!("{ref_base_first}-{del_len}{del_seq}");
    }
    alt_base
}

pub(crate) fn format_hap_bound(
    hap_start_1based: i64,
    hap_end_1based: i64,
    is_truncated: &[String],
) -> String {
    let hap_start = hap_start_1based.to_string();
    let hap_end = hap_end_1based.to_string();
    if is_truncated.is_empty() {
        return format!("{hap_start}-{hap_end}");
    }
    if is_truncated.contains(&String::from("5p")) && !is_truncated.contains(&String::from("3p")) {
        return format!("{hap_start}truncated-{hap_end}");
    }
    if !is_truncated.contains(&String::from("5p")) && is_truncated.contains(&String::from("3p")) {
        return format!("{hap_start}-{hap_end}truncated");
    }
    format!("{hap_start}truncated-{hap_end}truncated")
}

/// Given a set of variants at a position, return whether this pos is reference only
pub(crate) fn test_ref_only(variant_observed: &HashSet<(String, String)>) -> Result<bool, DError> {
    let mut ref_only = false;
    if variant_observed.len() == 1 {
        let (ref_base, alt_base) = variant_observed
            .iter()
            .next()
            .ok_or_else(|| Exception::new("variant_observed was unexpectedly empty"))?;
        if ref_base == alt_base {
            ref_only = true;
        }
    } else if variant_observed.len() == 2 {
        let mut has_ref = false;
        let mut has_del = false;
        for (ref_base, alt_base) in variant_observed {
            if ref_base == alt_base {
                has_ref = true
            } else if *alt_base == "*" {
                has_del = true;
            }
        }
        if has_ref && has_del {
            ref_only = true;
        }
    }
    Ok(ref_only)
}

/// Return genotype given base counts at a site
pub(crate) fn get_gt(
    this_hap_call: &VariantInfoByHP,
    ref_base: &String,
    alt_base: &String,
) -> GenotypeAllele {
    if let Some(ref variant_call) = this_hap_call.base {
        let depth = this_hap_call.depth;
        if depth >= MIN_DEPTH && this_hap_call.nread as f32 >= (depth as f32) * 0.7 {
            if *variant_call == this_hap_call.ref_base {
                return GenotypeAllele::Phased(0);
            }
            if variant_call == alt_base && this_hap_call.ref_base == *ref_base {
                return GenotypeAllele::Phased(1);
            }
        }
    }
    GenotypeAllele::UnphasedMissing
}

/// Return whether a pileup-derived haplotype call is a confident reference call.
///
/// This mirrors the Python merge behavior for symbolic SV records, where a
/// non-SV haplotype should be emitted as `0` only if its regular pileup call
/// at the same anchor position is truly reference.
pub(crate) fn is_confident_ref_call(this_hap_call: &VariantInfoByHP) -> bool {
    if let Some(ref variant_call) = this_hap_call.base {
        let depth = this_hap_call.depth;
        depth >= MIN_DEPTH
            && this_hap_call.nread as f32 >= (depth as f32) * 0.7
            && *variant_call == this_hap_call.ref_base
    } else {
        false
    }
}

/// Return `DP` / `AD` strings for a symbolic SV record from the pileup-derived
/// call at the same anchor position.
///
/// For a confident reference call, the symbolic SV ALT is absent by definition,
/// so the ALT depth is reported as `0`.
pub(crate) fn symbolic_sv_dp_ad(this_hap_call: &VariantInfoByHP) -> (Vec<u8>, Vec<u8>) {
    let dp = this_hap_call.depth.to_string().into_bytes();
    let ref_count = *this_hap_call
        .bases_count
        .get(this_hap_call.ref_base.as_bytes())
        .unwrap_or(&0);
    let ad = format!("{ref_count},0").into_bytes();
    (dp, ad)
}

/// Get consensus among a set of bases
/// # Arguments
/// * `bases` - all bases at this position, each read at this position is a vec<u8>
/// * `pos` - position on repeat unit
/// * `ref_seq` - reference sequence
/// # Returns
/// * `VariantInfoByHP` - variant at this position from this set of reads
pub(crate) fn get_consensus_var(
    bases: Vec<Vec<u8>>,
    pos: i64,
    ref_seq: &[u8],
) -> Result<VariantInfoByHP, DError> {
    let depth = bases.len();
    let ref_base = vec![ref_base_at(ref_seq, pos)?];
    let ref_base_string = std::str::from_utf8(&ref_base)?.to_string();
    let none_var = VariantInfoByHP {
        base: None,
        ref_base: ref_base_string.clone(),
        depth,
        nread: 0,
        original_base: None,
        bases_count: HashMap::new(),
    };
    let mut bases_count: HashMap<Vec<u8>, usize> = HashMap::new();
    for base in bases {
        *bases_count.entry(base).or_default() += 1;
    }
    let mut bases_count2 = bases_count
        .iter()
        .map(|(x, y)| (x.clone(), *y))
        .collect::<Vec<_>>();
    bases_count2.sort_by(|a, b| {
        b.1.cmp(&a.1)
            .then(a.0.len().cmp(&b.0.len()))
            .then((b.0 == ref_base).cmp(&(a.0 == ref_base)))
            .then(a.0.cmp(&b.0))
    });
    if bases_count2.is_empty() {
        return Ok(none_var.clone());
    }

    let base_consensus = bases_count2[0].0.clone();
    let base_consensus_string = std::str::from_utf8(&base_consensus)?.to_string();
    if base_consensus.len() > 1 && base_consensus.contains(&b'-') {
        let re = Regex::new(r"\-\d+")?;
        let del_seq = re
            .split(&base_consensus_string)
            .collect::<Vec<&str>>()
            .last()
            .ok_or_else(|| {
                Exception::new(format!(
                    "Failed to parse deletion consensus '{}' with deletion regex split",
                    base_consensus_string
                ))
            })?
            .to_string();
        let mut new_ref_base = ref_base.clone();
        new_ref_base.extend_from_slice(del_seq.as_bytes());
        return Ok(VariantInfoByHP {
            base: Some(ref_base_string.clone()),
            ref_base: std::str::from_utf8(&new_ref_base)?.to_string(),
            depth,
            nread: bases_count2[0].1,
            original_base: Some(base_consensus_string),
            bases_count,
        });
    }
    if base_consensus.len() > 1 && base_consensus.contains(&b'+') {
        let re = Regex::new(r"\+\d+")?;
        let ins_seq = re
            .split(&base_consensus_string)
            .collect::<Vec<&str>>()
            .last()
            .ok_or_else(|| {
                Exception::new(format!(
                    "Failed to parse insertion consensus '{}' with insertion regex split",
                    base_consensus_string
                ))
            })?
            .to_string();
        let mut variant_base = ref_base.clone();
        variant_base.extend_from_slice(ins_seq.as_bytes());
        return Ok(VariantInfoByHP {
            base: Some(std::str::from_utf8(&variant_base)?.to_string()),
            ref_base: ref_base_string.clone(),
            depth,
            nread: bases_count2[0].1,
            original_base: Some(base_consensus_string),
            bases_count,
        });
    }
    if base_consensus_string != "x" && base_consensus_string != "-" {
        return Ok(VariantInfoByHP {
            base: Some(base_consensus_string.clone()),
            ref_base: ref_base_string.clone(),
            depth,
            nread: bases_count2[0].1,
            original_base: Some(base_consensus_string.clone()),
            bases_count,
        });
    }
    Ok(none_var.clone())
}

/// Return the raw position of the current alignment on the read
/// # Arguments
/// * `x` - reference to an alignment
fn raw_qpos<'a>(x: &'a bam::pileup::Alignment<'a>) -> usize {
    static_assertions::const_assert_eq!(
        std::mem::size_of::<&bam::pileup::Alignment<'_>>(),
        std::mem::size_of::<&htslib::bam_pileup1_t>()
    );
    let ptr = (x as *const bam::pileup::Alignment<'_>).cast::<&htslib::bam_pileup1_t>();
    unsafe { *ptr }.qpos as usize
}

/// Given a pileup at one position, return bases grouped by HPs
pub(crate) fn query_seq_pileup(
    x: &bam::pileup::Pileup,
    ref_seq: &[u8],
    aln2seq: &mut BTreeMap<String, Vec<u8>>,
    min_base_quality: u8,
    offset: usize,
    use_supplementary: bool,
) -> Result<BTreeMap<String, BTreeMap<String, Vec<u8>>>, DError> {
    let mut pileups_raw: BTreeMap<String, BTreeMap<String, Vec<u8>>> = BTreeMap::new();
    let pos = x.pos();
    for aln in x.alignments() {
        let query_pos_raw = raw_qpos(&aln);
        let record = aln.record();
        let updated_qname = Phaser::get_read_name_free(&record, use_supplementary);
        match record.aux(b"HP") {
            Ok(value) => {
                if let Aux::String(hp_value) = value {
                    let bq = record.qual().get(query_pos_raw).copied().unwrap_or(0);
                    let entry = aln2seq
                        .entry(updated_qname.to_string())
                        .or_insert_with(|| record.seq().as_bytes());
                    let seq: &[u8] = &entry[..];
                    let mut query_seq = Vec::new();

                    let ambiguous_bases = vec![b'-', b'N', b'<', b'>', b'*'];
                    if bq >= min_base_quality && !aln.is_refskip() {
                        let base = if !aln.is_del() {
                            seq.get(query_pos_raw).copied().unwrap_or(b'N')
                        } else {
                            b'*'
                        };
                        query_seq.push(base);
                        let pos = pos as usize;
                        if !ambiguous_bases.contains(&base) {
                            match aln.indel() {
                                bam::pileup::Indel::Ins(x) => {
                                    debug_assert!(x > 0);
                                    query_seq.push(b'+');
                                    query_seq.extend_from_slice(x.to_string().as_bytes());
                                    for j in 1..=(x as usize) {
                                        query_seq.push(seq[j + query_pos_raw]);
                                    }
                                }
                                bam::pileup::Indel::Del(x) => {
                                    debug_assert!(x > 0);
                                    query_seq.push(b'-');
                                    query_seq.extend_from_slice(x.to_string().as_bytes());
                                    for j in 1..=(x as usize) {
                                        query_seq.push(ref_base_at(
                                            ref_seq,
                                            (j as i64) + (pos as i64) - (offset as i64),
                                        )?);
                                    }
                                }
                                bam::pileup::Indel::None => {}
                            }
                        }
                        let query_seq_vec = query_seq.clone().to_vec();

                        pileups_raw
                            .entry(hp_value.to_string())
                            .or_default()
                            .entry(updated_qname.to_string())
                            .or_insert_with(|| query_seq_vec.clone());
                    }
                }
            }
            Err(_e) => {
                log::debug!("Missing HP tag for read {updated_qname}")
            }
        }
    }
    Ok(pileups_raw)
}

/// Return whether a no-call at `pos0` should be considered valid for emission
/// because it lies beyond a haplotype's confident boundary (including truncation rules).
pub(crate) fn no_call_beyond_haplotype_boundary(hap: &HapBoundForVcf, pos0: i64) -> bool {
    let hap_info_truncated = &hap.is_truncated;
    if hap_info_truncated.is_empty() {
        return true;
    }
    if hap_info_truncated.contains(&String::from("5p"))
        && !hap_info_truncated.contains(&String::from("3p"))
    {
        return pos0 > hap.start;
    }
    if !hap_info_truncated.contains(&String::from("5p"))
        && hap_info_truncated.contains(&String::from("3p"))
    {
        return pos0 < hap.end;
    }
    pos0 > hap.start && pos0 < hap.end
}

/// Build per-sample `GT`/`DP`/`AD` fields for one symbolic SV record.
///
/// Returns `(valid_gts, gts, dps, ads)` where `valid_gts` is used only for
/// write/no-write gating in no-call mode.
pub(crate) fn sv_sample_fields(
    variant_name: &str,
    start_1based: i64,
    pos_calls: &[Option<String>],
    hap_info: &[HapBoundForVcf],
    variants_info: &BTreeMap<i64, Vec<Option<VariantInfoByHP>>>,
) -> SampleFormatFields {
    let mut valid_gts = Vec::new();
    let mut gts = Vec::new();
    let mut dps = Vec::new();
    let mut ads = Vec::new();
    for (hap_index, hap) in hap_info.iter().enumerate() {
        let this_hap = &pos_calls[hap_index];
        let pileup_call = variants_info
            .get(&(start_1based - 1))
            .and_then(|calls| calls.get(hap_index))
            .and_then(|call| call.as_ref());
        let pileup_ref_call = pileup_call.is_some_and(is_confident_ref_call);
        if let Some(this_variant) = this_hap {
            if *this_variant == variant_name {
                gts.push(GenotypeAllele::Phased(1));
                valid_gts.push(GenotypeAllele::Phased(1));
                dps.push(String::from(".").into_bytes());
                ads.push(String::from(".").into_bytes());
            } else if pileup_ref_call {
                let (dp, ad) =
                    symbolic_sv_dp_ad(pileup_call.expect("pileup_ref_call implies pileup_call"));
                gts.push(GenotypeAllele::Phased(0));
                valid_gts.push(GenotypeAllele::Phased(0));
                dps.push(dp);
                ads.push(ad);
            } else {
                gts.push(GenotypeAllele::UnphasedMissing);
                dps.push(String::from(".").into_bytes());
                ads.push(String::from(".").into_bytes());
            }
        } else {
            let pos0 = start_1based - 1;
            if pileup_ref_call {
                let (dp, ad) =
                    symbolic_sv_dp_ad(pileup_call.expect("pileup_ref_call implies pileup_call"));
                gts.push(GenotypeAllele::Phased(0));
                valid_gts.push(GenotypeAllele::Phased(0));
                dps.push(dp);
                ads.push(ad);
            } else {
                gts.push(GenotypeAllele::UnphasedMissing);
                dps.push(String::from(".").into_bytes());
                ads.push(String::from(".").into_bytes());
            }
            if !pileup_ref_call && no_call_beyond_haplotype_boundary(hap, pos0) {
                valid_gts.push(GenotypeAllele::UnphasedMissing);
            }
        }
    }
    (valid_gts, gts, dps, ads)
}

/// Build per-sample `GT`/`DP`/`AD` fields for one small-variant (non-symbolic) record.
///
/// Returns `(valid_gts, gts, dps, ads)` where `valid_gts` is used only for
/// write/no-write gating in no-call mode.
pub(crate) fn small_variant_sample_fields(
    pos: i64,
    ref_base: &str,
    alt_base: &str,
    pos_calls: &[Option<VariantInfoByHP>],
    hap_info: &[HapBoundForVcf],
) -> SampleFormatFields {
    let mut valid_gts = Vec::new();
    let mut gts = Vec::new();
    let mut dps = Vec::new();
    let mut ads = Vec::new();
    for (hap_index, hap) in hap_info.iter().enumerate() {
        let this_hap = &pos_calls[hap_index];
        if let Some(this_hap_call) = this_hap {
            let this_dp = this_hap_call.depth.to_string();
            let this_counter = this_hap_call.clone().bases_count;
            let ori_ref = &ref_base[0..1];
            let converted_alt = convert_alt_record(ref_base.to_string(), alt_base.to_string());
            let ref_base_count = *this_counter.get(ori_ref.as_bytes()).unwrap_or(&0);
            let alt_base_count = *this_counter.get(converted_alt.as_bytes()).unwrap_or(&0);
            if ref_base != alt_base {
                ads.push(format!("{ref_base_count},{alt_base_count}").into_bytes());
            } else {
                let alt_count = this_hap_call.depth - ref_base_count;
                ads.push(format!("{ref_base_count},{alt_count}").into_bytes());
            }
            dps.push(this_dp.into_bytes());
            let this_gt = get_gt(this_hap_call, &ref_base.to_string(), &alt_base.to_string());
            gts.push(this_gt);
            valid_gts.push(this_gt);
        } else {
            gts.push(GenotypeAllele::UnphasedMissing);
            dps.push(String::from(".").into_bytes());
            ads.push(String::from(".").into_bytes());
            if no_call_beyond_haplotype_boundary(hap, pos) {
                valid_gts.push(GenotypeAllele::UnphasedMissing);
            }
        }
    }
    (valid_gts, gts, dps, ads)
}

/// Pad sample fields when writing two-region genes so output columns match
/// the full haplotype/sample list order.
pub(crate) fn pad_sample_formats_for_counter(
    counter: usize,
    haps_ids: &[String],
    haps_ids1: &[String],
    haps_ids2: &[String],
    gts: &mut Vec<GenotypeAllele>,
    dps: &mut Vec<Vec<u8>>,
    ads: &mut Vec<Vec<u8>>,
) {
    if counter == 0 && haps_ids.len() > haps_ids1.len() {
        for _i in 0..haps_ids2.len() {
            gts.push(GenotypeAllele::UnphasedMissing);
            dps.push(String::from(".").into_bytes());
            ads.push(String::from(".").into_bytes());
        }
    } else if counter > 0 {
        for _i in 0..haps_ids1.len() {
            gts.insert(0, GenotypeAllele::UnphasedMissing);
            dps.insert(0, String::from(".").into_bytes());
            ads.insert(0, String::from(".").into_bytes());
        }
    }
}
