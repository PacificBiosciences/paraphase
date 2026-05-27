use crate::io::json::GeneCall;
use crate::toolkit::util::{DError, DResult};
use serde_json::Value;
use std::collections::BTreeMap;

type FusionMap = BTreeMap<String, BTreeMap<String, Value>>;

/// Parse `region_specific_info["fusions_called"]` into a typed nested map.
///
/// Returns `None` when the JSON shape does not match the expected
/// object-of-objects schema.
pub(crate) fn parse_fusions_called(value: &serde_json::Value) -> Option<FusionMap> {
    let object = value.as_object()?;
    let mut parsed = BTreeMap::new();
    for (k, v) in object {
        let inner_obj = v.as_object()?;
        let mut inner = BTreeMap::new();
        for (a, b) in inner_obj {
            inner.insert(a.to_string(), b.clone());
        }
        parsed.insert(k.to_string(), inner);
    }
    Some(parsed)
}

/// Parse and merge `fusions_called` maps from two gene calls.
///
/// Returns:
/// - `Ok(Some(_))` when both fusion maps exist and are valid.
/// - `Ok(None)` when either fusion map is missing.
/// - `Err(())` when both exist but one has an unexpected JSON shape.
pub(crate) fn merged_fusions_called(
    gene1_fusion: Option<&serde_json::Value>,
    gene2_fusion: Option<&serde_json::Value>,
) -> Result<Option<FusionMap>, ()> {
    let (Some(gene1_fusion), Some(gene2_fusion)) = (gene1_fusion, gene2_fusion) else {
        return Ok(None);
    };
    let gene1 = parse_fusions_called(gene1_fusion).ok_or(())?;
    let gene2 = parse_fusions_called(gene2_fusion).ok_or(())?;
    let mut merged = gene1;
    merged.extend(gene2);
    Ok(Some(merged))
}

/// Compute total copy number for CFHclust and whether to clear two-copy haplotypes.
#[must_use]
pub(crate) fn cfhclust_total_cn(
    gene1_cn: Option<i32>,
    gene2_cn: Option<i32>,
    gene1_num_haps: usize,
    gene2_num_haps: usize,
    has_fusions: bool,
) -> Option<(i32, bool)> {
    let (Some(gene1_cn), Some(gene2_cn)) = (gene1_cn, gene2_cn) else {
        return None;
    };
    let mut total_cn = gene1_cn.min(gene2_cn);
    let mut clear_two_copy_haps = false;
    let gene1_num_haps = gene1_num_haps as i32;
    let gene2_num_haps = gene2_num_haps as i32;
    if has_fusions && gene1_num_haps >= 2 && gene2_num_haps >= 2 {
        total_cn = total_cn.min(gene1_num_haps);
        total_cn = total_cn.min(gene2_num_haps);
        if total_cn < gene1_cn || total_cn < gene2_cn {
            clear_two_copy_haps = true;
        }
    }
    Some((total_cn, clear_two_copy_haps))
}

/// Shared logic for single-copy calls corrected by a partner gene with more haplotypes.
#[must_use]
pub(crate) fn adjusted_single_copy_call(
    phasing_results: &BTreeMap<String, GeneCall>,
    primary_gene: &str,
    partner_gene: &str,
    cn_field: &str,
) -> Option<GeneCall> {
    let primary_call = phasing_results.get(primary_gene)?;
    let partner_call = phasing_results.get(partner_gene)?;
    let primary_num_haps = primary_call.final_haplotypes.len();
    let partner_num_haps = partner_call.final_haplotypes.len();
    let primary_cn = primary_call
        .region_specific_info
        .get(cn_field)
        .and_then(serde_json::Value::as_i64)?;
    if primary_cn == 1 && partner_num_haps > primary_num_haps {
        let mut new_call = primary_call.clone();
        new_call
            .region_specific_info
            .insert(cn_field.to_string(), Option::<i32>::None.into());
        Some(new_call)
    } else {
        None
    }
}

/// Insert `new_call` into `phasing_results` for `gene` when present.
pub(crate) fn insert_adjusted_call(
    phasing_results: &mut BTreeMap<String, GeneCall>,
    gene: &str,
    new_call: Option<GeneCall>,
) {
    if let Some(call) = new_call {
        phasing_results.insert(gene.to_string(), call);
    }
}

/// Build the adjusted `smn1` call when `SERF1A` suggests the current CN=1 call is unreliable.
#[must_use]
pub(crate) fn adjusted_smn1_call(phasing_results: &BTreeMap<String, GeneCall>) -> Option<GeneCall> {
    adjusted_single_copy_call(phasing_results, "smn1", "SERF1A", "smn1_cn")
}

/// Build the adjusted `ncf1` call when `GTF2I` suggests the current CN=1 call is unreliable.
#[must_use]
pub(crate) fn adjusted_ncf1_call(phasing_results: &BTreeMap<String, GeneCall>) -> Option<GeneCall> {
    adjusted_single_copy_call(phasing_results, "ncf1", "GTF2I", "gene_cn")
}

/// Build the adjusted `TNXB` call when TNXB copy number exceeds `rccx`.
#[must_use]
pub(crate) fn adjusted_tnxb_call(phasing_results: &BTreeMap<String, GeneCall>) -> Option<GeneCall> {
    let gene1_call = phasing_results.get("TNXB")?;
    let gene2_call = phasing_results.get("rccx")?;
    let (Some(gene1_cn), Some(gene2_cn)) = (gene1_call.total_cn, gene2_call.total_cn) else {
        return None;
    };
    if gene1_cn > gene2_cn {
        let mut new_call = gene1_call.clone();
        new_call.total_cn = None;
        new_call.two_copy_haplotypes = Vec::new();
        Some(new_call)
    } else {
        None
    }
}

/// Build the derived `CFHclust` call from `CFH` and `CFHR3`.
///
/// # Errors
/// Returns an error if fusion structures are present but cannot be converted back to JSON.
pub(crate) fn build_cfhclust_call(
    phasing_results: &BTreeMap<String, GeneCall>,
) -> Result<Option<GeneCall>, DError> {
    let Some(gene1_call) = phasing_results.get("CFH") else {
        return Ok(None);
    };
    let Some(gene2_call) = phasing_results.get("CFHR3") else {
        return Ok(None);
    };
    let mut has_fusions = false;
    let gene1_haps = &gene1_call.final_haplotypes;
    let gene2_haps = &gene2_call.final_haplotypes;
    let gene1_two_cp = &gene1_call.two_copy_haplotypes;
    let gene2_two_cp = &gene2_call.two_copy_haplotypes;
    let gene1_phase_region = &gene1_call.phase_region;
    let gene2_phase_region = &gene2_call.phase_region;
    let gene1_genes_in_region = gene1_call.genes_in_region.as_deref().unwrap_or("");
    let gene2_genes_in_region = gene2_call.genes_in_region.as_deref().unwrap_or("");
    let gene1_fusion = gene1_call.region_specific_info.get("fusions_called");
    let gene2_fusion = gene2_call.region_specific_info.get("fusions_called");
    let gene1_cn = gene1_call.total_cn;
    let gene2_cn = gene2_call.total_cn;

    let mut new_call = GeneCall {
        gene_name: String::from("CFHclust"),
        sample_sex: gene1_call.sample_sex.clone(),
        genome_depth: gene1_call.genome_depth,
        phase_region: format!("{},{}", gene1_phase_region, gene2_phase_region,),
        genes_in_region: {
            let mut genes = gene1_genes_in_region
                .split(',')
                .chain(gene2_genes_in_region.split(','))
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .map(ToString::to_string)
                .collect::<Vec<_>>();
            genes.sort();
            genes.dedup();
            if genes.is_empty() {
                None
            } else {
                Some(genes.join(","))
            }
        },
        ..Default::default()
    };

    if !gene1_haps.is_empty() && !gene2_haps.is_empty() {
        new_call.final_haplotypes = gene1_haps.clone();
        new_call.final_haplotypes.extend(gene2_haps.clone());

        new_call.two_copy_haplotypes = gene1_two_cp.clone();
        new_call.two_copy_haplotypes.extend(gene2_two_cp.clone());
    }

    match merged_fusions_called(gene1_fusion, gene2_fusion) {
        Ok(Some(new_fusion)) => {
            if !new_fusion.is_empty() {
                has_fusions = true;
            }
            new_call.region_specific_info.insert(
                String::from("fusions_called"),
                serde_json::to_value(&new_fusion)?,
            );
        }
        Ok(None) => {}
        Err(()) => {
            log::warn!(
                "Skipping CFHclust fusion merge because `fusions_called` has an unsupported JSON shape."
            );
        }
    }

    if let Some((total_cn, clear_two_copy_haps)) = cfhclust_total_cn(
        gene1_cn,
        gene2_cn,
        gene1_haps.len(),
        gene2_haps.len(),
        has_fusions,
    ) {
        if clear_two_copy_haps {
            new_call.two_copy_haplotypes = Vec::new();
        }
        new_call.total_cn = Some(total_cn);
    }
    Ok(Some(new_call))
}

/// Apply cross-gene post-processing adjustments to per-gene calls.
///
/// This consolidates copy-number corrections and derived cluster calls that
/// depend on multiple gene results.
///
/// # Errors
/// Returns an error if `CFHclust` fusion maps cannot be serialized.
pub fn update_calls_after_per_gene_analysis(
    phasing_results: &mut BTreeMap<String, GeneCall>,
) -> DResult {
    insert_adjusted_call(phasing_results, "smn1", adjusted_smn1_call(phasing_results));
    insert_adjusted_call(phasing_results, "ncf1", adjusted_ncf1_call(phasing_results));
    insert_adjusted_call(phasing_results, "TNXB", adjusted_tnxb_call(phasing_results));

    if let Some(new_call) = build_cfhclust_call(phasing_results)? {
        phasing_results.insert(String::from("CFHclust"), new_call);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn make_call(name: &str) -> GeneCall {
        GeneCall {
            gene_name: name.to_string(),
            ..Default::default()
        }
    }

    #[test]
    fn merged_fusions_called_handles_missing_and_invalid_shapes() {
        let missing = merged_fusions_called(None, None).expect("missing fusion maps are valid");
        assert!(missing.is_none());

        let invalid = merged_fusions_called(Some(&json!(null)), Some(&json!({"x": {"a": "b"}})));
        assert!(invalid.is_err());
    }

    #[test]
    fn cfhclust_total_cn_clears_two_copy_haps_when_cn_reduced_by_fusions() {
        let out = cfhclust_total_cn(Some(3), Some(2), 2, 2, true)
            .expect("both CN values present should produce a total");
        assert_eq!(out, (2, true));
    }

    #[test]
    fn adjusted_single_copy_call_returns_none_when_partner_not_higher() {
        let mut calls = BTreeMap::<String, GeneCall>::new();
        let mut a = make_call("smn1");
        a.final_haplotypes
            .insert(String::from("h1"), String::from("smn1_hap1"));
        a.region_specific_info
            .insert(String::from("smn1_cn"), json!(1));
        calls.insert(String::from("smn1"), a);

        let mut b = make_call("SERF1A");
        b.final_haplotypes
            .insert(String::from("h1"), String::from("SERF1A_hap1"));
        calls.insert(String::from("SERF1A"), b);

        let out = adjusted_single_copy_call(&calls, "smn1", "SERF1A", "smn1_cn");
        assert!(out.is_none());
    }

    #[test]
    fn insert_adjusted_call_inserts_when_present() {
        let mut calls = BTreeMap::<String, GeneCall>::new();
        let call = make_call("smn1");
        insert_adjusted_call(&mut calls, "smn1", Some(call));
        assert!(calls.contains_key("smn1"));
    }

    #[test]
    fn update_calls_sets_smn1_cn_to_null_when_serf1_has_more_haps() -> DResult {
        let mut calls = BTreeMap::<String, GeneCall>::new();
        let mut serf1 = make_call("SERF1A");
        serf1
            .final_haplotypes
            .insert(String::from("a"), String::from("SERF1A_hap1"));
        serf1
            .final_haplotypes
            .insert(String::from("b"), String::from("SERF1A_hap2"));
        calls.insert(String::from("SERF1A"), serf1);

        let mut smn1 = make_call("smn1");
        smn1.final_haplotypes
            .insert(String::from("a"), String::from("smn1_hap1"));
        smn1.region_specific_info
            .insert(String::from("smn1_cn"), json!(1));
        calls.insert(String::from("smn1"), smn1);

        update_calls_after_per_gene_analysis(&mut calls)?;

        let updated = calls
            .get("smn1")
            .expect("smn1 should still be present after update");
        assert_eq!(
            updated.region_specific_info.get("smn1_cn"),
            Some(&json!(null))
        );
        Ok(())
    }

    #[test]
    fn update_calls_sets_ncf1_cn_to_null_when_gtf2i_has_more_haps() -> DResult {
        let mut calls = BTreeMap::<String, GeneCall>::new();
        let mut gtf2i = make_call("GTF2I");
        gtf2i
            .final_haplotypes
            .insert(String::from("a"), String::from("GTF2I_hap1"));
        gtf2i
            .final_haplotypes
            .insert(String::from("b"), String::from("GTF2I_hap2"));
        calls.insert(String::from("GTF2I"), gtf2i);

        let mut ncf1 = make_call("ncf1");
        ncf1.final_haplotypes
            .insert(String::from("a"), String::from("ncf1_hap1"));
        ncf1.region_specific_info
            .insert(String::from("gene_cn"), json!(1));
        calls.insert(String::from("ncf1"), ncf1);

        update_calls_after_per_gene_analysis(&mut calls)?;

        let updated = calls
            .get("ncf1")
            .expect("ncf1 should still be present after update");
        assert_eq!(
            updated.region_specific_info.get("gene_cn"),
            Some(&json!(null))
        );
        Ok(())
    }

    #[test]
    fn update_calls_clears_tnxb_cn_when_exceeding_rccx() -> DResult {
        let mut calls = BTreeMap::<String, GeneCall>::new();
        let mut tnx = make_call("TNXB");
        tnx.total_cn = Some(3);
        tnx.two_copy_haplotypes = vec![String::from("TNXB_hap1")];
        calls.insert(String::from("TNXB"), tnx);

        let mut rccx = make_call("rccx");
        rccx.total_cn = Some(2);
        calls.insert(String::from("rccx"), rccx);

        update_calls_after_per_gene_analysis(&mut calls)?;

        let updated = calls
            .get("TNXB")
            .expect("TNXB should still be present after update");
        assert_eq!(updated.total_cn, None);
        assert!(updated.two_copy_haplotypes.is_empty());
        Ok(())
    }

    #[test]
    fn update_calls_builds_cfhclust_from_cfh_and_cfhr3() -> DResult {
        let mut calls = BTreeMap::<String, GeneCall>::new();

        let mut cfh = make_call("CFH");
        cfh.phase_region = String::from("chr1:1-10");
        cfh.total_cn = Some(2);
        cfh.sample_sex = Some(String::from("male"));
        cfh.final_haplotypes
            .insert(String::from("h1"), String::from("CFH_hap1"));
        cfh.final_haplotypes
            .insert(String::from("h2"), String::from("CFH_hap2"));
        cfh.region_specific_info.insert(
            String::from("fusions_called"),
            json!({"f1": {"type": "deletion", "sequence": "1212", "breakpoint": [[10, 20], [30, 40]]}}),
        );
        calls.insert(String::from("CFH"), cfh);

        let mut cfhr3 = make_call("CFHR3");
        cfhr3.phase_region = String::from("chr1:11-20");
        cfhr3.total_cn = Some(1);
        cfhr3
            .final_haplotypes
            .insert(String::from("h3"), String::from("CFHR3_hap1"));
        cfhr3.region_specific_info.insert(
            String::from("fusions_called"),
            json!({"f2": {"type": null, "sequence": "2121", "breakpoint": [[50, 60], [70, 80]]}}),
        );
        calls.insert(String::from("CFHR3"), cfhr3);

        update_calls_after_per_gene_analysis(&mut calls)?;

        let merged = calls
            .get("CFHclust")
            .expect("CFHclust should be created when CFH and CFHR3 exist");
        assert_eq!(merged.gene_name, "CFHclust");
        assert_eq!(merged.phase_region, "chr1:1-10,chr1:11-20");
        assert_eq!(merged.total_cn, Some(1));
        assert_eq!(merged.final_haplotypes.len(), 3);
        assert_eq!(
            merged.region_specific_info.get("fusions_called"),
            Some(&json!({
                "f1": {"type": "deletion", "sequence": "1212", "breakpoint": [[10, 20], [30, 40]]},
                "f2": {"type": null, "sequence": "2121", "breakpoint": [[50, 60], [70, 80]]}
            }))
        );
        Ok(())
    }
}
