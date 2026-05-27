use crate::assembly::assembly_result::AssembledPaths;
use crate::config::{self, Region as RegionConfig};
use crate::toolkit::util::DError;

use vstr::VString;

use std::collections::BTreeMap;

use super::errors::{missing_required, should_be_found_debug};
use super::validate::expect_object;
use super::{Error, JsonMap, ParaphaseOutput, PhasingSite, ReadAlignmentId, ReadFingerprintMap};

/// Require that an object-field entry is a string.
///
/// # Errors
/// Returns an error with field/key context when the value is not a string.
fn expect_string_in_object_field<'a>(
    entry: &'a serde_json::Value,
    field: &str,
    key: &str,
) -> Result<&'a str, DError> {
    entry
        .as_str()
        .ok_or_else(|| DError::from(Error::new(format!("{field}['{key}'] should be a string"))))
}

/// Require that an object-field entry is a string, including the found value in the message.
///
/// # Errors
/// Returns an error with field/key context when the value is not a string.
fn expect_string_in_object_field_with_found<'a>(
    entry: &'a serde_json::Value,
    field: &str,
    key: &str,
) -> Result<&'a str, DError> {
    entry.as_str().ok_or_else(|| {
        DError::from(Error::new(format!(
            "{field}['{key}'] should be a string, found {entry:?}"
        )))
    })
}

/// Require that an array-field entry is a string.
///
/// # Errors
/// Returns an error with field/index context when the value is not a string.
fn expect_string_in_array_field<'a>(
    entry: &'a serde_json::Value,
    field: &str,
    idx: usize,
) -> Result<&'a str, DError> {
    entry.as_str().ok_or_else(|| {
        DError::from(Error::new(format!(
            "{field}[{idx}] should be a string, found {entry:?}"
        )))
    })
}

/// Parse `read_details` into a read->hap map.
///
/// # Errors
/// Returns an error if the field is missing or malformed.
fn parse_read_to_hap(value: &JsonMap) -> Result<ReadFingerprintMap, DError> {
    let read_details = value
        .get("read_details")
        .ok_or_else(|| missing_required("read_details"))?;
    let mut ret = ReadFingerprintMap::new();
    match read_details {
        serde_json::Value::Object(read_details) => {
            for (read_name, hap) in read_details {
                ret.insert(
                    ReadAlignmentId::from_name(read_name),
                    VString::from(expect_string_in_object_field(
                        hap,
                        "read_details",
                        read_name,
                    )?),
                );
            }
        }
        serde_json::Value::Null => {}
        other => {
            return Err(
                should_be_found_debug("read_details", "a JSON object or null", other).into(),
            );
        }
    }
    Ok(ret)
}

/// Parse `final_haplotypes`.
///
/// Returns `None` when field is absent/null.
///
/// # Errors
/// Returns an error if the field exists but has invalid type or values.
fn parse_final_haps(value: &JsonMap) -> Result<Option<BTreeMap<VString, String>>, DError> {
    Ok(match value.get("final_haplotypes") {
        None | Some(serde_json::Value::Null) => None,
        Some(serde_json::Value::Object(dict)) => {
            let mut ret = BTreeMap::<VString, String>::new();
            for (key, val) in dict {
                let val = expect_string_in_object_field_with_found(val, "final_haplotypes", key)?;
                ret.insert(VString::from(key), val.into());
            }
            Some(ret)
        }
        Some(other) => {
            return Err(
                should_be_found_debug("final_haplotypes", "an object or null", other).into(),
            );
        }
    })
}

/// Parse `assembled_haplotypes`.
///
/// Returns `None` when field is absent/null.
///
/// # Errors
/// Returns an error if the field exists but has invalid type or values.
fn parse_hap_asm(value: &JsonMap) -> Result<Option<AssembledPaths>, DError> {
    Ok(match value.get("assembled_haplotypes") {
        None => None,
        Some(serde_json::Value::Null) => None,
        Some(serde_json::Value::Array(arr)) => {
            let mut seqs = Vec::with_capacity(arr.len());
            for (idx, entry) in arr.iter().enumerate() {
                let hap = expect_string_in_array_field(entry, "assembled_haplotypes", idx)?;
                seqs.push(hap);
            }
            Some(AssembledPaths::from_seqs(seqs))
        }
        Some(other) => {
            return Err(
                should_be_found_debug("assembled_haplotypes", "an array or null", other).into(),
            );
        }
    })
}

/// Parse `sites_for_phasing`.
///
/// Returns an empty vector when field is absent/null.
///
/// # Errors
/// Returns an error if the field exists but has invalid type or values.
fn parse_sites_for_phasing(value: &JsonMap) -> Result<Vec<PhasingSite>, DError> {
    match value.get("sites_for_phasing") {
        None | Some(serde_json::Value::Null) => Ok(Vec::new()),
        Some(serde_json::Value::Array(items)) => items
            .iter()
            .enumerate()
            .map(|(idx, item)| {
                let raw = expect_string_in_array_field(item, "sites_for_phasing", idx)?;
                PhasingSite::try_new(raw)
            })
            .collect::<Result<Vec<_>, DError>>(),
        Some(other) => {
            Err(should_be_found_debug("sites_for_phasing", "an array or null", other).into())
        }
    }
}

/// Decode a `ParaphaseOutput` from JSON.
///
/// # Errors
/// Returns an error if `value` is not the expected object schema.
pub(super) fn parse_paraphase_output_from_json(
    value: &serde_json::Value,
    gene: String,
    config: Option<&RegionConfig>,
) -> Result<ParaphaseOutput, DError> {
    let dict = expect_object(value, "ParaphaseOutput::from_json")?;
    let read_to_hap = parse_read_to_hap(dict)?;
    let assembled_haps = parse_hap_asm(dict)?;
    let final_haps = parse_final_haps(dict)?;
    let sites_for_phasing = parse_sites_for_phasing(dict)?;
    let pivot_site = config
        .unwrap_or(&config::CONFIG)
        .get(&gene[..])
        .and_then(|x| x.get("pivot_site").and_then(serde_yaml::Value::as_i64));
    Ok(ParaphaseOutput {
        pivot_site,
        gene,
        read_to_hap,
        assembled_haps,
        final_haps,
        sites_for_phasing,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn decode_minimal_ok() -> Result<(), DError> {
        let value = json!({
            "read_details": {
                "r1": "121x"
            }
        });
        let out = parse_paraphase_output_from_json(&value, String::from("UNKNOWN_GENE"), None)?;
        assert_eq!(out.gene, "UNKNOWN_GENE");
        assert_eq!(out.read_to_hap.len(), 1);
        assert_eq!(
            out.read_to_hap
                .get(&ReadAlignmentId::from_name("r1"))
                .map(|x| x.as_slice()),
            Some(b"121x".as_slice())
        );
        assert!(out.assembled_haps.is_none());
        assert!(out.final_haps.is_none());
        assert!(out.sites_for_phasing.is_empty());
        Ok(())
    }

    #[test]
    fn decode_optional_fields_ok() -> Result<(), DError> {
        let value = json!({
            "read_details": {
                "r1": "12"
            },
            "assembled_haplotypes": ["111", "222"],
            "final_haplotypes": {
                "111": "GENE_hap1"
            },
            "sites_for_phasing": ["100_A_T", "200_del5"]
        });
        let out = parse_paraphase_output_from_json(&value, String::from("GENE"), None)?;
        assert_eq!(
            out.assembled_haps,
            Some(AssembledPaths::from_seqs(["111", "222"]))
        );
        assert_eq!(out.final_haps.as_ref().map(BTreeMap::len), Some(1));
        assert_eq!(out.sites_for_phasing.len(), 2);
        assert_eq!(out.sites_for_phasing[0].to_string(), "100_A_T");
        assert_eq!(out.sites_for_phasing[1].deletion, "del5");
        Ok(())
    }

    #[test]
    fn decode_requires_read_details() {
        let value = json!({});
        let err = parse_paraphase_output_from_json(&value, String::from("GENE"), None)
            .expect_err("missing read_details should fail")
            .to_string();
        assert!(
            err.contains("Missing required field 'read_details'"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn decode_rejects_bad_assembled_haplotypes() {
        let value = json!({
            "read_details": {
                "r1": "12"
            },
            "assembled_haplotypes": [1]
        });
        let err = parse_paraphase_output_from_json(&value, String::from("GENE"), None)
            .expect_err("bad assembled_haplotypes should fail")
            .to_string();
        assert!(
            err.contains("assembled_haplotypes[0] should be a string"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn decode_rejects_bad_sites_for_phasing_item() {
        let value = json!({
            "read_details": {
                "r1": "12"
            },
            "sites_for_phasing": [1]
        });
        let err = parse_paraphase_output_from_json(&value, String::from("GENE"), None)
            .expect_err("bad sites_for_phasing should fail")
            .to_string();
        assert!(
            err.contains("sites_for_phasing[0] should be a string"),
            "unexpected error: {err}"
        );
    }
}
