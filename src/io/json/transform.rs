use crate::config::Region as RegionConfig;
use crate::toolkit::util::DError;

use std::collections::BTreeMap;

use super::decode::parse_paraphase_output_from_json;
use super::errors::not_dictionary_for;
use super::{Error, JsonMap, ParsedParaphaseOutputJSON};

/// Build a parsed output wrapper from a JSON value.
///
/// # Errors
/// Returns an error if the root is not a JSON object.
pub(super) fn parse_parsed_output_from_json(
    object: serde_json::Value,
    config: Option<&RegionConfig>,
) -> Result<ParsedParaphaseOutputJSON, DError> {
    let serde_json::Value::Object(object) = object else {
        return Err(not_dictionary_for("ParsedParaphaseOutputJSON").into());
    };
    log::debug!("Parsing JSON object into internal representation: {object:?}");
    parse_parsed_output_from_object(object, config)
}

/// Build a parsed output wrapper from a JSON object.
///
/// # Errors
/// Returns an error if any gene block fails schema decoding.
pub(super) fn parse_parsed_output_from_object(
    object: JsonMap,
    config: Option<&RegionConfig>,
) -> Result<ParsedParaphaseOutputJSON, DError> {
    let mut gene_data = BTreeMap::new();
    for (key, value) in object.iter().map(|(key, value)| {
        (
            String::from(key),
            parse_paraphase_output_from_json(value, String::from(key), config),
        )
    }) {
        let parsed = value
            .map_err(|e| DError::from(Error::new(format!("Failed parsing gene '{key}': {e}"))))?;
        gene_data.insert(key, parsed);
    }
    Ok(ParsedParaphaseOutputJSON { object, gene_data })
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn parse_parsed_output_from_json_rejects_non_object() {
        let err = parse_parsed_output_from_json(json!(["x"]), None)
            .expect_err("root array should fail")
            .to_string();
        assert!(
            err.contains("Not a dictionary for ParsedParaphaseOutputJSON"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn parse_parsed_output_from_object_wraps_gene_context() {
        let object = serde_json::Map::from_iter([(
            String::from("GENE1"),
            json!({
                "read_details": 1
            }),
        )]);
        let err = parse_parsed_output_from_object(object, None)
            .expect_err("invalid gene payload should fail")
            .to_string();
        assert!(
            err.contains("Failed parsing gene 'GENE1'"),
            "unexpected error: {err}"
        );
        assert!(
            err.contains("read_details should be a JSON object"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn parse_parsed_output_from_object_ok() -> Result<(), DError> {
        let object = serde_json::Map::from_iter([(
            String::from("GENE1"),
            json!({
                "read_details": {
                    "r1": "12"
                }
            }),
        )]);
        let parsed = parse_parsed_output_from_object(object.clone(), None)?;
        assert_eq!(parsed.object, object);
        assert_eq!(parsed.gene_data.len(), 1);
        assert!(parsed.gene_data.contains_key("GENE1"));
        Ok(())
    }
}
