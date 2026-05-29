use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

type SchemaError = simple_error::SimpleError;

/// Parse YAML bytes and require a mapping at the root.
///
/// # Errors
/// Returns an error if YAML parsing fails or if the root is not a mapping.
pub(crate) fn parse_yaml_root_mapping(
    data: &[u8],
    context: &str,
) -> Result<serde_yaml::Mapping, SchemaError> {
    let value = serde_yaml::from_slice::<serde_yaml::Value>(data)
        .map_err(|e| SchemaError::new(format!("{context}: could not parse yaml: {e}")))?;
    let serde_yaml::Value::Mapping(mapping) = value else {
        return Err(SchemaError::new(format!(
            "{context}: yaml root should be a mapping"
        )));
    };
    Ok(mapping)
}

/// Extract a string key from a YAML mapping key node.
///
/// # Errors
/// Returns an error if the key is not a string.
pub(crate) fn mapping_key_str<'a>(
    key: &'a serde_yaml::Value,
    context: &str,
) -> Result<&'a str, SchemaError> {
    key.as_str().ok_or_else(|| {
        SchemaError::new(format!(
            "{context}: mapping key should be a string, found {key:?}"
        ))
    })
}

/// Require that a YAML value is a sequence.
///
/// # Errors
/// Returns an error if `value` is not an array-like YAML sequence.
pub(crate) fn value_as_sequence<'a>(
    value: &'a serde_yaml::Value,
    path: &str,
) -> Result<&'a serde_yaml::Sequence, SchemaError> {
    value
        .as_sequence()
        .ok_or_else(|| SchemaError::new(format!("{path}: expected an array, found {value:?}")))
}

/// Convert a YAML sequence into a set of strings.
///
/// # Errors
/// Returns an error if any element in the sequence is not a string.
pub(crate) fn sequence_to_string_set(
    values: &serde_yaml::Sequence,
    path: &str,
) -> Result<BTreeSet<String>, SchemaError> {
    values
        .iter()
        .enumerate()
        .map(|(idx, value)| {
            value.as_str().map(str::to_string).ok_or_else(|| {
                SchemaError::new(format!("{path}[{idx}]: expected a string, found {value:?}"))
            })
        })
        .collect::<Result<BTreeSet<_>, SchemaError>>()
}

/// Convert a YAML sequence into an ordered list of strings.
///
/// # Errors
/// Returns an error if any element in the sequence is not a string.
pub(crate) fn sequence_to_string_vec(
    values: &serde_yaml::Sequence,
    path: &str,
) -> Result<Vec<String>, SchemaError> {
    values
        .iter()
        .enumerate()
        .map(|(idx, value)| {
            value.as_str().map(str::to_string).ok_or_else(|| {
                SchemaError::new(format!("{path}[{idx}]: expected a string, found {value:?}"))
            })
        })
        .collect::<Result<Vec<_>, SchemaError>>()
}

/// Convert a YAML mapping into a `BTreeMap<String, Value>`.
///
/// # Errors
/// Returns an error if any key in the mapping is not a string.
pub(crate) fn mapping_to_string_keyed_map(
    mapping: &serde_yaml::Mapping,
    context: &str,
) -> Result<BTreeMap<String, serde_yaml::Value>, SchemaError> {
    let mut out = BTreeMap::new();
    for (key, value) in mapping {
        let key = mapping_key_str(key, context)?;
        out.insert(key.to_string(), value.clone());
    }
    Ok(out)
}

/// Parse YAML where the root is a mapping whose values are themselves mappings
/// and convert both mapping levels to string-keyed `BTreeMap`s.
///
/// # Errors
/// Returns an error if root parsing fails, root is not a mapping, any key is not
/// a string, or any top-level value is not a mapping.
pub(crate) fn parse_yaml_nested_string_keyed_mappings(
    data: &[u8],
    context: &str,
) -> Result<BTreeMap<String, BTreeMap<String, serde_yaml::Value>>, SchemaError> {
    let root = parse_yaml_root_mapping(data, context)?;
    let mut out = BTreeMap::new();
    for (key, value) in root {
        let key = mapping_key_str(&key, context)?;
        let inner = value.as_mapping().ok_or_else(|| {
            SchemaError::new(format!(
                "{context}: entry '{key}' should be a mapping, found {value:?}"
            ))
        })?;
        let inner_context = format!("{context} entry '{key}'");
        out.insert(
            key.to_string(),
            mapping_to_string_keyed_map(inner, &inner_context)?,
        );
    }
    Ok(out)
}

/// Shared loader policy for config data:
/// try user-provided input first, fall back to embedded bytes, then to a final default.
///
/// # Behavior
/// - If parsing `input` fails, logs a warning and retries with `embedded`.
/// - If parsing `embedded` also fails, logs a warning and returns `default()`.
pub(crate) fn load_with_embedded_fallback<T, E, Parse, DefaultFn>(
    context: &str,
    input: Option<&[u8]>,
    embedded: &[u8],
    try_load: Parse,
    default: DefaultFn,
) -> T
where
    E: fmt::Display,
    Parse: Fn(Option<&[u8]>) -> Result<T, E>,
    DefaultFn: FnOnce() -> T,
{
    match try_load(input) {
        Ok(value) => value,
        Err(e) => {
            log::warn!("Failed to parse {context}: {e}. Using embedded config fallback.");
            match try_load(Some(embedded)) {
                Ok(value) => value,
                Err(e2) => {
                    log::warn!(
                        "Failed to parse embedded {context}: {e2}. Falling back to default."
                    );
                    default()
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_yaml_root_mapping_rejects_non_mapping() {
        let err = parse_yaml_root_mapping(b"- a\n- b\n", "cfg")
            .unwrap_err()
            .to_string();
        assert!(err.contains("yaml root should be a mapping"));
    }

    #[test]
    fn sequence_to_string_set_rejects_non_string() {
        let seq = serde_yaml::from_str::<serde_yaml::Value>("[\"a\", 1]")
            .unwrap()
            .as_sequence()
            .unwrap()
            .clone();
        let err = sequence_to_string_set(&seq, "field")
            .unwrap_err()
            .to_string();
        assert!(err.contains("field[1]"));
    }

    #[test]
    fn sequence_to_string_vec_rejects_non_string() {
        let seq = serde_yaml::from_str::<serde_yaml::Value>("[\"a\", 1]")
            .unwrap()
            .as_sequence()
            .unwrap()
            .clone();
        let err = sequence_to_string_vec(&seq, "field")
            .unwrap_err()
            .to_string();
        assert!(err.contains("field[1]"));
    }

    #[test]
    fn parse_yaml_nested_string_keyed_mappings_rejects_non_mapping_entry() {
        let err = parse_yaml_nested_string_keyed_mappings(b"foo: 1\n", "region")
            .unwrap_err()
            .to_string();
        assert!(err.contains("foo"));
        assert!(err.contains("mapping"));
    }
}
