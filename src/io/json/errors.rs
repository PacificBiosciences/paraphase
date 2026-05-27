use super::Error;

pub(super) fn missing_required(field: &str) -> Error {
    Error::new(format!("Missing required field '{field}'"))
}

pub(super) fn should_be_found_debug(
    field: &str,
    expected: &str,
    found: &serde_json::Value,
) -> Error {
    Error::new(format!("{field} should be {expected}, found {found:?}"))
}

pub(super) fn expected_json_object(context: &str, found: &serde_json::Value) -> Error {
    Error::new(format!("{context}: expected JSON object, found {found:?}"))
}

pub(super) fn not_dictionary_for(type_name: &str) -> Error {
    Error::new(format!("Not a dictionary for {type_name}"))
}
