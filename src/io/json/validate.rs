use crate::toolkit::util::DError;

use super::errors::expected_json_object;

/// Require that `value` is a JSON object and return it.
///
/// # Errors
/// Returns an error if `value` is not an object.
pub(super) fn expect_object<'a>(
    value: &'a serde_json::Value,
    context: &str,
) -> Result<&'a serde_json::Map<String, serde_json::Value>, DError> {
    let serde_json::Value::Object(dict) = value else {
        return Err(expected_json_object(context, value).into());
    };
    Ok(dict)
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn expect_object_ok() -> Result<(), DError> {
        let value = json!({"k": 1});
        let dict = expect_object(&value, "test")?;
        assert_eq!(dict.get("k"), Some(&json!(1)));
        Ok(())
    }

    #[test]
    fn expect_object_err_mentions_context() {
        let value = json!(["not", "an", "object"]);
        let err = expect_object(&value, "validate::tests")
            .expect_err("non-object should fail")
            .to_string();
        assert!(
            err.contains("validate::tests: expected JSON object"),
            "unexpected error: {err}"
        );
    }
}
