use crate::toolkit::util::DError;

use std::path::PathBuf;

use super::Error;

/// Load a JSON value from a path.
///
/// Supports plain `.json` files and `.xz` compressed JSON files.
/// Compressed inputs are streamed through `xz -dc`.
///
/// # Errors
/// Returns an error if the file can't be read, decompressed, or parsed.
pub(super) fn load_json_value_from_path(
    path: impl Into<PathBuf>,
) -> Result<serde_json::Value, DError> {
    let path = path.into();
    if path
        .extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case("xz"))
    {
        let mut child = std::process::Command::new("xz")
            .arg("-dc")
            .arg(&path)
            .stdout(std::process::Stdio::piped())
            .spawn()?;
        let stdout = child
            .stdout
            .take()
            .ok_or_else(|| Error::new("Failed to capture xz stdout"))?;

        let parsed = serde_json::from_reader(std::io::BufReader::new(stdout))?;
        let status = child.wait()?;
        if !status.success() {
            return Err(Error::new(format!("xz -dc failed with status {status}")).into());
        }
        Ok(parsed)
    } else {
        Ok(serde_json::from_reader(std::io::BufReader::new(
            std::fs::File::open(path)?,
        ))?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temp_json_path(stem: &str) -> PathBuf {
        let mut p = std::env::temp_dir();
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos();
        p.push(format!(
            "paraphase-{stem}-{}-{nonce}.json",
            std::process::id()
        ));
        p
    }

    #[test]
    fn load_json_value_from_path_plain_ok() -> Result<(), DError> {
        let path = temp_json_path("source-ok");
        std::fs::write(&path, br#"{"GENE":{"read_details":{"r1":"12"}}}"#)?;
        let value = load_json_value_from_path(&path)?;
        assert!(value.get("GENE").is_some(), "unexpected value: {value:?}");
        let _ = std::fs::remove_file(path);
        Ok(())
    }

    #[test]
    fn load_json_value_from_path_plain_bad_json() -> Result<(), DError> {
        let path = temp_json_path("source-bad");
        std::fs::write(&path, b"{")?;
        let err = load_json_value_from_path(&path)
            .expect_err("invalid json should fail")
            .to_string();
        assert!(
            err.contains("EOF") || err.contains("eof"),
            "unexpected error: {err}"
        );
        let _ = std::fs::remove_file(path);
        Ok(())
    }
}
