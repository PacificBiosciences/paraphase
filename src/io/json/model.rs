use crate::assembly::assembly_result::AssembledPaths;
use crate::config::Region as RegionConfig;
use crate::toolkit::util::DError;

use vstr::VString;

use std::collections::BTreeMap;

use super::decode;
use super::source;
use super::transform;
use super::{JsonMap, PhasingSite, ReadFingerprintMap};

/// This struct contains the inputs to `VariantGraph` (`read_to_hap`),
/// and, if present, the assembled haps in result json being parsed.
#[derive(Debug, Clone, Default)]
pub struct ParaphaseOutput {
    pub pivot_site: Option<i64>,
    pub gene: String,
    pub read_to_hap: ReadFingerprintMap,
    pub assembled_haps: Option<AssembledPaths>,
    pub final_haps: Option<BTreeMap<VString, String>>,
    pub sites_for_phasing: Vec<PhasingSite>,
}

impl ParaphaseOutput {
    #[must_use]
    /// Resolve the pivot site index within `sites_for_phasing`.
    ///
    /// Uses the provided `pivot_site` when set, otherwise falls back to the
    /// parsed `self.pivot_site`. Returns `-1` when no matching site is found.
    pub fn pivot_index(&self, pivot_site: Option<i64>) -> i64 {
        let pivot_site = pivot_site.or(self.pivot_site);
        pivot_site
            .and_then(|pos| {
                self.sites_for_phasing
                    .iter()
                    .position(|x| x.pos == pos)
                    .map(|x| x as i64)
            })
            .unwrap_or(-1)
    }

    #[must_use]
    /// Create a minimal `ParaphaseOutput` with default metadata.
    pub fn new(read_to_hap: ReadFingerprintMap) -> Self {
        Self {
            read_to_hap,
            gene: String::from("NoGene"),
            ..Default::default()
        }
    }

    /// Parse a `ParaphaseOutput` from a gene-level JSON value.
    pub fn from_json(
        value: &serde_json::Value,
        gene: String,
        config: Option<&RegionConfig>,
    ) -> Result<Self, DError> {
        decode::parse_paraphase_output_from_json(value, gene, config)
    }

    /// Read a `ParaphaseOutput` struct from a `std::io::Read` object.
    ///
    /// # Errors
    /// 1. Malformatted json - could not parse.
    /// 2. Malformatted paraphase config - expected a json dictionary object.
    pub fn from_reader(
        reader: impl std::io::Read,
        config: Option<&RegionConfig>,
    ) -> Result<Self, DError> {
        Self::from_json(
            &serde_json::from_reader(reader)?,
            String::from("NoGene"),
            config,
        )
    }

    /// Read a `ParaphaseOutput` struct from a `&std::path::Path` object.
    /// Calls `from_reader` after opening.
    ///
    /// # Errors
    /// 1. Malformatted json - could not parse.
    /// 2. Malformatted paraphase config - expected a json dictionary object.
    pub fn from_path(
        path: &std::path::Path,
        config: Option<&RegionConfig>,
    ) -> Result<Self, DError> {
        Self::from_reader(std::io::BufReader::new(std::fs::File::open(path)?), config)
    }
}

/// Yields the json object from which it was parsed
/// as well as a map from the gene names to the assembly inputs.
/// The `ParaphaseOutput` can then be fed to `VariantGraph` for assembly.
#[derive(Debug, Clone)]
pub struct ParsedParaphaseOutputJSON {
    pub object: serde_json::Map<String, serde_json::Value>,
    pub gene_data: BTreeMap<String, ParaphaseOutput>,
}

impl ParsedParaphaseOutputJSON {
    /// Parse wrapped Paraphase JSON from a reader.
    ///
    /// # Errors
    /// 1. Mal-formatted json.
    pub fn from_reader(
        reader: impl std::io::Read,
        config: Option<&RegionConfig>,
    ) -> Result<Self, DError> {
        Self::from_json(serde_json::from_reader(reader)?, config)
    }

    /// Parse wrapped Paraphase JSON from a path.
    ///
    /// Supports plain JSON and `.xz`-compressed JSON inputs.
    ///
    /// # Errors
    /// 1. Fail to read from input path.
    /// 2. If ends with `.xz`, fail to decompress with `xz` executable.
    /// 3. Mal-formatted json.
    pub fn from_path(
        path: impl Into<std::path::PathBuf>,
        config: Option<&RegionConfig>,
    ) -> Result<Self, DError> {
        Self::from_json(source::load_json_value_from_path(path)?, config)
    }

    /// Parse wrapped Paraphase JSON from a generic JSON value.
    ///
    /// # Errors
    /// 1. Unexpected json format - expected a dictionary type but found something else.
    pub fn from_json(
        object: serde_json::Value,
        config: Option<&RegionConfig>,
    ) -> Result<Self, DError> {
        transform::parse_parsed_output_from_json(object, config)
    }

    /// From `serde_json::Object`, which is a dictionary, build a specific config.
    ///
    /// # Errors
    /// 1. Unexpected json format.
    pub fn from_object(object: JsonMap, config: Option<&RegionConfig>) -> Result<Self, DError> {
        transform::parse_parsed_output_from_object(object, config)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn paraphase_output_pivot_index_uses_arg_then_default() {
        let out = ParaphaseOutput {
            pivot_site: Some(100),
            gene: String::from("GENE"),
            read_to_hap: ReadFingerprintMap::new(),
            assembled_haps: None,
            final_haps: None,
            sites_for_phasing: vec![
                PhasingSite::try_new("100_A_T").expect("valid site"),
                PhasingSite::try_new("200_G_C").expect("valid site"),
            ],
        };
        assert_eq!(out.pivot_index(None), 0);
        assert_eq!(out.pivot_index(Some(200)), 1);
        assert_eq!(out.pivot_index(Some(300)), -1);
    }

    #[test]
    fn paraphase_output_from_reader_ok() -> Result<(), DError> {
        let raw = br#"{
            "read_details": {
                "r1": "121x"
            },
            "sites_for_phasing": ["100_A_T"]
        }"#;
        let out = ParaphaseOutput::from_reader(std::io::Cursor::new(raw.as_slice()), None)?;
        assert_eq!(out.gene, "NoGene");
        assert_eq!(out.read_to_hap.len(), 1);
        assert_eq!(out.sites_for_phasing.len(), 1);
        Ok(())
    }

    #[test]
    fn parsed_paraphase_output_from_json_ok() -> Result<(), DError> {
        let root = json!({
            "GENE1": {
                "read_details": {
                    "r1": "12"
                }
            }
        });
        let parsed = ParsedParaphaseOutputJSON::from_json(root, None)?;
        assert_eq!(parsed.gene_data.len(), 1);
        assert!(parsed.gene_data.contains_key("GENE1"));
        Ok(())
    }
}
