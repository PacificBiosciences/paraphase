use std::collections::BTreeSet;

use crate::config::schema;
use crate::toolkit::util::DError;

/// Configuration specifying logic by gene name.
#[derive(Clone, Debug, Default)]
pub struct Config {
    pub no_genome_depth_genes: BTreeSet<String>, // Don't check genomic depth when calling
    pub no_vcf_genes: BTreeSet<String>,          // Don't write a vcf
    pub genes_to_call: BTreeSet<String>,         // Specify which genes to analyze in this run
    pub check_sex_genes: BTreeSet<String>,       // Which genes require checking sample sex.
    pub two_reference_regions_genes: BTreeSet<String>, // Which genes have multiple reference regions to resolve.
    pub priority_genes: Vec<String>,                   // Which genes should be submitted first.
}

/// Internal decode representation for `gene::Config`.
#[derive(Debug, Default)]
struct RawGeneConfig {
    no_genome_depth_genes: BTreeSet<String>,
    no_vcf_genes: BTreeSet<String>,
    genes_to_call: BTreeSet<String>,
    check_sex_genes: BTreeSet<String>,
    two_reference_regions_genes: BTreeSet<String>,
    priority_genes: Vec<String>,
    genome_depth_genes: BTreeSet<String>,
}

impl RawGeneConfig {
    /// Decode `RawGeneConfig` directly from YAML bytes.
    ///
    /// # Errors
    /// Returns an error if the YAML schema is invalid (wrong root type, unknown
    /// fields, or wrong field element types).
    fn from_yaml(data: &[u8]) -> Result<Self, DError> {
        let mapping = schema::parse_yaml_root_mapping(data, "gene config")?;
        let mut raw = Self::default();
        for (key, value) in mapping {
            let k = schema::mapping_key_str(&key, "gene config")?;
            let path = format!("Gene config field '{k}'");
            let v = schema::value_as_sequence(&value, &path)?;
            match k {
                "genes_to_call" => raw.genes_to_call = schema::sequence_to_string_set(v, k)?,
                "check_sex_genes" => raw.check_sex_genes = schema::sequence_to_string_set(v, k)?,
                "no_vcf_genes" => raw.no_vcf_genes = schema::sequence_to_string_set(v, k)?,
                "genome_depth_genes" => {
                    raw.genome_depth_genes = schema::sequence_to_string_set(v, k)?
                }
                "no_genome_depth_genes" => {
                    raw.no_genome_depth_genes = schema::sequence_to_string_set(v, k)?
                }
                "two_reference_regions_genes" => {
                    raw.two_reference_regions_genes = schema::sequence_to_string_set(v, k)?
                }
                "priority_genes" => raw.priority_genes = schema::sequence_to_string_vec(v, k)?,
                _ => return Err(format!("Unexpected field {k} in gene config").into()),
            }
        }
        Ok(raw)
    }
}

impl From<&RawGeneConfig> for Config {
    fn from(raw: &RawGeneConfig) -> Self {
        Self {
            no_genome_depth_genes: raw.no_genome_depth_genes.clone(),
            no_vcf_genes: raw.no_vcf_genes.clone(),
            genes_to_call: raw.genes_to_call.clone(),
            check_sex_genes: raw.check_sex_genes.clone(),
            two_reference_regions_genes: raw.two_reference_regions_genes.clone(),
            priority_genes: raw.priority_genes.clone(),
        }
    }
}

impl Config {
    /// Builds `Config` from a file location.
    ///
    /// # Errors
    ///
    /// Returns `Err(config::Error)` on error.
    ///
    /// Malformatted yaml.
    /// Missing file at given path or not a file.
    pub fn try_from_path(x: impl Into<std::path::PathBuf>) -> Result<Self, DError> {
        Self::try_load(Some(&std::fs::read(x.into())?))
    }

    /// Parse a `gene::Config` from a path.
    /// Used to provide specific configurations.
    ///
    #[must_use]
    pub fn from_path(x: impl Into<std::path::PathBuf>) -> Self {
        let path = x.into();
        Self::try_from_path(path.clone()).unwrap_or_else(|e| {
            log::warn!(
                "Failed to load gene config from path '{}': {}. Using embedded config fallback.",
                path.display(),
                e
            );
            Self::load(None)
        })
    }

    /// Parse config from optional YAML bytes.
    ///
    /// Uses embedded defaults when `data` is `None`.
    ///
    /// # Errors
    /// Returns an error when YAML/schema validation fails.
    pub fn try_load(data: Option<&[u8]>) -> Result<Self, DError> {
        load(data)
    }

    /// Parse config with fallback to embedded defaults on error.
    #[must_use]
    pub fn load(data: Option<&[u8]>) -> Self {
        const DATA: &[u8] =
            std::include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/genes.yaml"));
        schema::load_with_embedded_fallback(
            "gene config",
            data,
            DATA,
            Self::try_load,
            Self::default,
        )
    }

    /// Whether a gene should use background (genome-wide) depth correction.
    #[must_use]
    pub fn uses_depth(&self, x: &str) -> bool {
        !self.no_genome_depth_genes.contains(x)
    }
}

/// Parse gene confug from `&[u8]` slice.
/// To read from a file, load the file or mmap it.
/// If `None` is provided, use a genes yaml embedded in the executable.
///
/// # Errors
/// If the slice in invalid yaml.
///
/// # Panics
/// • If the parsed yaml is not a dictionary at the root level.
/// • If the keys are not string.
/// • If the values are not arrays of strings.
/// • If there is an unexpected key in the dictionary.
fn load(data: Option<&[u8]>) -> Result<Config, DError> {
    const DATA: &[u8] =
        std::include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/genes.yaml"));
    let data = data.unwrap_or(DATA);
    let raw = RawGeneConfig::from_yaml(data)?;
    let mut ret = Config::from(&raw);

    let all_genes_to_call = if ret.genes_to_call.is_empty() {
        [
            &ret.check_sex_genes,
            &ret.no_vcf_genes,
            &ret.two_reference_regions_genes,
        ]
        .iter()
        .flat_map(|x| x.iter())
        .cloned()
        .collect::<BTreeSet<_>>()
    } else {
        ret.genes_to_call.clone()
    };

    // Account for old paraphase using genome_depth_genes and new paraphase using no_genome_depth_genes.
    // Use set subtraction to get parity.
    if !raw.genome_depth_genes.is_empty() {
        ret.no_genome_depth_genes = all_genes_to_call
            .difference(&raw.genome_depth_genes)
            .cloned()
            .collect::<_>();
    }

    Ok(ret)
}

#[cfg(test)]
mod tests {
    #[test]
    fn gene_config_rejects_non_string_array_items() {
        let data = br#"
no_vcf_genes:
  - CFH
  - 42
"#;
        let err = super::load(Some(data)).unwrap_err().to_string();
        assert!(err.contains("no_vcf_genes[1]"));
    }

    #[test]
    fn gene_config_preserves_priority_gene_order() {
        let data = br#"
priority_genes:
  - FAM86B1
  - ANKRD20A1
  - CLEC18C
"#;
        let config = super::load(Some(data)).expect("priority genes should parse");
        assert_eq!(
            config.priority_genes,
            vec![
                String::from("FAM86B1"),
                String::from("ANKRD20A1"),
                String::from("CLEC18C")
            ]
        );
    }
}
