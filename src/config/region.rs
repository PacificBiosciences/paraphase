use crate::config::schema;
use crate::toolkit::range::I64 as Range64;
use crate::toolkit::util::{consumes_qry, consumes_ref, DError};

use itertools::Itertools;
use minimap2::ffi;
use rust_htslib::{
    bam::record::{Cigar, CigarString},
    faidx,
};

type ConfigError = simple_error::SimpleError;

use std::collections::BTreeMap;

/// This config corresponds to `${build}/config.yaml`.
/// It is a key-value store from gene/region ID to analysis parameters.
#[derive(Debug, Clone)]
pub struct Config(pub BTreeMap<String, Locus>);

/// This config corresponds to a gene-specific portion of `${build}/config.yaml`.
/// It is a key-value store from feature to parameter.
#[derive(Debug, Clone)]
pub struct Locus(pub BTreeMap<String, serde_yaml::Value>);

/// Internal decode representation for `region::Config`.
#[derive(Debug, Clone)]
struct RawRegionConfig(BTreeMap<String, RawLocus>);

/// Internal decode representation for a single locus entry.
#[derive(Debug, Clone)]
struct RawLocus(BTreeMap<String, serde_yaml::Value>);

impl RawRegionConfig {
    /// Decode `RawRegionConfig` directly from YAML bytes.
    ///
    /// # Errors
    /// Returns an error if the YAML schema is invalid (wrong root type, wrong
    /// nested shape, or non-string mapping keys).
    fn from_yaml(data: &[u8]) -> Result<Self, ConfigError> {
        let data = schema::parse_yaml_nested_string_keyed_mappings(data, "region config")
            .map_err(|e| ConfigError::new(e.to_string()))?;
        Ok(Self(
            data.into_iter()
                .map(|(k, v)| (k, RawLocus(v)))
                .collect::<BTreeMap<_, _>>(),
        ))
    }
}

impl From<RawRegionConfig> for Config {
    fn from(value: RawRegionConfig) -> Self {
        Config(
            value
                .0
                .into_iter()
                .map(|(k, v)| (k, Locus(v.0)))
                .collect::<BTreeMap<_, _>>(),
        )
    }
}

impl Config {
    /// Builds Config from a raw text slice.
    #[must_use]
    pub fn from_slice(x: &[u8]) -> Self {
        load(Some(x))
    }

    /// Builds `Config` from a raw text slice.
    #[must_use]
    pub fn load(x: Option<&[u8]>) -> Self {
        try_load(x).unwrap_or_else(|e| {
            log::warn!("Failed to parse region config: {e}. Using built-in default fallback.");
            load(None)
        })
    }

    /// Builds `Config` from a raw text slice.
    /// # Errors
    /// Errors on malformatted yaml.
    pub fn try_load(x: Option<&[u8]>) -> Result<Self, ConfigError> {
        try_load(x)
    }

    /// Builds `Config` from a file location.
    ///
    /// # Errors
    ///
    /// Returns `Err(config::Error)` on error.
    ///
    /// Malformatted yaml.
    /// Missing file at given path or not a file.
    pub fn try_from_path(x: impl Into<std::path::PathBuf>) -> Result<Self, DError> {
        Ok(try_load(Some(&std::fs::read(x.into())?))?)
    }

    #[must_use]
    pub fn try_from_slice(x: &[u8]) -> Self {
        try_load(Some(x)).unwrap_or_else(|e| {
            log::warn!(
                "Failed to parse yaml from provided slice ({:?} bytes): {e:?}. Using built-in default fallback.",
                x.len()
            );
            load(None)
        })
    }

    /// Parse a `region::Config` from a path.
    /// Used to provide specific configurations.
    /// We also bake in hg19 and hg38 configurations into the binary for convenience.
    #[must_use]
    pub fn from_path(x: impl Into<std::path::PathBuf>) -> Self {
        let path = x.into();
        Self::try_from_path(path.clone()).unwrap_or_else(|e| {
            log::warn!(
                "Failed to read region config from path '{}': {}. Using built-in default fallback.",
                path.display(),
                e
            );
            load(None)
        })
    }
}

impl std::default::Default for Config {
    fn default() -> Self {
        load(None)
    }
}

impl std::ops::Deref for Config {
    type Target = BTreeMap<String, Locus>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::ops::DerefMut for Config {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl std::ops::Deref for Locus {
    type Target = BTreeMap<String, serde_yaml::Value>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::ops::DerefMut for Locus {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

// From minimap2-rs
fn cigar_to_cigarstr(cigar: &[(u32, u8)]) -> CigarString {
    CigarString(
        cigar
            .iter()
            .map(|(len, op)| match op {
                0 => Cigar::Match(*len),
                1 => Cigar::Ins(*len),
                2 => Cigar::Del(*len),
                3 => Cigar::RefSkip(*len),
                4 => Cigar::SoftClip(*len),
                5 => Cigar::HardClip(*len),
                6 => Cigar::Pad(*len),
                7 => Cigar::Equal(*len),
                8 => Cigar::Diff(*len),
                _ => {
                    log::warn!("Unexpected cigar operation code {op}; using soft-clip fallback.");
                    Cigar::SoftClip(*len)
                }
            })
            .collect::<Vec<Cigar>>(),
    )
}

///
/// Loads stats from a byte slice.
///
/// # Panics
///
///
/// Panics if the first level of keys in the dictionary are not all string types.
/// These keys should all be gene/region names.
///
/// Panics if the second level of keys in the dictionary are not all string types.
///
/// # Errors
/// Throws a `ConfigError` if the config.yaml was not a dictionary.
pub fn try_load(data: Option<&[u8]>) -> Result<Config, ConfigError> {
    const DATA: &[u8] =
        std::include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/38/config.yaml"));
    let data = data.unwrap_or(DATA);
    Ok(RawRegionConfig::from_yaml(data)?.into())
}

///
#[must_use]
pub fn load(data: Option<&[u8]>) -> Config {
    const DATA: &[u8] =
        std::include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/38/config.yaml"));
    schema::load_with_embedded_fallback("region config", data, DATA, try_load, || {
        Config(BTreeMap::new())
    })
}

#[cfg(test)]
mod tests {
    #[test]
    fn region_config_rejects_non_mapping_locus() {
        let data = br#"
smn1: "not a mapping"
"#;
        let err = super::try_load(Some(data)).unwrap_err().to_string();
        assert!(err.contains("smn1"));
        assert!(err.contains("mapping"));
    }
}

impl Locus {
    #[must_use]
    /// Build a locus wrapper from raw key/value YAML content.
    pub fn new(data: BTreeMap<String, serde_yaml::Value>) -> Self {
        Self(data)
    }

    /// Extract the noisy regions from a json.
    /// # Panics
    /// 1. malformatted noisy regions - should be an array of arrays of integers.
    /// 2. `noisy_region` field should be a sequence/array.
    ///
    /// ```
    /// use paraphase::toolkit::range::I64 as Range64;
    /// let mut mapping: std::collections::BTreeMap<String, serde_yaml::Value> =
    ///     std::collections::BTreeMap::new();
    /// let mut regions =
    ///     serde_yaml::from_str("[[100, 200], [500, 600]]").expect("Failed to parse json");
    /// mapping.insert("noisy_region".into(), regions);
    /// let mut config = paraphase::config::Locus::new(mapping);
    /// let regions = config.extract_noisy_regions();
    /// assert_eq!(regions.len(), 2);
    /// assert_eq!(regions[0], Range64::new(99, 200));
    /// assert_eq!(regions[1], Range64::new(499, 600));
    ///
    /// //  Panics if regions are not a sequence of sequences.
    /// let mut mapping: std::collections::BTreeMap<String, serde_yaml::Value> =
    ///     std::collections::BTreeMap::new();
    /// let mut regions = serde_yaml::from_str("[100, 200]").expect("Failed to parse json");
    /// mapping.insert("noisy_region".into(), regions);
    /// let mut config = paraphase::config::Locus::new(mapping.clone());
    /// let result = std::panic::catch_unwind(|| config.extract_noisy_regions());
    /// assert!(result.is_err());
    ///
    /// let mut regions = serde_yaml::from_str("hello world").expect("Failed to parse json");
    /// mapping.insert("noisy_region".into(), regions);
    /// let mut config = paraphase::config::Locus::new(mapping);
    /// let result = std::panic::catch_unwind(|| config.extract_noisy_regions());
    /// assert!(result.is_err());
    /// ```
    #[must_use]
    pub fn extract_noisy_regions(&self) -> Vec<Range64> {
        self.try_extract_noisy_regions()
            .unwrap_or_else(|e| panic!("{e}"))
    }

    /// Fallibly parse the optional `noisy_region` field.
    ///
    /// # Errors
    /// 1. If `noisy_region` exists but is not a list of `[start, stop]` integer pairs.
    pub fn try_extract_noisy_regions(&self) -> Result<Vec<Range64>, DError> {
        let Some(noisy_region) = self.get("noisy_region") else {
            return Ok(Vec::new());
        };
        let intervals = noisy_region
            .as_sequence()
            .ok_or_else(|| ConfigError::new("noisy_region should be a sequence"))?;
        intervals
            .iter()
            .enumerate()
            .map(|(idx, item)| {
                let interval = item.as_sequence().ok_or_else(|| {
                    ConfigError::new(format!(
                        "noisy_region[{idx}] should be a sequence of length 2, found {item:?}"
                    ))
                })?;
                if interval.len() != 2 {
                    return Err(ConfigError::new(format!(
                        "noisy_region[{idx}] should have length 2, found {}",
                        interval.len()
                    ))
                    .into());
                }
                let start = interval[0].as_i64().ok_or_else(|| {
                    ConfigError::new(format!(
                        "noisy_region[{idx}][0] should be an integer, found {:?}",
                        interval[0]
                    ))
                })?;
                let stop = interval[1].as_i64().ok_or_else(|| {
                    ConfigError::new(format!(
                        "noisy_region[{idx}][1] should be an integer, found {:?}",
                        interval[1]
                    ))
                })?;
                Ok(Range64::new(start - 1, stop))
            })
            .collect()
    }

    /// Access `extract_regions` field from input yaml.
    /// This identifies which regions to extract from.
    #[must_use]
    pub fn extract_regions(&self, strip_chr: bool) -> Vec<String> {
        if !strip_chr {
            self.get("extract_regions")
                .and_then(|x| x.as_str())
                .map(|region| {
                    region
                        .split_whitespace()
                        .map(std::string::ToString::to_string)
                        .collect::<Vec<_>>()
                })
                .unwrap_or_default()
        } else {
            self.get("extract_regions")
                .and_then(|x| x.as_str())
                .map(|region| {
                    region
                        .split_whitespace()
                        .map(std::string::ToString::to_string)
                        .map(|a| a.strip_prefix("chr").unwrap_or(&a).to_string())
                        .collect::<Vec<_>>()
                })
                .unwrap_or_default()
        }
    }

    /// Sets the chain bandwidth for alignment.
    /// The default is 500 (minimap2's default), but if `use_r2k` is set in the config, we use 2000 instead.
    #[must_use]
    pub fn chain_bandwidth(&self) -> Option<i32> {
        self.get("use_r2k").map(|_x| 2000)
    }

    /// This function generates a `Vec<String>` of extra arguments to pass to minimap2 CLI calls.
    /// It is generally not used, as we have switched to the `minimap2-rs` crate.
    #[must_use]
    pub fn extra_mm2_options(&self) -> Vec<String> {
        let mut ret = vec![];
        if self.get("use_r2k").is_some() {
            ret.push("-r2k".into());
        }
        ret
    }

    /// Access `use_supplementary` field from input yaml.
    /// This determines if a alignments for a given region which are supplementary should be considered.
    #[must_use]
    pub fn use_supplementary(&self) -> bool {
        self.get("is_tandem")
            .or(self.get("use_supplementary"))
            .is_some()
    }

    /// Whether to phase hapltypes into alleles
    #[must_use]
    pub fn to_phase(&self) -> bool {
        self.get("is_tandem").or(self.get("to_phase")).is_some()
    }

    /// Whether to call fusion
    #[must_use]
    pub fn call_fusion(&self) -> Option<&str> {
        self.get("call_fusion").and_then(|x| x.as_str())
    }

    /// Access `gene2_region` field from input yaml.
    #[must_use]
    pub fn gene2_region(&self, strip_chr: bool) -> Option<&str> {
        let gene2 = self.get("gene2_region").and_then(|x| x.as_str())?;
        let gene2 = if strip_chr {
            gene2.strip_prefix("chr")?
        } else {
            gene2
        };
        Some(gene2)
    }

    /// Access `position_match` field from input yaml.
    #[must_use]
    pub fn position_match(&self) -> Option<&str> {
        self.get("gene_position_match").and_then(|x| x.as_str())
    }

    /// Parse optional depth-normalization intervals from `depth_region`.
    ///
    /// Intervals are returned as-is in config coordinate space.
    #[must_use]
    pub fn depth_region(&self) -> Vec<Range64> {
        self.get("depth_region")
            .map(|x| {
                if let Some(seq) = x.as_sequence() {
                    seq.iter()
                        .filter_map(|x| {
                            let interval = x.as_sequence()?;
                            let (start, stop) = interval.iter().next_tuple()?;
                            let start = start.as_i64()?;
                            let stop = stop.as_i64()?;
                            Some(Range64::new(start, stop))
                        })
                        .collect::<Vec<_>>()
                } else {
                    Vec::new()
                }
            })
            .unwrap_or_default()
    }
} // impl Locus

pub type MatchMap = BTreeMap<i64, i64>;

/// We use 0-based coordinates here instead of 1-based Paraphase.
#[derive(Clone, PartialEq, Debug)]
pub struct GenePositionCorrelation {
    pub matches: MatchMap,
}

impl GenePositionCorrelation {
    /// Build from a cached text file.
    ///```
    /// use paraphase::config::region::GenePositionCorrelation;
    /// use paraphase::toolkit::util::test_file;
    /// let input = test_file("position_match_output.txt");
    /// let correlation = GenePositionCorrelation::from_text_file(&input).unwrap();
    /// assert_eq!(
    ///     correlation.matches.get(&(70963975 - 1)),
    ///     Some(&(70088528 - 1))
    /// );
    /// ```
    pub fn from_text_file(x: impl AsRef<std::path::Path>) -> Result<Self, DError> {
        use std::io::BufRead;
        let input = std::fs::File::open(x.as_ref()).map_err(|e| {
            ConfigError::new(format!(
                "Failed to open file at {:?}. Error: {e:?}",
                x.as_ref()
            ))
        })?;
        let reader = std::io::BufReader::new(input);
        let mut matches = MatchMap::default();
        for line in reader.lines() {
            let line = line?;
            let (from, to) = line
                .split_whitespace()
                .map(str::parse::<i64>)
                .next_tuple()
                .ok_or(ConfigError::new(format!("Mal-formatted line {line}")))?;
            matches.insert(from? - 1, to? - 1); // Convert 1-based to 0-based.
        }
        Ok(Self { matches })
    }

    /// Build a cross-region base-position mapping by aligning two reference windows.
    ///
    /// Coordinates in the returned map are zero-based in both regions.
    ///```
    /// if let Ok(hg38) = std::env::var("HG38") {
    ///     use paraphase::config::region::GenePositionCorrelation;
    ///     use paraphase::toolkit::util::test_file;
    ///
    ///     let faidx = rust_htslib::faidx::Reader::from_path(hg38).unwrap();
    ///     let region1 = "chr5:70890000-71100000"; // SMN1
    ///     let region2 = "chr5:70040526-70088546"; // SMN2
    ///
    ///     // First, check identity.
    ///     let matches =
    ///         GenePositionCorrelation::from_ref_regions(region1, region1, &faidx, None).unwrap();
    ///     assert!(matches.matches.iter().all(|(x, y)| x == y));
    ///
    ///     // Now, check parity with paraphase/minimap2.
    ///     let matches = GenePositionCorrelation::from_ref_regions(region1, region2, &faidx, None)
    ///         .unwrap()
    ///         .into_inner();
    ///     let input = test_file("position_match_output.txt");
    ///
    ///     let from_python = GenePositionCorrelation::from_text_file(&input)
    ///         .unwrap()
    ///         .into_inner();
    ///     let mut failures = Vec::new();
    ///
    ///     // Make sure all python coordinates are in ours.
    ///     assert!(from_python.keys().all(|x| matches.contains_key(x)));
    ///     for (k, v) in from_python.iter() {
    ///         if *v != matches[k] {
    ///             failures.push((*v, matches[v]))
    ///         }
    ///     }
    ///     assert!(failures.is_empty());
    ///
    ///     // Make sure all rust coordinates are in python.
    ///     assert!(matches.keys().all(|x| from_python.contains_key(x)));
    ///     for (k, v) in matches.iter() {
    ///         if *v != from_python[k] {
    ///             failures.push((*v, from_python[v]))
    ///         }
    ///     }
    ///     assert!(failures.is_empty());
    ///     assert_eq!(matches, from_python);
    /// }
    /// ```
    pub fn from_ref_regions(
        region1: &str,
        region2: &str,
        fa: &faidx::Reader,
        chain_bandwidth: Option<i32>,
    ) -> Result<Self, DError> {
        use rust_htslib::bam::record::Cigar::{Diff, Equal, Match};
        let parsed = |x: &str| -> Result<(String, usize, usize), DError> {
            let mut toks = x.split_terminator(':');
            let chr = toks.next().ok_or_else(|| {
                ConfigError::new(format!(
                    "Could not find chromosome in reference region string {x}"
                ))
            })?;
            let (start, end) = toks
                .next()
                .ok_or_else(|| {
                    ConfigError::new("Could not find coordinates in reference region string {x}")
                })?
                .split_terminator('-')
                .next_tuple()
                .ok_or_else(|| {
                    ConfigError::new("Could not find start/stop in reference region string {x}")
                })?;
            let start = start.parse::<i64>()? - 1;
            let end = end.parse::<i64>()? - 1;
            Ok((chr.to_owned(), start as usize, end as usize))
        };
        let parsed1 = parsed(region1)?;
        let parsed2 = parsed(region2)?;
        let seq1 = fa.fetch_seq(&parsed1.0, parsed1.1, parsed1.2)?;
        let seq2 = fa.fetch_seq(&parsed2.0, parsed2.1, parsed2.2)?;
        /*
        let tempdir = tempfile::TempDir::new()?;
        let seqname_1 = itertools::intersperse(region1.split_terminator(|x| matches!(*x, ':' | '-')), "_").collect::<String>();
        let tmpfa = tempdir.path().join(format!("{seqname_1}.fa"));
        */
        let mut matches = MatchMap::default();
        let mut dna_aligner = minimap2::Aligner::builder()
            .with_cigar()
            .with_seq_and_id(&seq1, region1.as_bytes())?;
        dna_aligner.mapopt.flag |= ffi::MM_F_CIGAR as i64;
        dna_aligner.mapopt.flag |= ffi::MM_F_EQX as i64;
        dna_aligner.mapopt.bw = chain_bandwidth.unwrap_or(500);
        let mappings = dna_aligner.map(
            &seq2, /* output_cigar= */ true, /* output_md= */ true,
            /* max_frag_len= */ None,
            /* extra_flags= */ None, //Some(realign::EXTRA_FLAGS_SLICE),
            None,
        )?;
        let get_offset = |x: &str, num_idx: usize| -> Result<u32, ConfigError> {
            Ok(x.split_terminator(':')
                .nth(1)
                .and_then(|x| x.split_terminator('-').nth(num_idx))
                .and_then(|x| x.parse::<u32>().ok())
                .ok_or_else(|| {
                    ConfigError::new(format!(
                        "Failed to get offset from string + num_idx: {x}/{num_idx}"
                    ))
                })?
                - 1)
        };
        let ref_offset = get_offset(region1, 0)?;
        let qry_offset = get_offset(region2, 0)?;
        let qry_end = get_offset(region2, 1)?;
        log::trace!(
            "Coordinate offsets for correlation: ref_offset={ref_offset}, qry_offset={qry_offset}, qry_end={qry_end}"
        );
        for mapping in mappings.iter().filter(|x| x.is_primary) {
            let cigar = cigar_to_cigarstr(
                mapping
                    .alignment
                    .as_ref()
                    .ok_or_else(|| {
                        ConfigError::new(format!("Missing alignment in mapping {mapping:?}"))
                    })?
                    .cigar
                    .as_ref()
                    .ok_or_else(|| {
                        ConfigError::new(format!("Missing cigar in alignment {mapping:?}"))
                    })?,
            );

            let reverse = mapping.strand == minimap2::Strand::Reverse;
            log::trace!("Processing primary mapping; reverse_strand={reverse}");

            if !reverse {
                let mut ref_idx = mapping.target_start as u32;
                let mut qry_idx = mapping.query_start as u32;
                for op in &cigar {
                    let len = op.len();
                    log::trace!("Applying CIGAR op {op:?} at ref_idx={ref_idx}, qry_idx={qry_idx}");
                    match op {
                        Diff(len) | Equal(len) | Match(len) => {
                            for idx in 0..*len {
                                let offset_ref_pos = ref_offset + ref_idx + idx;
                                let offset_qry_pos = qry_offset + qry_idx + idx;
                                matches
                                    .entry(offset_ref_pos.into())
                                    .or_insert_with(|| offset_qry_pos.into());
                                log::trace!(
                                    "Mapped local index {idx}: ref_pos={offset_ref_pos}, qry_pos={offset_qry_pos}"
                                );
                            }
                        }
                        _ => {}
                    };
                    if consumes_ref(*op) {
                        ref_idx += len;
                    }
                    if consumes_qry(*op) {
                        qry_idx += len;
                    }
                    log::trace!(
                        "Completed CIGAR op {op:?}; updated indices ref_idx={ref_idx}, qry_idx={qry_idx}"
                    );
                }
            } else {
                let mut ref_idx = mapping.target_start as u32;
                let mut qry_idx = mapping.query_end as u32;
                for op in &cigar {
                    let len = op.len();
                    log::trace!("Applying CIGAR op {op:?} at ref_idx={ref_idx}, qry_idx={qry_idx}");
                    match op {
                        Diff(len) | Equal(len) | Match(len) => {
                            for idx in 0..*len {
                                let offset_ref_pos = ref_offset + ref_idx + idx;
                                let offset_qry_pos = qry_offset + qry_idx - idx - 1;
                                matches
                                    .entry(offset_ref_pos.into())
                                    .or_insert_with(|| offset_qry_pos.into());
                                log::trace!(
                                    "Mapped local index {idx}: ref_pos={offset_ref_pos}, qry_pos={offset_qry_pos}"
                                );
                            }
                        }
                        _ => {}
                    };
                    if consumes_ref(*op) {
                        ref_idx += len;
                    }
                    if consumes_qry(*op) {
                        qry_idx -= len;
                    }
                    log::trace!(
                        "Completed CIGAR op {op:?}; updated indices ref_idx={ref_idx}, qry_idx={qry_idx}"
                    );
                }
            }
        }

        Ok(Self { matches })
    }

    /// Extract the internal position map from the `GenePositionCorrelation` object.
    /// Consumes the calling object.
    #[must_use]
    pub fn into_inner(self) -> MatchMap {
        self.matches
    }
}
