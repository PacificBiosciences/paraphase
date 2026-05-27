//! Shared VCF data structures.
//!
//! These types define the intermediate representation used by VCF collection
//! and merging (haplotype boundaries, per-position calls, and I/O tuple helpers).

use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};

#[derive(Clone, Debug)]
pub struct HapBoundForVcf {
    /// Output haplotype/sample name for this boundary record.
    pub hap_name: String,
    /// 0-based
    pub start: i64,
    /// use only unique reads for variant calling between start and start_strict
    pub start_strict: i64,
    pub end: i64,
    /// use only unique reads for variant calling between end and end_strict
    pub end_strict: i64,
    pub is_truncated: Vec<String>,
}

#[derive(Clone, Debug)]
pub struct HapVariantInfo {
    /// Per-position small-variant observations aligned to `hap_info` ordering.
    pub(crate) variants_info: BTreeMap<i64, Vec<Option<VariantInfoByHP>>>,
    /// Per-position symbolic SV observations aligned to `hap_info` ordering.
    pub(crate) sv_variants: BTreeMap<i64, Vec<Option<String>>>,
    /// Ordered haplotype boundary metadata used for VCF sample columns.
    pub(crate) hap_info: Vec<HapBoundForVcf>,
}

/// Input/output tuple for per-region BAM processing used by VCF writing.
pub struct IOTuple(pub PathBuf, pub PathBuf, pub String, pub bool);

impl IOTuple {
    #[must_use]
    /// Source BAM path to read alignments from.
    pub fn source_bam(&self) -> &Path {
        &self.0
    }
    #[must_use]
    /// Destination BAM path for rewritten/filtered records.
    pub fn dest_bam(&self) -> &Path {
        &self.1
    }
    #[must_use]
    /// Chromosome/contig name associated with this I/O task.
    pub fn chromosome_name(&self) -> &str {
        &self.2
    }
    #[must_use]
    /// Whether this tuple refers to secondary-region (`gene2`) processing.
    pub fn is_gene2(&self) -> bool {
        self.3
    }
}

impl std::fmt::Debug for IOTuple {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "IOTuple{{Source: {:?}. Dest: {:?}. Chrom: {:?}. Primary or secondary gene: {}",
            self.source_bam(),
            self.dest_bam(),
            self.chromosome_name(),
            if self.is_gene2() {
                "secondary"
            } else {
                "primary"
            }
        )
    }
}

#[derive(Clone, Debug)]
/// Variant at a position from a set of reads
pub struct VariantInfoByHP {
    /// variant bases
    pub base: Option<String>,
    /// reference bases
    pub ref_base: String,
    /// total depth
    pub depth: usize,
    /// number of reads supporting the variant
    pub nread: usize,
    /// original consensus call format (with + or - for indels)
    pub original_base: Option<String>,
    pub bases_count: HashMap<Vec<u8>, usize>,
}
