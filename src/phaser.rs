use crate::config::{
    self, region::MatchMap, Gene as GeneConfig, Locus as LocusConfig, Region as RegionConfig,
};
use crate::depth::{Result as GenomeDepthResult, Sex};
use crate::io::json::{GeneCall, ReadFingerprintMap};
use crate::realign::RealignSettings;
use crate::toolkit::deletion::{BigDeletionSettings, Datum as DeletionDatum};
use crate::toolkit::low_complexity::LowConfidenceSites;
use crate::toolkit::range;
use crate::toolkit::site_selection::{CandidateSite, Settings as SiteSelectionSettings};
use crate::toolkit::util::DError;

use itertools::Itertools;
use vstr::VString;

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Debug;
use std::path::PathBuf;
use std::str::FromStr;

pub(crate) const MIN_BASE_QUALITY: u8 = 25u8;
pub type Exception = simple_error::SimpleError;

#[path = "phaser/big_deletions.rs"]
mod big_deletions;
#[path = "phaser/core_utils.rs"]
mod core_utils;
#[path = "phaser/depth_based_operations.rs"]
mod depth_based_operations;
#[path = "phaser/hap_support.rs"]
mod hap_support;
#[path = "phaser/hap_variants.rs"]
mod hap_variants;
#[path = "phaser/io_paths.rs"]
mod io_paths;
#[path = "phaser/read_abstraction.rs"]
mod read_abstraction;
#[path = "phaser/region_depth.rs"]
pub(crate) mod region_depth;
#[path = "phaser/report_call.rs"]
mod report_call;
#[path = "phaser/variant_sites.rs"]
mod variant_sites;
#[path = "phaser/workflow.rs"]
mod workflow;

pub(crate) use core_utils::base_qual;
pub use hap_support::ReadSupport;
pub use hap_variants::{get_start_end, Assignment, Genotype, HapInfo, HapInfoForJson};
pub use io_paths::{build_faidx, FaiBuildError};
pub use read_abstraction::{check_del, fiveprime_clip_length, threeprime_clip_length};

type HapToReads = BTreeMap<VString, Vec<String>>;
type ReadToPossibleHaps = BTreeMap<String, Vec<VString>>;

/// Tuple containing assemblies, read assignments, raw read haps, and counts for each haplotype.
#[derive(Debug, Default, Clone)]
pub struct PhasedResult {
    pub assemblies: crate::assembly::assembly_result::AssemblyResult, // final haps, main haps, hcn
    pub uniquely_supporting_reads: HapToReads,
    pub nonuniquely_supporting_reads: ReadToPossibleHaps,
    pub raw_read_haps: ReadFingerprintMap,
    pub read_counts: (BTreeMap<VString, i32>, BTreeMap<VString, f64>),
}

/// `PhasedResult`, but converted to easily-displayed types for JSON output.
#[derive(Debug, Default, Clone, serde::Deserialize, serde::Serialize)]
pub struct PhasedResultForJson {
    pub assemblies: crate::assembly::assembly_result::AssemblyResultForJson, // final haps, main haps, hcn
    pub uniquely_supporting_reads: BTreeMap<String, Vec<String>>,
    pub nonuniquely_supporting_reads: BTreeMap<String, Vec<String>>,
    pub raw_read_haps: BTreeMap<String, String>,
    pub read_counts: (BTreeMap<String, i32>, BTreeMap<String, f64>),
}

impl PhasedResultForJson {
    /// Generate a struct which can be printed via `serde` from a `PhasedResult` object.
    #[must_use]
    pub fn new(x: &PhasedResult) -> Self {
        let assemblies =
            crate::assembly::assembly_result::AssemblyResultForJson::new(&x.assemblies);
        let uniquely_supporting_reads = x
            .uniquely_supporting_reads
            .iter()
            .map(|(hap, reads)| (hap.to_string(), reads.clone()))
            .collect::<BTreeMap<String, _>>();
        let nonuniquely_supporting_reads = x
            .nonuniquely_supporting_reads
            .iter()
            .map(|(read, haps)| {
                (
                    read.clone(),
                    haps.iter()
                        .map(std::string::ToString::to_string)
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<BTreeMap<String, _>>();
        let raw_read_haps = x
            .raw_read_haps
            .iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect::<BTreeMap<String, _>>();
        let read_counts = (
            x.read_counts
                .0
                .iter()
                .map(|x| (x.0.to_string(), *x.1))
                .collect::<BTreeMap<String, i32>>(),
            x.read_counts
                .1
                .iter()
                .map(|x| (x.0.to_string(), *x.1))
                .collect::<BTreeMap<String, f64>>(),
        );
        Self {
            assemblies,
            uniquely_supporting_reads,
            nonuniquely_supporting_reads,
            raw_read_haps,
            read_counts,
        }
    }
}

/// Sample-level settings
///
/// # Members
/// `genome_reference` - `PathBuf`. Path to genome ref
/// `genome_bam` - `PathBuf`. Path to genome-aligned, indexed bam.
/// `outdir` - `PathBuf`. Path to output directory. Must be writeable to store output bams/vcfs.
/// `sample_id`: `String`,
/// `min_base_quality`: `u8`,
/// `depth`: `Option<GenomeDepthResult>`,
/// `sample_sex`: `Sex`,
/// `homopolymer_window_size`: `Option<usize>`,
/// `big_deletion_settings`: `BigDeletionSettings`,
/// `max_number_deletions`: `u32`,
/// `site_selection_settings`: `site_selection::Settings`,

#[derive(Debug, Clone)]
pub struct Settings {
    pub genome_reference: PathBuf,
    pub genome_bam: PathBuf, // Path to sorted, indexed genome-aligned bam
    pub gene_name: String,
    pub outdir: PathBuf,
    pub sample_id: String,
    pub min_base_quality: u8,
    pub depth: Option<GenomeDepthResult>,
    pub sample_sex: Sex,
    pub homopolymer_window_size: Option<usize>,
    pub big_deletion_settings: BigDeletionSettings,
    pub site_selection_settings: SiteSelectionSettings,
    pub region_config: RegionConfig,
    pub max_number_deletions: u32,
    pub allow_low_coverage: bool,
    pub min_hap_support: i32,
    pub genome: String,
    pub min_variant_frequency: Option<f64>,
    pub min_haplotype_frequency: f64,
    pub targeted: bool,
}

impl Settings {
    ///
    /// Construct settings from sample id, gene name, and paths to aligned bam + reference + output directory.
    ///
    /// Also requires a `paraphase::config::region::Config`, from which site selection parameters are extracted.
    ///
    /// # Arguments
    /// `sample_id`: String representation for sample.
    /// `genome_data`: Tuple (`genome_ref`, `genome_bam`) for genomic reference + reads aligned to this reference.
    /// `outdir`: Directory for storing outputs + intermediate results.
    /// `gene_name`: gene name as a string. Must be a key in the `region_config` file.
    /// `depth`: genomic depth.
    /// `sex`: sample sex. `Sex::{Male, Female, Other}`.
    pub fn new(
        sample_id: impl Into<String>,
        genome_data: (impl Into<PathBuf>, impl Into<PathBuf>),
        outdir: impl Into<PathBuf>,
        gene_name: impl Into<String>,
        region_config: &config::region::Config,
        depth: Option<GenomeDepthResult>,
        sample_sex: Option<Sex>,
        genome: impl Into<String>,
        min_variant_frequency: Option<f64>,
        min_haplotype_frequency: f64,
        targeted: bool,
    ) -> Self {
        let sample_id = sample_id.into();
        let genome = genome.into();
        let (genome_reference, genome_bam) = genome_data;
        let genome_reference = genome_reference.into();
        let genome_bam = genome_bam.into();
        let outdir = outdir.into();
        let gene_name = gene_name.into();
        let sample_sex = sample_sex.unwrap_or(Sex::Other);
        let site_selection_settings = region_config
            .get(&gene_name)
            .map(SiteSelectionSettings::new_from_settings)
            .unwrap_or_else(|| {
                log::warn!(
                    "Gene name '{}' missing in region config; using default site-selection settings.",
                    gene_name
                );
                SiteSelectionSettings::default()
            });
        let region_config = region_config.clone();
        Self {
            genome_reference,
            genome_bam,
            gene_name,
            outdir,
            sample_id,
            depth,
            sample_sex,
            site_selection_settings,
            region_config,
            genome,
            min_variant_frequency,
            min_haplotype_frequency,
            targeted,
            ..Default::default()
        }
    }
}

impl std::default::Default for Settings {
    fn default() -> Self {
        Self {
            sample_sex: Sex::Other,
            min_base_quality: MIN_BASE_QUALITY,
            sample_id: String::from("Unknown"),
            genome_reference: PathBuf::default(),
            genome_bam: PathBuf::default(),
            gene_name: String::default(),
            outdir: PathBuf::default(),
            depth: None,
            homopolymer_window_size: None,
            max_number_deletions: 2,
            big_deletion_settings: BigDeletionSettings::default(),
            site_selection_settings: SiteSelectionSettings::default(),
            allow_low_coverage: false,
            min_hap_support: 4,
            region_config: config::Region::default(),
            genome: String::from("38"),
            min_variant_frequency: None,
            min_haplotype_frequency: 0.03,
            targeted: false,
        }
    }
}

/// Bit flags for Phaser, which lets us pack 4 bools into one u8.
#[repr(u8)]
pub enum FlagBits {
    ToPhase = 1,
    UseSupplementary = 2,
    IsReverse = 4,
    ExpectCN2 = 8,
}

/// Merged config for both locus and genes.
#[derive(Debug, Clone)]
pub struct MergedConfig {
    pub locus: LocusConfig,
    pub gene: GeneConfig,
}

/// Core Phaser
#[derive(Debug)]
pub struct Phaser {
    // Main settings: reference + aligned bam paths, output folder.
    pub settings: Settings,
    pub realign_settings: RealignSettings,

    // Flags
    pub flag: u8,

    // Analysis parameters for specific gene
    pub config: MergedConfig,

    pub realign_region: String,

    pub left_boundary: Option<i64>,
    pub right_boundary: Option<i64>,

    pub gene_start: Option<i64>,
    pub gene_end: Option<i64>,

    pub add_sites: Vec<CandidateSite>,
    pub clip_3p_positions: Vec<i64>,
    pub clip_5p_positions: Vec<i64>,
    pub noisy_regions: Vec<range::I64>,
    pub pivot_site: Option<i64>,

    // Runtime details
    pub low_complexity_sites: LowConfidenceSites,
    pub het_sites: Vec<CandidateSite>,
    pub init_het_sites: Vec<CandidateSite>,
    pub het_sites_no_phasing: Vec<CandidateSite>,
    pub hom_sites: Vec<CandidateSite>,
    pub candidate_sites: BTreeSet<CandidateSite>,
    pub matches: MatchMap, // match in paraphase - maps between coordinates which match between gene and pseudogene.
    pub region_avg_depth: Vec<(f32, f32)>,

    // Deletion regions: usually 0-2 // Could be parameterized later, but using defaults for now.
    pub del_data: Vec<DeletionDatum>,
}

impl Phaser {
    /// Creates a Phaser structure from settings, gene name, and region config.
    pub fn new(
        settings: Settings,
        gene_config: Option<GeneConfig>,
        site_selection_settings: Option<SiteSelectionSettings>,
        realign_settings: Option<RealignSettings>,
    ) -> Result<Self, DError> {
        Self::try_new(
            settings,
            gene_config,
            site_selection_settings,
            realign_settings,
        )
    }

    /// Fallible new construction of Phaser.
    pub fn try_new(
        settings: Settings,
        gene_config: Option<GeneConfig>,
        site_selection_settings: Option<SiteSelectionSettings>,
        realign_settings: Option<RealignSettings>,
    ) -> Result<Self, DError> {
        try_new_impl(
            settings,
            gene_config,
            site_selection_settings,
            realign_settings,
        )
    }

    fn run_gene_specific_workflows(&mut self) -> Option<Result<GeneCall, DError>> {
        match self.settings.gene_name.to_ascii_lowercase().as_str() {
            "cfc1" => Some(self.run_cfc1()),
            "hba" => Some(self.run_hba()),
            "strc" => Some(self.run_strc()),
            "ncf1" => Some(self.run_ncf1()),
            "pms2" => Some(self.run_pms2()),
            "neb" => Some(self.run_neb()),
            "ikbkg" => Some(self.run_ikbkg()),
            "f8" => Some(self.run_f8()),
            "opn1lw" => Some(self.run_opn1lw()),
            "rccx" => Some(self.run_rccx()),
            "smn1" => Some(self.run_smn1()),
            _ => None,
        }
    }

    pub fn run(&mut self) -> Result<GeneCall, DError> {
        log::trace!(
            "Starting Phaser run: gene={}, sample={}, settings={:?}",
            self.gene_name(),
            self.sample_id(),
            self.settings
        );
        if let Some(gene_specific_result) = self.run_gene_specific_workflows() {
            return gene_specific_result;
        }
        self.run_default_gene_pipeline()
    }
}

fn try_new_impl(
    mut settings: Settings,
    gene_config: Option<GeneConfig>,
    site_selection_settings: Option<SiteSelectionSettings>,
    realign_settings: Option<RealignSettings>,
) -> Result<Phaser, DError> {
    let mut locus_config = settings
        .region_config
        .get(&settings.gene_name)
        .or_else(|| settings.region_config.get(&settings.gene_name))
        .ok_or(Exception::new(format!(
            "Failed to get gene config for region {} from region config {:?}",
            settings.gene_name, settings.region_config
        )))?
        .clone();
    if let Some(mut site_settings) = site_selection_settings {
        site_settings.update_from_settings(&locus_config);
        settings.site_selection_settings = site_settings;
    }
    let gene_config = gene_config.unwrap_or_default();
    let realign_settings = realign_settings
        .unwrap_or_default()
        .update_from_locus(&locus_config);
    locus_config.insert("gene_name".into(), settings.gene_name.clone().into());
    let realign_region = if settings.genome != "37" {
        locus_config
            .get("realign_region")
            .and_then(|x| x.as_str())
            .ok_or_else(|| {
                Exception::new(format!(
                    "Missing realign_region in config for gene '{}'",
                    settings.gene_name
                ))
            })?
            .to_string()
    } else {
        locus_config
            .get("realign_region")
            .and_then(|x| x.as_str())
            .ok_or_else(|| {
                Exception::new(format!(
                    "Missing realign_region in config for gene '{}'",
                    settings.gene_name
                ))
            })?
            .strip_prefix("chr")
            .ok_or_else(|| {
                Exception::new(format!(
                    "Expected realign_region to start with 'chr' for genome 37, found '{:?}'",
                    locus_config
                        .get("realign_region")
                        .and_then(serde_yaml::Value::as_str)
                ))
            })?
            .to_string()
    };

    let get_int_field = |key: &str| -> Option<i64> {
        locus_config.get(key).and_then(|x| {
            x.as_str()
                .and_then(|s| s.parse::<i64>().ok())
                .or(x.as_i64())
        })
    };

    let get_sorted_int_vec_field = |key: &str| -> Vec<i64> {
        locus_config
            .get(key)
            .and_then(|x| x.as_sequence())
            .map(|x| {
                x.iter()
                    .filter_map(|x| x.as_i64().map(|v| v - 1))
                    .sorted()
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default()
    };

    let (left_boundary, right_boundary, gene_start, gene_end) =
        ["left_boundary", "right_boundary", "gene_start", "gene_end"]
            .into_iter()
            .map(get_int_field)
            .next_tuple()
            .ok_or(Exception::new(
                "The impossible happened: wrong number of fields in a compile-time array.",
            ))?;
    let pivot_site = get_int_field("pivot_site");
    let add_sites = locus_config
        .get("add_sites")
        .and_then(|x| x.as_sequence())
        .map(|x| {
            x.iter()
                .filter_map(|x| x.as_str().map(CandidateSite::from_str))
                .flatten()
                .collect::<Vec<_>>()
        })
        .unwrap_or_default();
    let clip_3p_positions = get_sorted_int_vec_field("clip_3p_positions");
    let clip_5p_positions = get_sorted_int_vec_field("clip_5p_positions");
    assert!(clip_3p_positions.windows(2).all(|x| x[1] >= x[0]));
    assert!(clip_5p_positions.windows(2).all(|x| x[1] >= x[0]));

    let expect_cn2 = locus_config
        .get("expect_cn2")
        .and_then(serde_yaml::Value::as_bool)
        .unwrap_or(false);
    let is_reverse = locus_config
        .get("is_reverse")
        .and_then(serde_yaml::Value::as_bool)
        .unwrap_or(false);
    let use_supplementary = locus_config.use_supplementary();
    let to_phase = ["in_tandem", "to_phase"]
        .map(String::from)
        .iter()
        .any(|x| locus_config.contains_key(x));
    let flag = [to_phase, use_supplementary, is_reverse, expect_cn2]
        .into_iter()
        .zip([
            FlagBits::ToPhase,
            FlagBits::UseSupplementary,
            FlagBits::IsReverse,
            FlagBits::ExpectCN2,
        ])
        .fold(0u8, |mut acc, x| {
            let (set, value) = x;
            if set {
                acc |= value as u8;
            }
            acc
        });

    let noisy_regions = locus_config
        .try_extract_noisy_regions()
        .map_err(|e| Exception::new(format!("Failed to parse noisy_region: {e}")))?;
    let locus = locus_config;
    let gene = gene_config;

    let mut ret = Phaser {
        settings: settings.clone(),
        realign_settings,
        config: MergedConfig { locus, gene },
        realign_region,
        left_boundary,
        right_boundary,
        add_sites,
        clip_3p_positions,
        clip_5p_positions,
        gene_start,
        gene_end,
        pivot_site,
        candidate_sites: BTreeSet::new(),
        del_data: vec![],
        het_sites: vec![],
        init_het_sites: vec![],
        het_sites_no_phasing: vec![],
        hom_sites: vec![],
        low_complexity_sites: Default::default(),
        matches: Default::default(),
        noisy_regions,
        flag,
        region_avg_depth: vec![],
    };

    ret.matches =
        if let Some(gene2_region) = ret.locus_config().gene2_region(settings.genome == "37") {
            let faidx = ret.make_faidx()?;
            config::region::GenePositionCorrelation::from_ref_regions(
                &ret.realign_region,
                gene2_region,
                &faidx,
                ret.locus_config().chain_bandwidth(),
            )?
            .into_inner()
        } else {
            Default::default()
        };

    Ok(ret)
}
