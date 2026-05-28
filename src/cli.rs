use crate::toolkit::util::FULL_VERSION;
use chrono::Datelike;
use clap::Parser;
use std::path::{Path, PathBuf};

// Clap does not like boxed errors, use strings just in cli
pub type Result<T> = std::result::Result<T, String>;

#[derive(Parser, Clone)]
#[clap(
    name="paraphase",
    author,
    version=&**FULL_VERSION,
    long_about = None,
    disable_help_subcommand = true,
    about,
    after_help = format!("Copyright (C) 2004-{}     Pacific Biosciences of California, Inc.
This program comes with ABSOLUTELY NO WARRANTY; it is intended for
Research Use Only and not for use in diagnostic procedures.", chrono::Utc::now().year()),
help_template = "{name} {version}\n{author}\n{about-section}\n{usage-heading}\n    {usage}\n\n{all-args}{after-help}",
)]
pub struct Settings {
    #[clap(required = true)]
    #[clap(short = 'b')]
    #[clap(long = "bam")]
    #[clap(help = "BAM file with aligned HiFi reads")]
    #[clap(value_name = "BAM")]
    #[arg(value_parser = check_file_exists)]
    pub bam: PathBuf,

    #[clap(required = true)]
    #[clap(short = 'r')]
    #[clap(long = "reference")]
    #[clap(help = "Path to reference genome FASTA file")]
    #[clap(value_name = "FASTA")]
    #[arg(value_parser = check_file_exists)]
    pub reference: PathBuf,

    #[clap(required = true)]
    #[clap(short = 'o')]
    #[clap(long = "out")]
    #[clap(help = "Output directory")]
    #[clap(value_name = "outdir")]
    pub outdir: PathBuf,

    #[clap(short, long)]
    #[clap(help = "Prefix of output files for a single sample.\n\
If not provided, prefix is extracted from the header of the input BAM.")]
    pub prefix: Option<String>,

    #[clap(long, short, default_value = "")]
    #[clap(
        help = "Optionally specify which regions(s) to run (separated by comma).\n\
If not provided, all regions are run.\n\
The full set of accepted regions are defined in the config file."
    )]
    pub gene: String,

    #[clap(long, short)]
    #[clap(
        help = "Optional path to a user-defined config file listing the full set of regions to analyze.\n\
By default paraphase uses the config file in data/38/config.yaml."
    )]
    #[arg(value_parser = check_file_exists)]
    pub config: Option<PathBuf>,

    #[clap(long = "genome")]
    #[clap(
        help = "Optionally specify which genome reference build the input BAM files are aligned against.\n\
Accepted values are 19, 37, chm13, and 38."
    )]
    #[clap(default_value = "38")]
    pub genome: String,

    #[clap(short = 't')]
    #[clap(long = "threads")]
    #[clap(help = "Number of threads to use. Defaults to number available.")]
    #[clap(value_name = "THREADS")]
    #[arg(value_parser = threads_in_range)]
    pub num_threads: Option<usize>,

    #[clap(long, action)]
    #[clap(help = "If specified, paraphase will not assume depth is uniform across the genome.")]
    pub targeted: bool,

    #[clap(long)]
    #[clap(
        help = "Minimum frequency for a variant to be used for phasing. Works with targeted mode.\n\
The cutoff for variant-supporting reads is max(5, total_depth * min_frequency).\n\
total_depth is the combined depth of all paralogs in a paralog group.\n\
Default: 0.11."
    )]
    pub min_variant_frequency: Option<f64>,

    #[clap(long)]
    #[clap(
        help = "Minimum frequency of unique supporting reads for a haplotype. Works with targeted mode.\n\
The cutoff for haplotype-supporting reads is max(4, total_depth * min_frequency).\n\
total_depth is the combined depth of all paralogs in a paralog group."
    )]
    #[clap(default_value = "0.03")]
    pub min_haplotype_frequency: f64,

    #[clap(long, action)]
    #[clap(help = "If specified, paraphase will not write VCFs.")]
    pub novcf: bool,

    #[clap(long, action)]
    #[clap(
        help = "If specified, paraphase will write no-call sites in the VCFs, marked with LowQual filter."
    )]
    pub write_nocalls_in_vcf: bool,

    #[clap(long, action)]
    #[clap(
        help = "If specified, variant calls are made against the main gene only.\n\
By default, for SMN1, PMS2, STRC, NCF1, and IKBKG, haplotypes are assigned to gene or\n\
paralog/pseudogene, and variants are called against gene or paralog/pseudogene, respectively."
    )]
    pub gene1only: bool,

    #[cfg(feature = "pprof")]
    #[clap(help_heading("Advanced"))]
    #[clap(long)]
    #[clap(help = "Profiling sample frequency in Hz (samples/second). If unset, no profiling.")]
    pub profile_freq: Option<i32>,

    #[clap(help_heading("Advanced"))]
    #[clap(short = 'v')]
    #[clap(long = "verbose")]
    #[clap(action = clap::ArgAction::Count)]
    #[clap(help = "Verbose output.\n\
`-v` enables debug-level logs.\n\
`-vv` (or higher) enables trace-level logs.")]
    pub verbosity: u8,

    #[clap(help_heading("Advanced"))]
    #[clap(short = 'q')]
    #[clap(long = "quiet", action)]
    #[clap(help = "Quiet output (errors only). If set, this overrides `-v` / `--verbose`.")]
    pub quiet: bool,
}

impl Settings {
    /// Resolve effective log level from `-v` / `--verbose`.
    pub fn log_level(&self) -> log::LevelFilter {
        use log::LevelFilter;
        if self.quiet {
            return LevelFilter::Error;
        }
        match self.verbosity {
            0 => LevelFilter::Info,
            1 => LevelFilter::Debug,
            _ => LevelFilter::Trace,
        }
    }

    /// Normalize CLI settings after parse (currently fills default thread count).
    pub fn check(mut self) -> Self {
        let default_threads = std::thread::available_parallelism()
            .map(std::num::NonZeroUsize::get)
            .unwrap_or(1);
        self.num_threads = Some(self.num_threads.unwrap_or(default_threads));
        self
    }
}

/// Initialize process-wide logger from CLI settings.
pub fn setup_logging(settings: &Settings) {
    env_logger::builder()
        .format_timestamp_millis()
        .filter_level(settings.log_level())
        .init();
}

/// Clap value parser: require an existing path.
fn check_file_exists(s: &str) -> Result<PathBuf> {
    let path = Path::new(s);
    if !path.exists() {
        Err(format!("File does not exist: {}", path.display()))
    } else {
        Ok(path.to_path_buf())
    }
}

/// Clap value parser: parse and validate `threads >= 1`.
fn threads_in_range(s: &str) -> Result<usize> {
    let thread: usize = s
        .parse()
        .map_err(|_| format!("`{}` is not a valid thread number", s))?;
    if thread >= 1 {
        Ok(thread)
    } else {
        Err("Number of threads must be at least 1".into())
    }
}
