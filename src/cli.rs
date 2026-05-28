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

    #[clap(long = "genome")]
    #[clap(
        help = "Optionally specify which genome reference build the input BAM files are aligned against. 
        Accepted values are 19, 37, chm13, and 38."
    )]
    #[clap(default_value = "38")]
    pub genome: String,

    #[clap(required = true)]
    #[clap(short = 'o')]
    #[clap(long = "out")]
    #[clap(help = "Output directory")]
    #[clap(value_name = "outdir")]
    pub outdir: PathBuf,

    /// Optional path to a user-defined config file listing the full set of regions to analyze.
    /// By default paraphase uses the config file in data/38/config.yaml
    #[clap(long, short)]
    #[arg(value_parser = check_file_exists)]
    pub config: Option<PathBuf>,

    /// Prefix of output files for a single sample.
    /// If not provided, prefix will be extracted from the header of the input BAM.
    #[clap(short, long)]
    pub prefix: Option<String>,

    /// Optionally specify which region(s) to run (separated by comma).
    /// Will run all regions if not specified.
    /// The full set of accepted regions are defined in the config file.
    #[clap(long, short, default_value = "")]
    pub gene: String,

    #[clap(short = 't')]
    #[clap(long = "threads")]
    #[clap(help = "Number of threads to use. Defaults to number available.")]
    #[clap(value_name = "THREADS")]
    #[arg(value_parser = threads_in_range)]
    pub num_threads: Option<usize>,

    #[clap(long)]
    #[clap(
        help = "Minimum frequency for a variant to be used for phasing. Works with the targeted mode.
        The cutoff for variant-supporting reads is determined by max(5, total_depth * min_frequency).
        Note that total_depth is the combined depth of all paralogs for a paralog group.
        Default is 0.11."
    )]
    pub min_variant_frequency: Option<f64>,

    #[clap(long)]
    #[clap(
        help = "Minimum frequency of unique supporting reads for a haplotype. Works with the targeted mode.
        The cutoff for haplotype-supporting reads is determined by max(4, total_depth * min_frequency).
        Note that total_depth is the combined depth of all paralogs for a paralog group.
        Default is 0.03."
    )]
    #[clap(default_value = "0.03")]
    pub min_haplotype_frequency: f64,

    #[clap(long, action)]
    #[clap(help = "If specified, paraphase will not assume depth is uniform across the genome.")]
    pub targeted: bool,

    #[clap(long, action)]
    #[clap(
        help = "If specified, variant calls will be made against the main gene only.
        By default, for SMN1, PMS2, STRC, NCF1 and IKBKG, haplotypes are assigned to gene or
        paralog/pseudogene, and variants are called against gene or paralog/pseudogene, respectively."
    )]
    pub gene1only: bool,

    #[clap(long, action)]
    #[clap(help = "If specified, paraphase will not write VCFs.")]
    pub novcf: bool,

    #[clap(long, action)]
    #[clap(
        help = "If specified, paraphase will write no-call sites in the VCFs, marked with LowQual filter."
    )]
    pub write_nocalls_in_vcf: bool,

    #[cfg(feature = "pprof")]
    /// Profiling sample frequency in Hz (samples/second). If unset, no profiling.
    #[clap(help_heading("Advanced"))]
    #[clap(long)]
    pub profile_freq: Option<i32>,

    /// Verbose output
    /// Apply once for debug-level messages
    /// Apply twice or more for trace-level messages
    /// If --log-level is set, this option takes precedence.
    /// Equivalence to --log-level:
    ///       => "info"
    ///    -v => "debug"
    ///   -vv => "trace"
    ///  -vvv => "trace"
    #[clap(help_heading("Advanced"))]
    #[clap(short = 'v')]
    #[clap(long = "verbose")]
    #[clap(action = clap::ArgAction::Count, verbatim_doc_comment)]
    #[clap(hide = true)]
    pub verbosity: u8,

    /// Values: "error", "warn", "info" (default), "debug", "trace".
    /// Higher verbosity or log levels emit more detailed diagnostic output.
    #[clap(help_heading("Advanced"))]
    #[clap(long, verbatim_doc_comment, default_value = "info")]
    pub log_level: String,
}

impl Settings {
    /// Resolve effective log level, with `-v`/`--verbose` taking precedence over `--log-level`.
    pub fn log_level(&self) -> log::LevelFilter {
        use log::LevelFilter;
        if self.verbosity > 0 {
            match self.verbosity {
                0 => LevelFilter::Info,
                1 => LevelFilter::Debug,
                _ => LevelFilter::Trace,
            }
        } else {
            match &self.log_level.clone().to_lowercase()[..] {
                "info" => LevelFilter::Info,
                "warn" => LevelFilter::Warn,
                "debug" => LevelFilter::Debug,
                "error" => LevelFilter::Error,
                "trace" => LevelFilter::Trace,
                _ => {
                    log::warn!(
                        "Unsupported log level '{}'; defaulting to 'info'.",
                        self.log_level
                    );
                    LevelFilter::Info
                }
            }
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
