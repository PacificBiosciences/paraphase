use crate::cli::Settings;
use crate::config;
use crate::config::gene::Config as GeneConfig;
use crate::config::region::Config as RegionConfig;
use crate::depth;
use crate::depth::Result as DepthResult;
use crate::depth::Settings as DepthSettings;
use crate::depth::Sex;
use crate::io::bam::BamWriter;
use crate::io::json::{write_outputs, GeneCall};
use crate::io::vcf::VcfWriter;
use crate::phaser;
use crate::toolkit::update_calls::update_calls_after_per_gene_analysis;
use crate::toolkit::util::{self, output_bam_header, DError, DResult};
use rust_htslib::bam::{self, Read};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Condvar, Mutex};

type GeneJobResult = Result<(GeneCall, PathBuf), String>;
type SharedJobResults = Arc<Mutex<Vec<Option<GeneJobResult>>>>;
type IndexedGene = (usize, String);

struct MergeHeapItem {
    tid: i32,
    pos: i64,
    source: usize,
    serial: u64,
    record: bam::Record,
}

impl PartialEq for MergeHeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.tid == other.tid
            && self.pos == other.pos
            && self.source == other.source
            && self.serial == other.serial
    }
}

impl Eq for MergeHeapItem {}

impl PartialOrd for MergeHeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for MergeHeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        (other.tid, other.pos, other.source, other.serial).cmp(&(
            self.tid,
            self.pos,
            self.source,
            self.serial,
        ))
    }
}

/// Source selector for region-config loading.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RegionConfigSource {
    CustomPath,
    Hg19,
    Chm13,
    Default,
}

#[must_use]
/// Fallback sample-name derivation from BAM filename stem.
pub fn sample_from_bam_name(bam: &PathBuf) -> String {
    bam.file_stem()
        .or_else(|| bam.file_name())
        .map(|name| name.to_string_lossy().to_string())
        .unwrap_or_else(|| String::from("sample"))
}

/// Parse the comma-delimited `--gene` argument into a list of requested gene names.
#[must_use]
pub fn parse_requested_genes(gene_arg: &str) -> Vec<String> {
    gene_arg
        .split_terminator(',')
        .map(std::borrow::ToOwned::to_owned)
        .collect::<Vec<_>>()
}

#[must_use]
/// Case-insensitive stable sort for gene names (with original-case tie break).
fn sort_gene_names_case_insensitive<I>(genes: I) -> Vec<String>
where
    I: IntoIterator<Item = String>,
{
    let mut genes = genes.into_iter().collect::<Vec<_>>();
    genes.sort_by_cached_key(|gene| (gene.to_ascii_lowercase(), gene.clone()));
    genes
}

/// Resolve requested genes against the region config and reject all-invalid requests.
///
/// # Errors
/// Returns an error when no valid genes remain after filtering.
pub fn resolve_requested_genes(
    config: &RegionConfig,
    gene_arg: &str,
) -> Result<Vec<String>, DError> {
    let requested = parse_requested_genes(gene_arg);
    let genes = if requested.is_empty() {
        sort_gene_names_case_insensitive(config.keys().cloned())
    } else {
        let mut valid = Vec::new();
        for gene in &requested {
            if config.contains_key(gene) {
                valid.push(gene.to_string());
            } else {
                log::warn!(
                    "Skipping unrecognized gene name '{gene}'; it is not present in the region config."
                );
            }
        }
        valid
    };
    if genes.is_empty() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!(
                "Please provide valid gene name(s). Accepted genes are {}",
                config.keys().cloned().collect::<Vec<_>>().join(",")
            ),
        )
        .into());
    }
    Ok(genes)
}

/// Return whether we should create a VCF output directory for this run.
#[must_use]
pub fn should_create_vcf_dir(
    novcf: bool,
    genes: &[String],
    no_vcf_genes: &BTreeSet<String>,
) -> bool {
    !novcf && genes.iter().any(|gene| !no_vcf_genes.contains(gene))
}

/// Build the VCF output directory path for a sample.
#[must_use]
pub fn vcf_output_dir(outdir: &Path, sample: &str) -> PathBuf {
    outdir.join(format!("{sample}_paraphase_vcfs"))
}

/// Build the tagged BAM output path for a sample.
#[must_use]
pub fn output_bam_path(outdir: &Path, sample: &str) -> PathBuf {
    outdir.join(format!("{sample}.paraphase.bam"))
}

/// Build the JSON output path for a sample.
#[must_use]
pub fn output_json_path(outdir: &Path, sample: &str) -> PathBuf {
    outdir.join(format!("{sample}.paraphase.json"))
}

/// Return whether VCF writing is enabled for one gene.
#[must_use]
pub fn should_write_gene_vcf(novcf: bool, gene: &str, no_vcf_genes: &BTreeSet<String>) -> bool {
    !novcf && !no_vcf_genes.contains(gene)
}

/// Compute genome-wide depth metrics using the CLI-selected genome/background source.
///
/// # Errors
/// Returns an error if depth calculation setup fails (e.g. BAM/bed loading).
pub fn compute_depth_result(args: &Settings) -> Result<DepthResult, DError> {
    let mut calculator = if args.genome == "19" {
        depth::Calculator::from_hg19(&args.bam, None)
    } else if args.genome == "37" {
        depth::Calculator::from_hg19(
            &args.bam,
            Some(DepthSettings {
                strip_chr: true,
                ..Default::default()
            }),
        )
    } else if args.genome == "chm13" {
        depth::Calculator::from_chm13(&args.bam, None)
    } else {
        log::debug!("Using HG38 background depth coordinates.");
        depth::Calculator::from_hg38(&args.bam, None)
    }?;
    Ok(calculator.compute())
}

/// Resolve the output sample name from `--prefix`, BAM read groups, or BAM filename.
#[must_use]
pub fn resolve_sample_name(args: &Settings) -> String {
    args.prefix.clone().unwrap_or_else(|| {
        util::sample_names_from_input(&args.bam)
            .and_then(|names| {
                log::debug!(
                    "Sample names discovered in BAM {}: {names:?}",
                    args.bam.display()
                );
                match names.len() {
                    1 => names.into_keys().next(),
                    _ => None,
                }
            })
            .unwrap_or_else(|| sample_from_bam_name(&args.bam))
    })
}

/// Decide whether genome depth should be used to correct per-gene depth calling.
#[must_use]
pub fn depth_for_correction(depth: DepthResult, targeted: bool) -> Option<DepthResult> {
    if targeted {
        log::info!(
            "Targeted mode is enabled; genome-wide coverage will not be used for depth correction."
        );
        None
    } else if depth.median < 10.0 || depth.median_absolute_difference > 0.25 {
        log::info!(
            "Genome-wide coverage is too low or too variable; skipping depth-based correction."
        );
        None
    } else {
        Some(depth)
    }
}

/// Decide which region configuration source should be used for this run.
#[must_use]
pub fn region_config_source(config_path: Option<&PathBuf>, genome: &str) -> RegionConfigSource {
    if config_path.is_some() {
        RegionConfigSource::CustomPath
    } else if genome == "19" || genome == "37" {
        RegionConfigSource::Hg19
    } else if genome == "chm13" {
        RegionConfigSource::Chm13
    } else {
        RegionConfigSource::Default
    }
}

/// Load region configuration from CLI arguments and genome selection.
///
/// # Errors
/// Returns an error when a user-provided config path cannot be loaded.
pub fn load_region_config(args: &Settings) -> Result<RegionConfig, DError> {
    Ok(
        match region_config_source(args.config.as_ref(), &args.genome) {
            RegionConfigSource::CustomPath => {
                let path = args
                    .config
                    .as_ref()
                    .expect("config path is present for CustomPath source");
                config::Region::try_from_path(path)?
            }
            RegionConfigSource::Hg19 => RegionConfig::load(Some(config::REGION_CONFIG_HG19)),
            RegionConfigSource::Chm13 => RegionConfig::load(Some(config::REGION_CONFIG_CHM13)),
            RegionConfigSource::Default => config::Region::default(),
        },
    )
}

#[must_use]
/// Build per-sample BAM shard directory path.
pub fn bam_shard_dir(outdir: &Path, sample: &str) -> PathBuf {
    outdir.join(format!("{sample}.paraphase.bam_shards"))
}

#[must_use]
/// Build the shard BAM path for one gene under the sample shard directory.
pub fn gene_bam_shard_path(outdir: &Path, sample: &str, gene: &str) -> PathBuf {
    bam_shard_dir(outdir, sample).join(format!("{sample}_{gene}.bam"))
}

/// Log the canonicalized CLI invocation and current working directory.
pub fn log_cli_invocation() {
    log::debug!(
        "CLI:'{}', called from {}",
        std::env::args()
            .map(|x| {
                if let Ok(path) = std::fs::canonicalize(std::path::PathBuf::from(&x)) {
                    path.to_string_lossy().to_string()
                } else {
                    x
                }
            })
            .collect::<Vec<String>>()
            .join(" "),
        std::env::current_dir()
            .map(|x| x.display().to_string())
            .unwrap_or_else(|_| String::from("<unknown_cwd>"))
    );
}

/// Create run output directories (`outdir` and optionally sample VCF directory).
///
/// # Errors
/// Returns an error if any directory cannot be created.
pub fn prepare_output_directories(
    outdir: &Path,
    sample: &str,
    novcf: bool,
    genes: &[String],
    no_vcf_genes: &BTreeSet<String>,
) -> DResult {
    std::fs::create_dir_all(outdir)?;
    std::fs::create_dir_all(bam_shard_dir(outdir, sample))?;
    if should_create_vcf_dir(novcf, genes, no_vcf_genes) {
        let vcf_dir = vcf_output_dir(outdir, sample);
        std::fs::create_dir_all(vcf_dir)?;
    }
    Ok(())
}

/// Write final phasing results JSON file for the sample.
///
/// # Errors
/// Returns an error if creating or writing the JSON output file fails.
pub fn write_json_output(
    phasing_results: &BTreeMap<String, GeneCall>,
    outdir: &Path,
    sample: &str,
) -> DResult {
    let output_path = output_json_path(outdir, sample);
    write_outputs(
        phasing_results,
        &mut std::io::BufWriter::new(std::fs::File::create(output_path)?),
    )?;
    Ok(())
}

/// Run all per-gene jobs in a Rayon pool and capture panics as error strings.
///
/// # Errors
/// Returns an error if the Rayon thread pool cannot be constructed.
#[allow(clippy::type_complexity)]
pub fn run_gene_jobs(
    genes: &[String],
    sample: &str,
    config: &RegionConfig,
    gene_config: &GeneConfig,
    args: &Settings,
    genome_depth: Option<DepthResult>,
    sample_sex: Sex,
) -> Result<Vec<Result<(GeneCall, PathBuf), String>>, DError> {
    let job_slots = args.num_threads.unwrap_or(1);
    log::debug!(
        "Starting per-region analysis with {} worker threads",
        job_slots
    );
    let (priority_jobs, nonpriority_jobs) = schedule_gene_jobs(genes, &gene_config.priority_genes);
    let priority_launch_target = priority_jobs.len().min(job_slots);
    let launch_gate = Arc::new((
        Mutex::new(PriorityLaunchState {
            started: 0,
            target: priority_launch_target,
            open: priority_launch_target == 0,
        }),
        Condvar::new(),
    ));
    let results = Arc::new(Mutex::new(vec![None; genes.len()]));

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(job_slots)
        .build()
        .map_err(|e| std::io::Error::other(format!("failed to build Rayon thread pool: {e}")))?;

    pool.scope_fifo(|scope| {
        for (original_index, gene) in priority_jobs {
            let ctx = JobContext {
                original_index,
                gene,
                sample: sample.to_string(),
                config: config.clone(),
                gene_config: gene_config.clone(),
                args: args.clone(),
                genome_depth,
                sample_sex,
                results: Arc::clone(&results),
                launch_gate: Arc::clone(&launch_gate),
            };
            scope.spawn_fifo(move |_| run_gene_job(ctx, true));
        }

        for (original_index, gene) in nonpriority_jobs {
            let ctx = JobContext {
                original_index,
                gene,
                sample: sample.to_string(),
                config: config.clone(),
                gene_config: gene_config.clone(),
                args: args.clone(),
                genome_depth,
                sample_sex,
                results: Arc::clone(&results),
                launch_gate: Arc::clone(&launch_gate),
            };
            scope.spawn_fifo(move |_| run_gene_job(ctx, false));
        }
    });

    take_job_results(results)
}

struct JobContext {
    original_index: usize,
    gene: String,
    sample: String,
    config: RegionConfig,
    gene_config: GeneConfig,
    args: Settings,
    genome_depth: Option<DepthResult>,
    sample_sex: Sex,
    results: SharedJobResults,
    launch_gate: Arc<(Mutex<PriorityLaunchState>, Condvar)>,
}

struct PriorityLaunchState {
    started: usize,
    target: usize,
    open: bool,
}

fn run_gene_job(ctx: JobContext, is_priority: bool) {
    if is_priority {
        mark_priority_job_started(&ctx.launch_gate);
    } else {
        wait_for_priority_launch(&ctx.launch_gate);
    }

    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        process_gene(
            ctx.gene,
            ctx.sample,
            ctx.config,
            ctx.gene_config,
            ctx.args,
            ctx.genome_depth,
            Some(ctx.sample_sex),
        )
        .map_err(|e| e.to_string())
    }));
    let result = match result {
        Ok(inner) => inner,
        Err(panic_payload) => Err(format!("Panic: {}", panic_payload_message(panic_payload))),
    };

    ctx.results.lock().expect("job results lock poisoned")[ctx.original_index] = Some(result);
}

/// Mark one priority job as launched and open the launch gate when enough
/// priority jobs have started.
fn mark_priority_job_started(launch_gate: &(Mutex<PriorityLaunchState>, Condvar)) {
    let (state, condvar) = launch_gate;
    let mut state = state.lock().expect("priority launch gate poisoned");
    state.started += 1;
    if !state.open && state.started >= state.target {
        state.open = true;
        condvar.notify_all();
    }
}

/// Block a non-priority job until the priority launch gate opens.
fn wait_for_priority_launch(launch_gate: &(Mutex<PriorityLaunchState>, Condvar)) {
    let (state, condvar) = launch_gate;
    let mut state = state.lock().expect("priority launch gate poisoned");
    while !state.open {
        state = condvar
            .wait(state)
            .expect("priority launch gate poisoned while waiting");
    }
}

/// Collect scoped job results and fail if any slot was left unset.
fn take_job_results(results: SharedJobResults) -> Result<Vec<GeneJobResult>, DError> {
    let mut results = results
        .lock()
        .map_err(|_| std::io::Error::other("job results lock poisoned"))?;
    results
        .drain(..)
        .map(|result| result.ok_or_else(|| std::io::Error::other("missing gene job result").into()))
        .collect()
}

/// Reorder genes for job submission so selected heavy genes enter the pool first.
///
/// Priority genes are emitted first in the exact config-defined order. All remaining genes
/// retain their original relative order afterward.
#[must_use]
pub fn schedule_gene_jobs(
    genes: &[String],
    priority_genes: &[String],
) -> (Vec<IndexedGene>, Vec<IndexedGene>) {
    let priority_order = priority_genes
        .iter()
        .enumerate()
        .map(|(idx, gene)| (gene.as_str(), idx))
        .collect::<BTreeMap<_, _>>();
    let mut priority = Vec::new();
    let mut nonpriority = Vec::new();
    for (original_index, gene) in genes.iter().cloned().enumerate() {
        if let Some(priority_index) = priority_order.get(gene.as_str()).copied() {
            priority.push((priority_index, original_index, gene));
        } else {
            nonpriority.push((original_index, gene));
        }
    }
    priority.sort_by_key(|(priority_index, original_index, _)| (*priority_index, *original_index));
    let priority = priority
        .into_iter()
        .map(|(_, original_index, gene)| (original_index, gene))
        .collect::<Vec<_>>();
    (priority, nonpriority)
}

/// Merge per-gene execution results into output structures and write BAM records.
///
/// # Errors
/// Returns an error if BAM record writing fails.
pub fn apply_gene_results(
    genes: &[String],
    results: &[Result<(GeneCall, PathBuf), String>],
    phasing_results: &mut BTreeMap<String, GeneCall>,
    bam_shards: &mut Vec<PathBuf>,
) -> DResult {
    for (gene, result) in genes.iter().zip(results.iter()) {
        match result {
            Ok((gene_call, bam_shard)) => {
                phasing_results.insert(gene_call.gene_name.clone(), gene_call.clone());
                bam_shards.push(bam_shard.clone());
            }
            Err(e) => {
                log::error!("Gene processing failed for {:?}: {}", gene, e);
            }
        }
    }
    Ok(())
}

/// Run the full per-sample analysis workflow from directory setup through final JSON output.
///
/// # Errors
/// Returns an error if any pipeline stage fails, including BAM/VCF/JSON IO and post-processing.
pub fn run_pipeline_for_sample(
    args: &Settings,
    config: &RegionConfig,
    gene_config: &GeneConfig,
    genes: &[String],
    sample: &str,
    genome_depth: Option<DepthResult>,
    sample_sex: Sex,
) -> DResult {
    let mut phasing_results = BTreeMap::<String, _>::new();
    let mut bam_shards = Vec::new();
    prepare_output_directories(
        args.outdir.as_path(),
        sample,
        args.novcf,
        genes,
        &gene_config.no_vcf_genes,
    )?;

    let results = run_gene_jobs(
        genes,
        sample,
        config,
        gene_config,
        args,
        genome_depth,
        sample_sex,
    )?;
    apply_gene_results(genes, &results, &mut phasing_results, &mut bam_shards)?;
    update_calls_after_per_gene_analysis(&mut phasing_results)?;
    merge_bam_shards(args, sample, &bam_shards)?;
    write_json_output(&phasing_results, args.outdir.as_path(), sample)?;
    Ok(())
}

/// Finalize profiling by building and writing the flamegraph when a profiler guard is present.
///
/// # Errors
/// Returns an error if writing `profile.flamegraph.svg` fails.
#[cfg(feature = "pprof")]
pub fn write_profile_flamegraph(guard: Option<pprof::ProfilerGuard<'_>>, outdir: &Path) -> DResult {
    if let Some(guard) = guard {
        match guard.report().build() {
            Ok(report) => {
                let file = std::fs::File::create(outdir.join("profile.flamegraph.svg"))?;
                report.flamegraph(file).map_err(|e| {
                    std::io::Error::other(format!("failed to write profile.flamegraph.svg: {e}"))
                })?;
            }
            Err(e) => {
                log::warn!("Profiling report build failed: {e}");
            }
        }
    }
    Ok(())
}

/// Conduct analysis for one gene and return the per-gene call plus tagged BAM records.
///
/// # Errors
/// Returns an error if phasing, BAM writing, VCF writing, or temp-dir operations fail.
pub fn process_gene(
    gene: String,
    sample: String,
    config: RegionConfig,
    gene_config: GeneConfig,
    args: Settings,
    genome_depth: Option<DepthResult>,
    sex: Option<Sex>,
) -> Result<(GeneCall, PathBuf), DError> {
    log::info!("Starting per-gene analysis for {}.", gene);
    let write_gene_vcf = should_write_gene_vcf(args.novcf, &gene, &gene_config.no_vcf_genes);
    let tmp_dir = tempfile::TempDir::new()?;
    let mut phaser = create_gene_phaser(
        &gene,
        &sample,
        &config,
        gene_config,
        &args,
        genome_depth,
        sex,
        tmp_dir.path(),
    )?;
    let res = phaser.run()?;
    let tagged_bams = collect_gene_bam_paths(&phaser, &res)?;
    let bam_shard = gene_bam_shard_path(args.outdir.as_path(), &sample, &gene);
    merge_tagged_gene_bams(&tagged_bams, &bam_shard)?;
    maybe_write_gene_vcf(&phaser, &res, &args, &sample, write_gene_vcf)?;
    tmp_dir.close()?;
    log::info!("Completed per-gene analysis for {gene}.");
    Ok((res, bam_shard))
}

/// Convert a panic payload into a human-readable message for logs.
#[must_use]
fn panic_payload_message(panic_payload: Box<dyn std::any::Any + Send>) -> String {
    panic_payload
        .downcast::<String>()
        .map(|s| *s)
        .or_else(|p| p.downcast::<&'static str>().map(|s| s.to_string()))
        .unwrap_or_else(|_| "Unknown panic (run with RUST_BACKTRACE=1 for details)".to_string())
}

/// Build a configured `Phaser` instance for one gene.
///
/// # Errors
/// Returns an error if phaser setup fails.
fn create_gene_phaser(
    gene: &str,
    sample: &str,
    config: &RegionConfig,
    gene_config: GeneConfig,
    args: &Settings,
    genome_depth: Option<DepthResult>,
    sex: Option<Sex>,
    tmp_dir: &Path,
) -> Result<phaser::Phaser, DError> {
    let settings = phaser::Settings::new(
        sample,
        (&args.reference, &args.bam),
        tmp_dir,
        gene.to_string(),
        config,
        genome_depth,
        sex,
        &args.genome,
        args.min_variant_frequency,
        args.min_haplotype_frequency,
        args.targeted,
    );

    phaser::Phaser::new(
        settings,
        Some(gene_config),
        None, // Option<SiteSelectionSettings>
        None, // Option<RealignSettings>
    )
}

/// Collect tagged BAM paths from a finished phaser result.
///
/// # Errors
/// Returns an error if BAM writing fails.
fn collect_gene_bam_paths(phaser: &phaser::Phaser, res: &GeneCall) -> Result<Vec<PathBuf>, DError> {
    let bam_writer = BamWriter::new(phaser, res);
    let _written_bams = bam_writer.write_bams()?;
    let mut bam_paths = Vec::new();

    let gene1_tagged = phaser.realigned_tagged_bam_path();
    let gene1_fallback = phaser.realigned_bam_path();
    if gene1_tagged.exists() {
        bam_paths.push(gene1_tagged);
    } else if gene1_fallback.exists() {
        bam_paths.push(gene1_fallback);
    }

    if phaser
        .locus_config()
        .gene2_region(phaser.settings.genome == "37")
        .is_some()
    {
        let gene2_tagged = phaser.realigned_tagged_gene2_bam_path();
        let gene2_fallback = phaser.realigned_gene2_bam_path();
        if gene2_tagged.exists() {
            bam_paths.push(gene2_tagged);
        } else if gene2_fallback.exists() {
            bam_paths.push(gene2_fallback);
        }
    }

    Ok(bam_paths)
}

fn next_bam_record(
    reader: &mut bam::Reader,
) -> Result<Option<bam::Record>, rust_htslib::errors::Error> {
    let mut record = bam::Record::new();
    if let Some(status) = reader.read(&mut record) {
        status?;
        return Ok(Some(record.clone()));
    }
    Ok(None)
}

fn merge_bam_readers_sorted(readers: &mut [bam::Reader], writer: &mut bam::Writer) -> DResult {
    let mut heap = BinaryHeap::new();
    let mut serial = 0u64;
    for (source, reader) in readers.iter_mut().enumerate() {
        if let Some(record) = next_bam_record(reader)? {
            heap.push(MergeHeapItem {
                tid: record.tid(),
                pos: record.pos(),
                source,
                serial,
                record,
            });
            serial += 1;
        }
    }

    while let Some(item) = heap.pop() {
        writer.write(&item.record)?;
        if let Some(record) = next_bam_record(&mut readers[item.source])? {
            heap.push(MergeHeapItem {
                tid: record.tid(),
                pos: record.pos(),
                source: item.source,
                serial,
                record,
            });
            serial += 1;
        }
    }
    Ok(())
}

fn merge_tagged_gene_bams(tagged_bams: &[PathBuf], output_bam: &Path) -> DResult {
    let existing_bams = tagged_bams
        .iter()
        .filter(|path| path.exists())
        .cloned()
        .collect::<Vec<_>>();
    let Some(first_bam) = existing_bams.first() else {
        return Err(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            format!(
                "No tagged BAM shards were created for {}",
                output_bam.display()
            ),
        )
        .into());
    };
    let first_reader = bam::Reader::from_path(first_bam)?;
    let mut writer = bam::Writer::from_path(
        output_bam,
        &output_bam_header(first_reader.header()),
        bam::Format::Bam,
    )?;
    drop(first_reader);
    let mut readers = existing_bams
        .iter()
        .map(bam::Reader::from_path)
        .collect::<Result<Vec<_>, _>>()?;
    merge_bam_readers_sorted(&mut readers, &mut writer)?;
    Ok(())
}

fn merge_bam_shards(args: &Settings, sample: &str, bam_shards: &[PathBuf]) -> DResult {
    let reader = bam::Reader::from_path(&args.bam)?;
    let output_bam = output_bam_path(args.outdir.as_path(), sample);
    let mut writer = bam::Writer::from_path(
        &output_bam,
        &output_bam_header(reader.header()),
        bam::Format::Bam,
    )?;
    let mut readers = bam_shards
        .iter()
        .map(bam::Reader::from_path)
        .collect::<Result<Vec<_>, _>>()?;
    merge_bam_readers_sorted(&mut readers, &mut writer)?;
    drop(writer);
    drop(readers);
    for shard_path in bam_shards {
        if shard_path.exists() {
            std::fs::remove_file(shard_path)?;
        }
    }
    let shard_dir = bam_shard_dir(args.outdir.as_path(), sample);
    if shard_dir.exists() {
        std::fs::remove_dir_all(shard_dir)?;
    }
    bam::index::build(&output_bam, None, bam::index::Type::Bai, 1)?;
    Ok(())
}

/// Write a per-gene VCF when enabled for this run and gene.
///
/// # Errors
/// Returns an error if VCF writing fails.
fn maybe_write_gene_vcf(
    phaser: &phaser::Phaser,
    res: &GeneCall,
    args: &Settings,
    sample: &str,
    enabled: bool,
) -> DResult {
    if enabled {
        let vcf_dir = &vcf_output_dir(&args.outdir, sample);
        let vcf_writer = VcfWriter::new(phaser, res, args.write_nocalls_in_vcf, args.gene1only);
        vcf_writer.write_vcf(vcf_dir)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::region::Locus;
    use crate::depth::Sex;
    use std::collections::BTreeMap as Map;

    #[test]
    fn parse_requested_genes_handles_comma_list() {
        assert_eq!(
            parse_requested_genes("smn1,CFH"),
            vec![String::from("smn1"), String::from("CFH")]
        );
    }

    #[test]
    fn region_config_source_selects_expected_source() {
        assert_eq!(
            region_config_source(None, "38"),
            RegionConfigSource::Default
        );
        assert_eq!(region_config_source(None, "19"), RegionConfigSource::Hg19);
        assert_eq!(region_config_source(None, "37"), RegionConfigSource::Hg19);
        assert_eq!(
            region_config_source(None, "chm13"),
            RegionConfigSource::Chm13
        );
        assert_eq!(
            region_config_source(Some(&PathBuf::from("/tmp/x.yaml")), "38"),
            RegionConfigSource::CustomPath
        );
    }

    #[test]
    fn should_create_vcf_dir_respects_flags_and_gene_policy() {
        let genes = vec![String::from("smn1"), String::from("CFH")];
        let no_vcf = BTreeSet::from([String::from("smn1")]);
        assert!(should_create_vcf_dir(false, &genes, &no_vcf));
        assert!(!should_create_vcf_dir(true, &genes, &no_vcf));

        let all_no_vcf = BTreeSet::from([String::from("smn1"), String::from("CFH")]);
        assert!(!should_create_vcf_dir(false, &genes, &all_no_vcf));
    }

    #[test]
    fn should_write_gene_vcf_respects_global_and_per_gene_flags() {
        let no_vcf = BTreeSet::from([String::from("smn1")]);
        assert!(should_write_gene_vcf(false, "CFH", &no_vcf));
        assert!(!should_write_gene_vcf(false, "smn1", &no_vcf));
        assert!(!should_write_gene_vcf(true, "CFH", &no_vcf));
    }

    #[test]
    fn vcf_output_dir_uses_sample_suffix() {
        let out = vcf_output_dir(Path::new("/tmp/out"), "S1");
        assert_eq!(out, PathBuf::from("/tmp/out/S1_paraphase_vcfs"));
    }

    #[test]
    fn output_paths_use_sample_suffixes() {
        let outdir = Path::new("/tmp/out");
        assert_eq!(
            output_bam_path(outdir, "S1"),
            PathBuf::from("/tmp/out/S1.paraphase.bam")
        );
        assert_eq!(
            output_json_path(outdir, "S1"),
            PathBuf::from("/tmp/out/S1.paraphase.json")
        );
    }

    #[test]
    fn prepare_output_directories_respects_novcf() -> DResult {
        let tmp = tempfile::tempdir()?;
        let genes = vec![String::from("smn1")];
        let no_vcf_genes = BTreeSet::new();
        prepare_output_directories(tmp.path(), "S1", true, &genes, &no_vcf_genes)?;
        assert!(tmp.path().exists());
        assert!(!vcf_output_dir(tmp.path(), "S1").exists());
        Ok(())
    }

    #[test]
    fn prepare_output_directories_creates_vcf_dir_when_needed() -> DResult {
        let tmp = tempfile::tempdir()?;
        let genes = vec![String::from("smn1")];
        let no_vcf_genes = BTreeSet::new();
        prepare_output_directories(tmp.path(), "S1", false, &genes, &no_vcf_genes)?;
        assert!(vcf_output_dir(tmp.path(), "S1").exists());
        Ok(())
    }

    #[test]
    fn depth_for_correction_returns_none_in_targeted_mode() {
        let depth = DepthResult {
            median: 30.0,
            median_absolute_difference: 0.1,
            sex: Sex::Other,
        };
        assert!(depth_for_correction(depth, true).is_none());
    }

    #[test]
    fn depth_for_correction_rejects_low_or_variable_depth() {
        let low = DepthResult {
            median: 9.0,
            median_absolute_difference: 0.1,
            sex: Sex::Other,
        };
        assert!(depth_for_correction(low, false).is_none());

        let variable = DepthResult {
            median: 30.0,
            median_absolute_difference: 0.3,
            sex: Sex::Other,
        };
        assert!(depth_for_correction(variable, false).is_none());
    }

    #[test]
    fn depth_for_correction_keeps_passing_depth() {
        let depth = DepthResult {
            median: 30.0,
            median_absolute_difference: 0.1,
            sex: Sex::Other,
        };
        let out = depth_for_correction(depth, false).expect("passing depth should be retained");
        assert_eq!(out.median, depth.median);
        assert_eq!(
            out.median_absolute_difference,
            depth.median_absolute_difference
        );
    }

    #[test]
    fn resolve_requested_genes_rejects_all_invalid_input() {
        let config = RegionConfig(Map::from([(String::from("smn1"), Locus::new(Map::new()))]));
        let err = resolve_requested_genes(&config, "INVALID")
            .unwrap_err()
            .to_string();
        assert!(err.contains("Please provide valid gene name(s)"));
        assert!(err.contains("smn1"));
    }

    #[test]
    fn resolve_requested_genes_filters_invalid_and_keeps_valid_order() -> DResult {
        let config = RegionConfig(Map::from([
            (String::from("CFH"), Locus::new(Map::new())),
            (String::from("smn1"), Locus::new(Map::new())),
        ]));
        let genes = resolve_requested_genes(&config, "smn1,invalid,CFH")?;
        assert_eq!(genes, vec![String::from("smn1"), String::from("CFH")]);
        Ok(())
    }

    #[test]
    fn resolve_requested_genes_defaults_to_case_insensitive_order() -> DResult {
        let config = RegionConfig(Map::from([
            (String::from("smn1"), Locus::new(Map::new())),
            (String::from("CFH"), Locus::new(Map::new())),
            (String::from("neb"), Locus::new(Map::new())),
            (String::from("ANKRD20A1"), Locus::new(Map::new())),
        ]));
        let genes = resolve_requested_genes(&config, "")?;
        assert_eq!(
            genes,
            vec![
                String::from("ANKRD20A1"),
                String::from("CFH"),
                String::from("neb"),
                String::from("smn1")
            ]
        );
        Ok(())
    }

    #[test]
    fn schedule_gene_jobs_prioritizes_configured_genes_and_preserves_others() {
        let genes = vec![
            String::from("smn1"),
            String::from("CFHR3"),
            String::from("strc"),
            String::from("FAM86B1"),
            String::from("neb"),
        ];
        let (priority, nonpriority) = schedule_gene_jobs(
            &genes,
            &[
                String::from("FAM86B1"),
                String::from("ANKRD20A1"),
                String::from("CFHR3"),
            ],
        );
        assert_eq!(
            priority,
            vec![(3, String::from("FAM86B1")), (1, String::from("CFHR3"))]
        );
        assert_eq!(
            nonpriority,
            vec![
                (0, String::from("smn1")),
                (2, String::from("strc")),
                (4, String::from("neb"))
            ]
        );
    }
}
