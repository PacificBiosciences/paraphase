use clap::Parser;
use paraphase::cli::{setup_logging, Settings};
use paraphase::config::gene::Config as GeneConfig;
use paraphase::config::region::Config as RegionConfig;
use paraphase::config::Gene;
use paraphase::toolkit::pipeline::*;
use paraphase::toolkit::util::{DResult, FULL_VERSION_PROGRAM, GIT_DESCRIBE};

fn main() -> DResult {
    let args = Settings::parse().check();
    run(args)
}

/// Execute one full Paraphase run from validated CLI settings.
///
/// # Errors
/// Returns an error if configuration loading, pipeline execution, or output writing fails.
pub fn run(args: Settings) -> DResult {
    setup_logging(&args);

    #[cfg(feature = "pprof")]
    let guard = args
        .profile_freq
        .and_then(|freq| match pprof::ProfilerGuard::new(freq) {
            Ok(guard) => Some(guard),
            Err(e) => {
                log::warn!("Profiling disabled: failed to start profiler at frequency {freq}: {e}");
                None
            }
        });

    log::info!("Running: {} ({})", &**FULL_VERSION_PROGRAM, &**GIT_DESCRIBE);
    log_cli_invocation();
    let gene_config: GeneConfig = Gene::try_load(None)?;
    let config: RegionConfig = load_region_config(&args)?;

    let genes = resolve_requested_genes(&config, &args.gene)?;
    log::debug!("Genes: {genes:?} from args {:?}", args.gene);

    let depth = compute_depth_result(&args)?;

    log::debug!("Genome depth: {depth:?}");
    let sample_sex = depth.sex;
    let genome_depth = depth_for_correction(depth, args.targeted);

    let sample = resolve_sample_name(&args);
    run_pipeline_for_sample(
        &args,
        &config,
        &gene_config,
        &genes,
        &sample,
        genome_depth,
        sample_sex,
    )?;

    // Profile
    #[cfg(feature = "pprof")]
    write_profile_flamegraph(guard, args.outdir.as_path())?;

    log::info!("Completed Paraphase analysis.");
    Ok(())
}
