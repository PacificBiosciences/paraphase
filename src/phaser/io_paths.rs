use crate::phaser::{Exception as PhaserException, Phaser};
use crate::toolkit::util::{DError, DResult};

use std::io::Write;
use thiserror::Error;

#[derive(Clone, Debug, Error)]
#[error("{0}")]
pub struct FaiBuildError(pub(crate) String);

/// Build `.fai` index for a FASTA path.
///
/// # Errors
/// Returns [`FaiBuildError`] when htslib indexing fails.
pub fn build_faidx(path: impl Into<std::path::PathBuf>) -> DResult {
    let path = path.into();
    let os_path = std::ffi::CString::new(path.display().to_string())?;
    let rc = unsafe { rust_htslib::htslib::fai_build(os_path.as_ptr()) };
    let msg = match rc {
        -1 => "indexing failed",
        _ => return Ok(()),
    };
    Err(
        FaiBuildError(format!("rc: {rc}. msg: {msg}. Path: {path:?}"))
            .to_string()
            .into(),
    )
}

impl Phaser {
    /// Generate path for `realigned_bam` from output directory.
    #[must_use]
    pub fn realigned_bam_path(&self) -> std::path::PathBuf {
        let suffix = format!(
            "{}_{}_realigned.bam",
            self.settings.sample_id,
            self.gene_name()
        );
        self.settings.outdir.join(suffix)
    }

    /// Write a BED file with sampled variants to the output directory.
    ///
    /// Includes heterozygous, homozygous, and non-phasing het-site buckets.
    pub fn write_variant_bed(&self) -> DResult {
        let path = self.settings.outdir.join(format!(
            "{}_{}_sampled.0based.bed",
            self.settings.sample_id,
            self.gene_name()
        ));
        let mut writer = std::io::BufWriter::new(std::fs::File::create(path)?);
        let gene_name = self.gene_name();
        let chr = self.chr().ok_or_else(|| {
            PhaserException::new("Failed to get chromosome for output BED file generation.")
        })?;
        let mut site_id = 0usize;
        writeln!(writer, "#chr\tstart\tstop\tSiteId\tRef_Var")?;
        for site in &self.het_sites {
            site_id += 1;
            writeln!(
                writer,
                "{chr}\t{}\t{}\t{gene_name}_het_site_{site_id}\t{}_{}",
                site.pos,
                site.pos + site.reference_length() as i64,
                site.ref_seq,
                site.var_seq
            )?;
        }
        for site in &self.hom_sites {
            site_id += 1;
            writeln!(
                writer,
                "{chr}\t{}\t{}\t{gene_name}_hom_site_{site_id}\t{}_{}",
                site.pos,
                site.pos + site.reference_length() as i64,
                site.ref_seq,
                site.var_seq
            )?;
        }
        for site in &self.het_sites_no_phasing {
            site_id += 1;
            writeln!(
                writer,
                "{chr}\t{}\t{}\t{gene_name}_het_site_no_phasing_{site_id}\t{}_{}",
                site.pos,
                site.pos + site.reference_length() as i64,
                site.ref_seq,
                site.var_seq
            )?;
        }
        Ok(())
    }

    /// Generate path for `realigned_tagged_path` from output directory.
    #[must_use]
    pub fn realigned_tagged_bam_path(&self) -> std::path::PathBuf {
        let suffix = format!(
            "{}_{}_realigned_tagged.bam",
            self.settings.sample_id,
            self.gene_name(),
        );
        self.settings.outdir.join(suffix)
    }

    /// Generate path for `realigned_tagged_gene2_path` from output directory.
    #[must_use]
    pub fn realigned_tagged_gene2_bam_path(&self) -> std::path::PathBuf {
        let suffix = format!(
            "{}_{}_gene2_realigned_tagged.bam",
            self.settings.sample_id,
            self.gene_name(),
        );
        self.settings.outdir.join(suffix)
    }

    /// Generate path for `realigned_gene2_path` from output directory.
    #[must_use]
    pub fn realigned_gene2_bam_path(&self) -> std::path::PathBuf {
        let suffix = format!(
            "{}_{}_gene2_realigned.bam",
            self.settings.sample_id,
            self.gene_name(),
        );
        self.settings.outdir.join(suffix)
    }

    /// Generate path for `genome_bam` from output directory.
    #[must_use]
    pub fn genome_bam_path(&self) -> std::path::PathBuf {
        self.settings.genome_bam.clone()
    }

    /// Generate path for local reference.
    #[must_use]
    pub fn local_reference_path(&self) -> std::path::PathBuf {
        let suffix = format!("{}_ref.fa", self.gene_name());
        let path = self.settings.outdir.join(suffix);
        log::debug!("Primary reference FASTA path for locus analysis: {path:?}");
        path
    }

    /// Generate path for local reference for secondary region, aka 'gene2'.
    #[must_use]
    pub fn secondary_reference_path(&self) -> std::path::PathBuf {
        let suffix = format!("{}_gene2_ref.fa", self.gene_name());
        let path = self.settings.outdir.join(suffix);
        log::debug!("Secondary (gene2) reference FASTA path for locus analysis: {path:?}");
        path
    }

    /// Generate a local-reference FASTA for the primary region and build its `.fai`.
    pub fn generate_local_reference(&self) -> Result<std::path::PathBuf, DError> {
        log::trace!("Generating primary locus reference FASTA.");
        let faidx = self.make_faidx()?;
        let (chrom, start, stop) = self.parsed_nchr_0based().ok_or_else(|| {
            PhaserException::new(format!(
                "Malformed realign region '{}'; expected 'chr:start-end' format",
                self.realign_region
            ))
        })?;
        let seq = std::str::from_utf8(&faidx.fetch_seq(chrom, start as usize, stop as usize)?)?
            .to_ascii_uppercase();
        let dest = self.local_reference_path();
        let mut output = std::io::BufWriter::new(std::fs::File::create(&dest)?);
        let realign_old = self.realign_region_old().ok_or(PhaserException::new(
            "realign_region_old failed".to_string(),
        ))?;
        writeln!(output, ">{realign_old}\n{seq}")?;
        drop(output);
        build_faidx(&dest)?;
        Ok(dest)
    }

    /// Generate a local-reference FASTA for the secondary region (`gene2`) and index it.
    pub fn generate_secondary_reference(&self) -> Result<std::path::PathBuf, DError> {
        log::trace!("Generating secondary locus reference FASTA.");
        let faidx = self.make_faidx()?;
        let (chrom, start, stop) = self.parsed_nchr_secondary_0based().ok_or_else(|| {
            let region = self
                .secondary_region_old()
                .unwrap_or_else(|| "<missing>".to_string());
            PhaserException::new(format!(
                "Malformed secondary region '{region}'; expected 'chr:start-end' format"
            ))
        })?;
        let seq = std::str::from_utf8(&faidx.fetch_seq(chrom, start as usize, stop as usize)?)?
            .to_ascii_uppercase();
        let dest = self.secondary_reference_path();
        let mut output = std::io::BufWriter::new(std::fs::File::create(&dest)?);
        writeln!(
            output,
            ">{seq_name}\n{seq}",
            seq_name = self
                .secondary_region_old()
                .ok_or(PhaserException::new("secondary_region_old failed"))?
        )?;
        drop(output);
        build_faidx(&dest)?;
        Ok(dest)
    }

    /// Return secondary reference path and whether it had to be generated.
    pub fn secondary_reference(&self) -> Result<(std::path::PathBuf, bool), DError> {
        let path = self.secondary_reference_path();
        let generated = if !path.exists() {
            self.generate_secondary_reference()?;
            true
        } else {
            false
        };
        Ok((path, generated))
    }

    /// Return local reference path and whether it had to be generated.
    pub fn local_reference(&self) -> Result<(std::path::PathBuf, bool), DError> {
        let path = self.local_reference_path();
        let generated = if !path.exists() {
            self.generate_local_reference()?;
            true
        } else {
            false
        };
        Ok((path, generated))
    }
}
