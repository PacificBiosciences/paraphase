use crate::config::{Gene as GeneConfig, Locus as LocusConfig};
use crate::phaser::FlagBits as PhaserFlagBits;
use crate::phaser::{Exception, Phaser};
use crate::toolkit::util;

use itertools::Itertools;
use rust_htslib::{bam, bam::pileup, bam::Read, faidx};

lazy_static::lazy_static! {
    static ref NAME_REG: Option<regex::bytes::RegexSet> = regex::bytes::RegexSetBuilder::new(["ccs", "transcript", "molecule"])
        .unicode(false)
        .build()
        .ok();
}

/// Return base quality for a pileup alignment, treating deletions as max quality.
pub(crate) fn base_qual(x: &pileup::Alignment<'_>) -> u8 {
    let record = x.record();
    let qual = record.qual();
    qual.get(util::raw_qpos(x)).copied().unwrap_or(u8::MAX)
}

impl Phaser {
    /// Yields a reference to the `LocusConfig`.
    #[must_use]
    pub fn locus_config(&self) -> &LocusConfig {
        &self.config.locus
    }

    /// Yields a reference to the `GeneConfig`.
    #[must_use]
    pub fn gene_config(&self) -> &GeneConfig {
        &self.config.gene
    }

    /// Local-reference contig name (`chr_start_end`) for the primary region.
    #[must_use]
    pub fn local_chr(&self) -> Option<String> {
        self.realign_region_old()
    }

    #[must_use]
    /// Whether this locus is configured to perform phasing.
    pub fn to_phase(&self) -> bool {
        self.flag & (PhaserFlagBits::ToPhase as u8) != 0
    }

    #[must_use]
    /// Whether supplementary alignments are included in read processing.
    pub fn use_supplementary(&self) -> bool {
        self.flag & (PhaserFlagBits::UseSupplementary as u8) != 0
    }

    #[must_use]
    /// Whether the locus is evaluated on the reverse orientation.
    pub fn is_reverse(&self) -> bool {
        self.flag & (PhaserFlagBits::IsReverse as u8) != 0
    }

    #[must_use]
    /// Whether diploid copy-number (CN=2) is the expected baseline.
    pub fn expect_cn2(&self) -> bool {
        self.flag & (PhaserFlagBits::ExpectCN2 as u8) != 0
    }

    /// Gets gene name.
    #[must_use]
    pub fn gene_name(&self) -> &str {
        &self.settings.gene_name
    }

    /// Locate pivot site index inside current `het_sites`.
    ///
    /// Returns `None` when no pivot is configured or the pivot was filtered out.
    #[must_use]
    pub fn get_pivot_index(&self) -> Option<i64> {
        self.pivot_site_0based()
            .and_then(|pos| self.het_sites.iter().position(|x| x.pos == pos))
            .and_then(|pos| i64::try_from(pos).ok())
    }

    #[must_use]
    /// Sample identifier used in output naming and metadata.
    pub fn sample_id(&self) -> &str {
        &self.settings.sample_id
    }

    /// Convert `chr:pos1-pos2` format to `chr_pos1_pos2`.
    #[must_use]
    pub fn realign_region_old(&self) -> Option<String> {
        let (chr, start, stop) = self.parsed_nchr()?;
        Some(format!("{chr}_{start}_{stop}"))
    }

    /// Format secondary region in `chr_pos1_pos2` format.
    #[must_use]
    pub fn secondary_region_old(&self) -> Option<String> {
        self.locus_config()
            .gene2_region(self.settings.genome == "37")
            .and_then(|region| {
                let (chr, coords) = region.split_terminator(':').next_tuple()?;
                let (start, stop) = coords.split_terminator('-').next_tuple()?;
                Some(format!("{chr}_{start}_{stop}"))
            })
    }

    /// Parse `realign_region` into `(chr, start, stop)` using one-based coordinates.
    #[must_use]
    pub fn parsed_nchr(&self) -> Option<(&str, &str, &str)> {
        let (chr, coords) = self.realign_region.split_terminator(':').next_tuple()?;
        let (start, stop) = coords.split_terminator('-').next_tuple()?;
        Some((chr, start, stop))
    }

    /// Parse `realign_region` into `(chr, start, stop)` using zero-based coordinates.
    #[must_use]
    pub fn parsed_nchr_0based(&self) -> Option<(&str, i64, i64)> {
        let (chr, coords) = self.realign_region.split_terminator(':').next_tuple()?;
        let (start, stop) = coords
            .split_terminator('-')
            .filter_map(|x| x.parse::<i64>().ok())
            .next_tuple()?;
        Some((chr, start - 1, stop - 1))
    }

    /// Parse the configured secondary region into `(chr, start, stop)` with zero-based coordinates.
    #[must_use]
    pub fn parsed_nchr_secondary_0based(&self) -> Option<(&str, i64, i64)> {
        self.locus_config()
            .gene2_region(self.settings.genome == "37")
            .and_then(|region| {
                let (chr, coords) = region.split_terminator(':').next_tuple()?;
                let (start, stop) = coords
                    .split_terminator('-')
                    .filter_map(|x| x.parse::<i64>().ok())
                    .next_tuple()?;
                Some((chr, start - 1, stop - 1))
            })
    }

    /// Chromosome name parsed from the primary realign region string.
    #[must_use]
    pub fn chr(&self) -> Option<&str> {
        self.parsed_nchr().map(|x| x.0)
    }

    #[must_use]
    /// Resolve the chromosome TID in the genome BAM header for this locus.
    pub fn genome_tid(&self) -> Option<u32> {
        let bam = self.try_genome_bam().ok()?;
        self.chr().and_then(|chr| bam.header().tid(chr.as_bytes()))
    }

    /// Build FASTA reader from the configured genome reference path.
    pub fn make_faidx(&self) -> Result<faidx::Reader, rust_htslib::errors::Error> {
        faidx::Reader::from_path(&self.settings.genome_reference)
    }

    /// Build FASTA reader from the local reference path.
    pub fn make_local_faidx(&self) -> Result<faidx::Reader, Exception> {
        let res = faidx::Reader::from_path(
            self.local_reference()
                .map_err(|e| {
                    Exception::new(format!(
                        "Error in make_local_faidx local reference path: {e:?}. Path: {:?}",
                        self.local_reference()
                    ))
                })?
                .0,
        )
        .map_err(|e| {
            Exception::new(format!(
                "Error in make_local_faidx from_path: {e:?}. Path: {:?}",
                self.local_reference()
            ))
        });
        res
    }

    /// Build FASTA reader from the secondary (gene2) local reference path.
    pub fn make_local_faidx_gene2(&self) -> Result<faidx::Reader, Exception> {
        let res = faidx::Reader::from_path(
            self.secondary_reference()
                .map_err(|e| {
                    Exception::new(format!(
                        "Error in make_local_faidx_gene2 local reference path: {e:?}. Path: {:?}",
                        self.secondary_reference()
                    ))
                })?
                .0,
        )
        .map_err(|e| {
            Exception::new(format!(
                "Error in make_local_faidx_gene2 from_path: {e:?}. Path: {:?}",
                self.secondary_reference()
            ))
        });
        res
    }

    /// Provides zero-based start offset for the realign region.
    pub fn try_offset(&self) -> Result<i64, Exception> {
        let (_chr, start, _stop) = self
            .parsed_nchr()
            .ok_or_else(|| Exception::new("Failed to parse region string"))?;
        let start = start.parse::<i64>().map_err(|e| {
            Exception::new(format!("Failed to parse start coordinate '{start}': {e}"))
        })?;
        Ok(start - 1)
    }

    /// Provides zero-based start offset for the realign region.
    #[must_use]
    pub fn offset(&self) -> i64 {
        self.try_offset().unwrap_or_else(|e| {
            log::warn!("Failed to compute primary-region offset: {e}. Using fallback offset 0.");
            0
        })
    }

    /// Provides zero-based start offset for the secondary region.
    #[must_use]
    pub fn secondary_offset(&self) -> Option<i64> {
        self.parsed_nchr_secondary_0based().map(|x| x.1)
    }

    /// Right boundary for gene in one-based coordinates.
    pub fn try_right_boundary(&self) -> Result<i64, Exception> {
        if let Some(right) = self.right_boundary {
            Ok(right)
        } else {
            let (_chr, _start, stop) = self
                .parsed_nchr()
                .ok_or_else(|| Exception::new("Failed to parse region string"))?;
            let stop = stop.parse::<i64>().map_err(|e| {
                Exception::new(format!("Failed to parse right boundary '{stop}': {e}"))
            })?;
            Ok(stop)
        }
    }

    /// Right boundary for gene in one-based coordinates.
    #[must_use]
    pub fn right_boundary(&self) -> i64 {
        self.try_right_boundary().unwrap_or_else(|e| {
            log::warn!(
                "Failed to compute right boundary from locus coordinates: {e}. Using fallback value 0."
            );
            0
        })
    }

    /// Right boundary for gene in zero-based coordinates.
    #[must_use]
    pub fn right_boundary_0based(&self) -> i64 {
        self.right_boundary() - 1
    }

    /// Left boundary for gene in one-based coordinates.
    pub fn try_left_boundary(&self) -> Result<i64, Exception> {
        if let Some(left) = self.left_boundary {
            Ok(left)
        } else {
            let (_chr, start, _stop) = self
                .parsed_nchr()
                .ok_or_else(|| Exception::new("Failed to parse region string"))?;
            let start = start.parse::<i64>().map_err(|e| {
                Exception::new(format!("Failed to parse left boundary '{start}': {e}"))
            })?;
            Ok(start)
        }
    }

    /// Left boundary for gene in one-based coordinates.
    #[must_use]
    pub fn left_boundary(&self) -> i64 {
        self.try_left_boundary().unwrap_or_else(|e| {
            log::warn!("Failed to compute left boundary from locus coordinates: {e}. Using fallback value 0.");
            0
        })
    }

    /// Left boundary for gene in zero-based coordinates.
    #[must_use]
    pub fn left_boundary_0based(&self) -> i64 {
        self.left_boundary() - 1
    }

    #[must_use]
    /// Pivot position converted from one-based to zero-based coordinates.
    pub fn pivot_site_0based(&self) -> Option<i64> {
        self.pivot_site.map(|x| x - 1)
    }

    /// Coordinate for `gene_start` (one-based).
    #[must_use]
    pub fn gene_start(&self) -> i64 {
        self.gene_start.unwrap_or_else(|| self.left_boundary())
    }

    /// Coordinate for `gene_end` (one-based).
    #[must_use]
    pub fn gene_end(&self) -> i64 {
        self.gene_end.unwrap_or_else(|| self.right_boundary())
    }

    /// Open a `bam::IndexedReader` from `self.realigned_bam_path()`.
    pub fn try_realigned_bam(&self) -> Result<bam::IndexedReader, Exception> {
        util::read_indexed_bam(self.realigned_bam_path().display().to_string())
            .map_err(|e| Exception::new(format!("Failed to open realigned bam file: {e}")))
    }

    /// Open a `bam::IndexedReader` from `self.realigned_bam_path()`.
    #[must_use]
    pub fn realigned_bam(&self) -> Option<bam::IndexedReader> {
        self.try_realigned_bam()
            .map_err(|e| {
                log::warn!("Failed to open aligned-locus BAM from generated outputs: {e}");
                e
            })
            .ok()
    }

    /// Open a `bam::IndexedReader` from `self.genome_bam_path()`.
    pub fn try_genome_bam(&self) -> Result<bam::IndexedReader, Exception> {
        util::read_indexed_bam_with_reference(
            self.genome_bam_path().display().to_string(),
            &self.settings.genome_reference,
        )
            .map_err(|e| Exception::new(format!("Failed to open genome alignment file: {e}")))
    }

    /// Open a `bam::IndexedReader` from `self.genome_bam_path()`.
    #[must_use]
    pub fn genome_bam(&self) -> Option<bam::IndexedReader> {
        self.try_genome_bam()
            .map_err(|e| {
                log::warn!("Failed to open genome BAM input: {e}");
                e
            })
            .ok()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config;
    use crate::phaser;
    use crate::toolkit::site_selection::CandidateSite;
    use crate::toolkit::util::{self, DResult};
    use std::str::FromStr;

    #[test]
    fn get_pivot_index_ok() -> DResult {
        let outdir = tempfile::TempDir::new()?;
        let genome_bam = util::test_file("HG00733_smn1_realigned.bam");
        let genome_path = if let Ok(x) = std::env::var("HG38") {
            x.trim_end_matches(".mmi").to_string()
        } else {
            log::warn!(
                "Skipping core_utils test because the HG38 environment variable is not configured."
            );
            return Ok(());
        };
        let settings = phaser::Settings::new(
            "HG00733",
            (genome_path, genome_bam),
            outdir.path(),
            "smn1",
            &config::Region::try_load(None)?,
            None,
            None,
            String::from("38"),
            None,
            0.03,
            false,
        );

        let gene_config = config::Gene::try_load(None)?;
        let mut phaser = Phaser::new(settings, Some(gene_config), None, None)?;

        phaser.het_sites = vec![
            CandidateSite::from_str("70951940_A_C")?,
            CandidateSite::from_str("70951946_T_G")?,
        ];
        assert_eq!(phaser.get_pivot_index(), Some(1));

        phaser.het_sites = vec![
            CandidateSite::from_str("70951940_A_C")?,
            CandidateSite::from_str("70951946_T_G")?,
            CandidateSite::from_str("70951958_T_G")?,
        ];
        assert_eq!(phaser.get_pivot_index(), Some(1));

        phaser.het_sites = vec![
            CandidateSite::from_str("70951947_A_C")?,
            CandidateSite::from_str("70951949_T_G")?,
            CandidateSite::from_str("70951958_T_G")?,
        ];
        assert_eq!(phaser.get_pivot_index(), None);
        Ok(())
    }
}
