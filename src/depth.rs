use crate::assembly::graph_stats::compute_median;
use crate::toolkit::util::{self, DError};

use itertools::Itertools;
use rust_htslib::bam::{self, FetchDefinition::Region, Read};
use serde::{Deserialize, Serialize};

use std::path::Path;
use std::sync::Arc;

/// Reimplementation of `paraphase/genome_depth.py`.
///
/// `Settings` determines how we compute depth.
/// `Sex` is an enum of `Male/Female/Other`, important for X-linked traits.

/// Enum representing the patient sex.
#[derive(Copy, Debug, Clone, Deserialize, Serialize, PartialEq, Default)]
pub enum Sex {
    Male,
    Female,
    #[default]
    Other,
}

impl Sex {
    #[must_use]
    /// Return normalized string form used in serialized outputs.
    pub fn as_str(&self) -> &'static str {
        match self {
            Sex::Male => "male",
            Sex::Female => "female",
            Sex::Other => "unknown",
        }
    }
}

impl std::fmt::Display for Sex {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

/// `depths` vector, `x_depths`, and `y_depths`.
/// First is all depths, and the second two vectors are annotated by X and Y chromosomes.
type Tuple = (Vec<i32>, Vec<i32>, Vec<i32>);

#[derive(Copy, Clone, Debug)]
pub struct Settings {
    pub exclude_flag: u16,
    pub window_size: i64,
    pub one_based: bool,
    pub strip_chr: bool,
}

impl std::default::Default for Settings {
    fn default() -> Self {
        let exclude_flag = 0x704;
        let window_size = 1600;
        let one_based = true;
        let strip_chr = false;
        Self {
            exclude_flag,
            window_size,
            one_based,
            strip_chr,
        }
    }
}

/// Holds depth results.
/// All depths, depths for only chroms X and Y, the median depth, and the relative median absolute difference.
#[derive(Clone, Debug, Copy, Serialize, Deserialize)]
pub struct Result {
    pub median: f64,
    pub median_absolute_difference: f64,
    pub sex: Sex,
}

/// Denotes if a sample can be called based on depth of coverage.
#[derive(Copy, Debug, Clone, PartialEq, Eq)]
pub enum MedianDepthStatus {
    Passing,
    Failing,
}

const MIN_MEDIAN_COVERAGE: f64 = 10f64;
const MAX_ABS_DIFF: f64 = 0.25f64;

impl Result {
    /// Create `depth::Result` from `(all_depths, x_depths, y_depths)` tuple.
    #[must_use]
    pub fn new(depths: &Tuple) -> Self {
        let (depths, x_depths, y_depths) = depths;
        log::debug!(
            "Computing depth status from {}/{}/{}",
            depths.len(),
            x_depths.len(),
            y_depths.len()
        );
        let median = f64::from(compute_median(&depths[..]));
        let abs_diffs = depths
            .iter()
            .map(|x| (f64::from(*x) - median).abs())
            .sorted_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal))
            .collect::<Vec<f64>>();

        let median_absolute_difference = if median > 0.0 && median.is_finite() {
            compute_median(&abs_diffs[..]) as f64 / median
        } else {
            f64::NAN
        };

        let sex = Self::check_sex(&x_depths[..], &y_depths[..]);

        Self {
            median,
            median_absolute_difference,
            sex,
        }
    }

    /// Compute if the coverage is sufficient for a call.
    #[must_use]
    pub fn status(&self) -> MedianDepthStatus {
        if self.median < MIN_MEDIAN_COVERAGE || self.median_absolute_difference > MAX_ABS_DIFF {
            MedianDepthStatus::Failing
        } else {
            MedianDepthStatus::Passing
        }
    }

    /// Estimate the sex of the sample based on coverage on X and Y chromosomes.
    #[must_use]
    fn check_sex(x_depths: &[i32], y_depths: &[i32]) -> Sex {
        let median_x = compute_median(x_depths);
        let median_y = compute_median(y_depths);
        if !median_x.is_finite() || !median_y.is_finite() || median_x <= 0.0 {
            return Sex::Other;
        }
        let y_x_ratio = median_y / median_x;
        if y_x_ratio < 0.05 {
            return Sex::Female;
        } else if y_x_ratio > 0.10 {
            if median_x > 1.95 * median_y {
                return Sex::Female;
            }
            return Sex::Male;
        }
        Sex::Other
    }
}

lazy_static::lazy_static! {
    pub static ref GENOME_BACKGROUND_BYTES_38: &'static [u8] = std::include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/38/genome_region.bed"));
    pub static ref GENOME_BACKGROUND_BYTES_19: &'static [u8] = std::include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/19/genome_region.bed"));
    pub static ref GENOME_BACKGROUND_BYTES_13: &'static [u8] = std::include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/chm13/genome_region.bed"));
}

/// Computes median + mean absolute difference for depth across a set of genomic locations
/// using an indexed bam file + a set of coordinates.
///
/// Established behavior is a non-standard 1-based BED file without an end coordinate.
/// This can be altered in `Settings` using `one_based = false`.
#[derive(Debug)]
pub struct Calculator {
    pub bam: bam::IndexedReader,
    pub bed: Vec<String>,
    pub settings: Settings,
}

impl Calculator {
    /// Fallibly generate coverage calculation.
    /// # Errors
    /// 1. If bam cannot be opened.
    /// 2. If bed file is not valid utf-8.
    pub fn try_new(
        bam: impl AsRef<Path>,
        bed: &Path,
        settings: Option<Settings>,
    ) -> std::result::Result<Self, DError> {
        Self::try_from_bed_slice(bam, &std::fs::read(bed)?, settings)
    }

    /// Generate coverage calculation.
    /// # Errors
    /// 1. If bam cannot be opened.
    /// 2. If bed file cannot be read or is not valid utf-8.
    pub fn new(
        bam: impl AsRef<Path>,
        bed: &Path,
        settings: Option<Settings>,
    ) -> std::result::Result<Self, DError> {
        Self::try_new(bam, bed, settings)
    }

    /// Generate coverage calculator with bundled hg38 background regions.
    pub fn from_hg38(
        bam: impl AsRef<Path>,
        settings: Option<Settings>,
    ) -> std::result::Result<Self, DError> {
        Self::try_from_bed_slice(bam, &GENOME_BACKGROUND_BYTES_38, settings)
    }

    /// Generate coverage calculator with bundled hg19 background regions.
    pub fn from_hg19(
        bam: impl AsRef<Path>,
        settings: Option<Settings>,
    ) -> std::result::Result<Self, DError> {
        Self::try_from_bed_slice(bam, &GENOME_BACKGROUND_BYTES_19, settings)
    }

    /// Generate coverage calculator with bundled CHM13 background regions.
    pub fn from_chm13(
        bam: impl AsRef<Path>,
        settings: Option<Settings>,
    ) -> std::result::Result<Self, DError> {
        Self::try_from_bed_slice(bam, &GENOME_BACKGROUND_BYTES_13, settings)
    }

    /// Generate coverage calculator with bundled hg38 background regions using a CRAM-aware reader.
    pub fn from_hg38_with_reference(
        bam: impl AsRef<Path>,
        reference: impl AsRef<Path>,
        settings: Option<Settings>,
    ) -> std::result::Result<Self, DError> {
        Self::try_from_bed_slice_with_reference(bam, reference, &GENOME_BACKGROUND_BYTES_38, settings)
    }

    /// Generate coverage calculator with bundled hg19 background regions using a CRAM-aware reader.
    pub fn from_hg19_with_reference(
        bam: impl AsRef<Path>,
        reference: impl AsRef<Path>,
        settings: Option<Settings>,
    ) -> std::result::Result<Self, DError> {
        Self::try_from_bed_slice_with_reference(bam, reference, &GENOME_BACKGROUND_BYTES_19, settings)
    }

    /// Generate coverage calculator with bundled CHM13 background regions using a CRAM-aware reader.
    pub fn from_chm13_with_reference(
        bam: impl AsRef<Path>,
        reference: impl AsRef<Path>,
        settings: Option<Settings>,
    ) -> std::result::Result<Self, DError> {
        Self::try_from_bed_slice_with_reference(bam, reference, &GENOME_BACKGROUND_BYTES_13, settings)
    }

    /// Generate coverage calculation from a text slice.
    /// # Panics
    /// 1. If bam cannot be opened.
    /// 2. If bed file is not valid utf-8.
    pub fn try_from_bed_slice(
        bam: impl AsRef<Path>,
        bed: &[u8],
        settings: Option<Settings>,
    ) -> std::result::Result<Self, DError> {
        let bam = util::read_indexed_bam(bam.as_ref().to_string_lossy().into_owned())?;
        let settings = settings.unwrap_or_default();
        let bed = std::str::from_utf8(bed)?;
        let bed = if !settings.strip_chr {
            bed.split_terminator('\n')
                .map(std::borrow::ToOwned::to_owned)
                .collect::<Vec<_>>()
        } else {
            bed.split_terminator('\n')
                .map(std::borrow::ToOwned::to_owned)
                .map(|x| x.strip_prefix("chr").unwrap_or(&x).to_string())
                .collect::<Vec<_>>()
        };
        Ok(Self { bam, bed, settings })
    }

    /// Generate coverage calculation from a text slice using a CRAM-aware reader.
    ///
    /// # Panics
    /// 1. If bam cannot be opened.
    /// 2. If bed file is not valid utf-8.
    pub fn try_from_bed_slice_with_reference(
        bam: impl AsRef<Path>,
        reference: impl AsRef<Path>,
        bed: &[u8],
        settings: Option<Settings>,
    ) -> std::result::Result<Self, DError> {
        let bam = util::read_indexed_bam_with_reference(
            bam.as_ref().to_string_lossy().into_owned(),
            reference,
        )?;
        let settings = settings.unwrap_or_default();
        let bed = std::str::from_utf8(bed)?;
        let bed = if !settings.strip_chr {
            bed.split_terminator('\n')
                .map(std::borrow::ToOwned::to_owned)
                .collect::<Vec<_>>()
        } else {
            bed.split_terminator('\n')
                .map(std::borrow::ToOwned::to_owned)
                .map(|x| x.strip_prefix("chr").unwrap_or(&x).to_string())
                .collect::<Vec<_>>()
        };
        Ok(Self { bam, bed, settings })
    }

    /// Calculate depth and convert to `Result`.
    ///
    /// # Errors
    /// Failed to open bam file.
    ///
    /// # Panics
    /// 1. Reference names not valid UTF-8.
    /// 2. Malformatted coordinate region file.
    /// 3. Failed to parse integers in coordinate region file.
    /// 4. Any panics from `count_pos`, including:
    ///    a. If the fetched region is not present in the indexed bam.
    ///    b. If depth at position is > `i32::MAX`.
    #[must_use]
    pub fn compute(&mut self) -> Result {
        Result::new(&self.depth())
    }

    /// Perform depth computation and return raw `(all, X, Y)` depth vectors.
    ///
    /// # Errors
    /// Failed to open bam file.
    ///
    /// # Panics
    /// 1. Reference names not valid UTF-8.
    /// 2. Malformatted coordinate region file.
    /// 3. Failed to parse integers in coordinate region file.
    /// 4. Any panics from `count_pos`, including:
    ///    a. If the fetched region is not present in the indexed bam.
    #[must_use]
    pub fn depth(&mut self) -> Tuple {
        compute_depth(
            &self.bed,
            &mut self.bam,
            Some(self.settings.exclude_flag),
            Some(self.settings.window_size),
            self.settings.one_based,
        )
    }
}

/// Count reads aligned to a given position.
/// This includes deleted bases. To eliminate them, we would have to run a full pileup.
/// # Panics
/// If the fetched region is not present in the indexed bam.
/// If depth at position is > `i32::MAX`.
#[must_use]
pub fn count_pos(pos: i64, tid: i32, index: &mut bam::IndexedReader, exclude_flag: u16) -> i32 {
    if let Err(e) = index.fetch(Region(tid, pos, pos + 1)) {
        log::warn!("Failed to fetch region tid={tid}, pos={pos}: {e:?}");
        return 0;
    }
    match i32::try_from(
        index
            .rc_records()
            .flatten()
            .filter(|x| (x.flags() & exclude_flag) == 0)
            .count(),
    ) {
        Ok(count) => count,
        Err(_) => {
            log::warn!("Read count exceeded i32::MAX at tid={tid}, pos={pos}; clamping.");
            i32::MAX
        }
    }
}

/// Compute depth
/// When indexing, we add chr to the chromosomes without it and remove chr from those which have it.
///
/// This way, our computation is agnostic to the presence or absence of a chr prefix.
/// There is a risk of erroneous results if someone chooses to use conflicting names, e.g. "X" and "chrX"
/// or "chrchrX" and "chrX". "chr" would also result in an empty chromosome name.
///
/// # Errors
/// If file at `path` cannot be opened.
///
/// # Panics
/// 1. Reference names not valid UTF-8.
/// 2. Malformatted coordinate region file.
/// 3. Failed to parse integers in coordinate region file.
/// 4. Any panics from `count_pos`, including:
///    a. If the fetched region is not present in the indexed bam.
///    b. If depth at position is > `i32::MAX`.
fn compute_depth(
    lines: &[String],
    bam: &mut bam::IndexedReader,
    exclude_flag: Option<u16>,
    window_size: Option<i64>,
    one_based: bool,
) -> Tuple {
    let exclude_flag = exclude_flag.unwrap_or(0x704u16);
    let window_size = window_size.unwrap_or(1600);
    /*
    let targets = bam
        .header()
        .target_names()
        .into_iter()
        .map(VString::from)
        .collect::<Vec<_>>();
    */
    let mut target_names = bam
        .header()
        .target_names()
        .iter()
        .map(|x| (*x).to_owned())
        .enumerate()
        .map(|(id, target)| (String::from_utf8_lossy(&target).into_owned(), id as i32))
        .collect::<std::collections::BTreeMap<String, i32>>();
    log::trace!("Reference target names: {target_names:?}");
    let to_add = target_names
        .iter()
        .filter(|x| x.0.starts_with("chr"))
        .map(|(target, tid)| (target.chars().skip(3).collect::<String>(), *tid))
        .collect::<Vec<(String, i32)>>();
    // Add the names without chr.
    for (target, tid) in to_add {
        target_names.insert(target, tid);
    }
    // And add names with chr
    let to_add = target_names
        .iter()
        .filter(|x| !x.0.starts_with("chr"))
        .map(|(target, tid)| (format!("chr{target}"), *tid))
        .collect::<Vec<(String, i32)>>();
    for (target, tid) in to_add {
        target_names.insert(target, tid);
    }
    let chrom2id = Arc::new(target_names);
    let mut auto_depths = Vec::new();
    let mut x_depths = Vec::new();
    let mut y_depths = Vec::new();
    let regions_to_count = lines
        .iter()
        .filter_map(|line| {
            let mut toks = line.split_terminator('\t');
            let Some(chrom) = toks.next() else {
                log::warn!("Skipping malformed BED line (missing chromosome): {line}");
                return None;
            };
            let Some(pos_str) = toks.next() else {
                log::warn!("Skipping malformed BED line (missing position): {line}");
                return None;
            };
            let Ok(mut pos) = pos_str.parse::<i64>() else {
                log::warn!("Skipping malformed BED line (invalid position): {line}");
                return None;
            };
            pos -= i64::from(one_based);
            // Subtract 1 if one-based.
            // Default config is 1-based.
            Some((chrom.to_owned(), pos))
        })
        .collect::<Vec<_>>();
    let mut region_category = Vec::new();
    for (chrom, pos) in &regions_to_count {
        if !chrom2id.contains_key(chrom) {
            log::warn!("Chromosome {chrom} was not found in BAM targets; skipping depth calculation for this entry.")
        } else {
            let Some(tid) = chrom2id.get(chrom) else {
                log::warn!("Chromosome {chrom} was not found in BAM targets; skipping depth calculation for this entry.");
                continue;
            };
            let category = if chrom.contains('X') {
                b'X'
            } else if chrom.contains('Y') {
                b'Y'
            } else {
                b'N'
            };
            region_category.push((category, *tid, *pos));
            region_category.push((category, *tid, *pos + window_size));
        }
    }
    for (category, tid, pos) in region_category {
        let res = count_pos(pos, tid, bam, exclude_flag);
        match category {
            b'X' => {
                x_depths.push(res);
            }
            b'Y' => {
                y_depths.push(res);
            }
            b'N' => {
                auto_depths.push(res);
            }
            _ => log::warn!("Unexpected depth category token {category}; expected one of X, Y, or N. Skipping entry."),
        }
    }
    log::debug!("Autosomal depth samples: {auto_depths:?}");
    log::debug!("X-chromosome depth samples: {x_depths:?}");
    log::debug!("Y-chromosome depth samples: {y_depths:?}");
    (auto_depths, x_depths, y_depths)
}

#[cfg(test)]
mod tests {
    use super::Calculator as GenomeDepth;
    use crate::toolkit::util::test_file;

    #[test]
    fn genome_depth_ok() {
        const EPS: f64 = 1e-7;

        let bed = test_file("fake.bed");
        let result = GenomeDepth::new(
            test_file("bams/HG00733.smn1.bam").display().to_string(),
            &bed,
            None,
        )
        .expect("Failed to build genome depth calculator")
        .compute();
        assert!((19.5f64 - result.median).abs() < EPS);
        assert!((0.487_179_5f64 - result.median_absolute_difference).abs() < EPS);
    }
}
