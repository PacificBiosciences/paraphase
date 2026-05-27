use crate::assembly::graph_stats::{compute_median, percentile_i32};
use crate::depth::count_pos;
use crate::phaser::Phaser;
use crate::toolkit::range;

use rust_htslib::{bam, bam::Read};
use vstr::{VStr, VString};

pub(crate) const NUM_SAMPLED_DEPTH_POSITIONS: usize = 100;
pub(crate) const PERCENTILE: i32 = 80;

impl Phaser {
    /// Computes median depth for a given region using a given stride for sampling.
    #[must_use]
    pub fn regional_depth(
        bam: &mut bam::IndexedReader,
        tid: i32,
        query: &[range::I64],
        num_intervals: Option<usize>,
        exclude_flag: Option<u16>,
        one_based: Option<bool>,
        percentile: Option<i32>,
    ) -> Vec<(f32, f32)> {
        log::trace!(
            "Depth query target context: tid={tid}, target_name={}, all_targets={:?}",
            VStr::from(bam.header().target_names()[tid as usize]),
            bam.header()
                .target_names()
                .into_iter()
                .map(VString::from)
                .collect::<Vec<_>>(),
        );
        let one_based = usize::from(one_based.unwrap_or(true));
        let exclude_flag = exclude_flag.unwrap_or(0x704u16);
        let percentile = percentile.unwrap_or(PERCENTILE);
        log::debug!(
            "Computing regional depth: query_count={}, num_intervals={num_intervals:?}",
            query.len()
        );
        let num_intervals = num_intervals.unwrap_or(NUM_SAMPLED_DEPTH_POSITIONS);
        let depths = query
            .iter()
            .map(|q| {
                let step_by = std::cmp::max(1usize, q.len() / num_intervals);
                log::debug!("Regional depth sampling stride: step_by={step_by}, requested_intervals={num_intervals}");
                let sampled_positions = (q.start..q.end)
                    .step_by(std::cmp::max(1usize, q.len() / num_intervals))
                    .collect::<Vec<_>>();
                log::debug!("Sampling regional depth at positions: {sampled_positions:?}");
                let depths = (q.start..q.end)
                    .step_by(std::cmp::max(1usize, q.len() / num_intervals))
                    .map(|pos| count_pos(pos - one_based as i64, tid, bam, exclude_flag))
                    .collect::<Vec<_>>();
                log::debug!("Sampled depth values: {depths:?}");
                let median = compute_median(&depths);
                let percentile = percentile_i32(&depths, percentile) as f32;
                (median, percentile)
            })
            .collect::<Vec<_>>();
        depths
    }

    /// Determines if the sample should be failed for this region for coverage reasons.
    #[must_use]
    pub fn coverage_passes(&mut self) -> bool {
        let Ok(mut bam) = self.try_realigned_bam() else {
            log::warn!(
                "coverage_passes: failed to open realigned bam; treating as failed coverage."
            );
            return false;
        };
        let Some(genome_tid) = self.genome_tid().map(|x| x as i32) else {
            log::warn!(
                "coverage_passes: missing genome tid in BAM header; treating as failed coverage."
            );
            return false;
        };
        let left = self.left_boundary_0based();
        let right = self.right_boundary_0based();
        // Match Python window-selection logic exactly:
        // 1) Start with full interval.
        // 2) If both clip lists exist, only use clipped interval when max(5p) < min(3p).
        // 3) Else if only one side exists, trim that side only.
        let mut check_region_start = left;
        let mut check_region_end = right;
        if !self.clip_5p_positions.is_empty() && !self.clip_3p_positions.is_empty() {
            let clip_5p_max = *self.clip_5p_positions.iter().max().unwrap_or(&left);
            let clip_3p_min = *self.clip_3p_positions.iter().min().unwrap_or(&right);
            if clip_5p_max < clip_3p_min {
                check_region_start = clip_5p_max;
                check_region_end = clip_3p_min;
            }
        } else if !self.clip_5p_positions.is_empty() {
            let clip_5p_max = *self.clip_5p_positions.iter().max().unwrap_or(&left);
            check_region_start = std::cmp::max(left, clip_5p_max);
        } else if !self.clip_3p_positions.is_empty() {
            let clip_3p_min = *self.clip_3p_positions.iter().min().unwrap_or(&right);
            check_region_end = std::cmp::min(right, clip_3p_min);
        }
        let depth1 = Self::regional_depth(
            &mut bam,
            genome_tid,
            &[range::I64::new(left, right)],
            None,
            Some(0),
            Some(false),
            Some(PERCENTILE),
        );
        self.region_avg_depth = depth1.clone();
        log::debug!(
            "Coverage window bounds after clip adjustment: left={}, right={}",
            check_region_start,
            check_region_end
        );
        if check_region_start != left || check_region_end != right {
            let depth2 = Self::regional_depth(
                &mut bam,
                genome_tid,
                &[range::I64::new(check_region_start, check_region_end)],
                None,
                Some(0),
                Some(false),
                Some(PERCENTILE),
            );
            log::debug!(
                "Coverage comparison full-vs-trimmed: full={:?}, trimmed={:?}",
                depth1[0],
                depth2[0]
            );
            if depth1[0].0.is_nan() || depth2[0].0 > depth1[0].0 {
                self.region_avg_depth = depth2;
            }
        }
        let (median, percentile) = self.region_avg_depth[0];
        log::debug!(
            "Coverage stats after depth sampling: median={median}, percentile_{PERCENTILE}={percentile}; thresholds: median>8 OR percentile>=50"
        );
        (median > 8. || percentile >= 50.) || self.settings.allow_low_coverage
    }
}
