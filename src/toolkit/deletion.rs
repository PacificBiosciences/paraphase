use crate::toolkit::range;

use std::collections::BTreeSet;

/// Tunable thresholds for large-deletion detection and read classification.
#[derive(Clone, Debug)]
pub struct BigDeletionSettings {
    pub min_size: i64,
    pub min_count: i64,
    pub padding: i64,
    pub min_extend: i64, // determines how far from a clip we consider a read to be 3'/5' clipped:
    pub min_clip_len: i64, /* how many bases need to be clipped for us to consider it a clipped read. */
    pub padding_negative_reads: i64, // a read that spans this many bases around deletion ends is considered negative for the deletion.
}

impl std::default::Default for BigDeletionSettings {
    fn default() -> Self {
        Self {
            min_size: 5000,
            min_count: 3,
            padding: 50,
            min_extend: 1000,
            min_clip_len: 300,
            padding_negative_reads: 300,
        }
    }
}

/// Replaces free-standing types in Phaser for deletions.
/// One struct of repeated data, will maintain one for del1/del2
#[derive(Debug, Clone, Default)]
pub struct Datum {
    pub del_reads: BTreeSet<String>,
    pub del_reads_partial: BTreeSet<String>,
    pub del_negative_reads: BTreeSet<String>,
    pub raw: range::I64,
    pub fivep_range: range::I64,
    pub threep_range: range::I64,
}

impl Datum {
    /// Build a deletion datum from a raw deletion span and either:
    /// 1. symmetric `padding` around breakpoints, or
    /// 2. explicit 5'/3' breakpoint intervals.
    ///
    /// If neither is provided, zero-width fallback breakpoint ranges are used.
    #[must_use]
    pub fn new(
        raw: range::I64,
        padding: Option<i64>,
        del_3p_pos1: Option<i64>,
        del_3p_pos2: Option<i64>,
        del_5p_pos1: Option<i64>,
        del_5p_pos2: Option<i64>,
    ) -> Self {
        let raw = if raw.end >= raw.start {
            raw
        } else {
            log::warn!("Invalid raw deletion range {raw:?}; normalizing by swapping start/end.");
            range::I64::new(raw.end, raw.start)
        };
        if let Some(padding_value) = padding {
            Self {
                raw: raw.clone(),
                threep_range: range::I64::new(raw.start - padding_value, raw.start + padding_value),
                fivep_range: range::I64::new(raw.end - padding_value, raw.end + padding_value),
                ..Default::default()
            }
        } else if let (Some(del_5p_pos1), Some(del_5p_pos2), Some(del_3p_pos1), Some(del_3p_pos2)) =
            (del_5p_pos1, del_5p_pos2, del_3p_pos1, del_3p_pos2)
        {
            Self {
                raw,
                fivep_range: range::I64::new(del_5p_pos1, del_5p_pos2),
                threep_range: range::I64::new(del_3p_pos1, del_3p_pos2),
                ..Default::default()
            }
        } else {
            let raw_start = raw.start;
            let raw_end = raw.end;
            log::warn!(
                "Datum::new did not receive explicit breakpoint ranges or padding; using zero-width fallback breakpoint ranges."
            );
            Self {
                raw,
                fivep_range: range::I64::new(raw_end, raw_end),
                threep_range: range::I64::new(raw_start, raw_start),
                ..Default::default()
            }
        }
    }

    /// Return Python-style deletion label (`<1-based-start>_del_<len>`).
    #[must_use]
    pub fn name(&self) -> String {
        format!("{}_del_{}", self.raw.start + 1, self.raw.len())
    }

    /// Return deletion length in bases.
    #[must_use]
    pub fn size(&self) -> i64 {
        self.raw.len() as i64
    }

    /// Return the configured 5' breakpoint search interval.
    #[must_use]
    pub fn fivep(&self) -> range::I64 {
        self.fivep_range.clone()
        //range::I64::new(self.raw.start - self.padding, self.raw.start + self.padding)
    }

    /// Return the configured 3' breakpoint search interval.
    #[must_use]
    pub fn threep(&self) -> range::I64 {
        self.threep_range.clone()
        //range::I64::new(self.raw.end - self.padding, self.raw.end + self.padding)
    }

    /// Return the raw half-open deletion interval.
    #[must_use]
    pub fn range(&self) -> range::I64 {
        self.raw.clone()
    }
}
