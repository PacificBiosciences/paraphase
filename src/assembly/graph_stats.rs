use itertools::Itertools;

/// Compute the median of a slice of integers.
/// Assigns a float, averages midpoints in case of an even-sized slice.
/// ```
/// use paraphase::assembly::graph_stats::compute_median;
/// let x = &[0, 3i32, 4];
/// assert_eq!(3., compute_median(x));
/// assert!(compute_median(&[] as &[i32]).is_nan());
/// let x = &[4i32, 2];
/// assert_eq!(compute_median(x), 3.);
/// let x = &[1i32, 2];
/// assert_eq!(compute_median(x), 1.5f32);
/// let x = &[1.0f32, 5.0, 2.0];
/// assert_eq!(compute_median(x), 2.0f32);
/// ```
pub trait ToMedianF32 {
    /// Convert numeric-like values to `f32` for median computation.
    fn to_median_f32(self) -> f32;
}

impl ToMedianF32 for f32 {
    fn to_median_f32(self) -> f32 {
        self
    }
}

impl ToMedianF32 for f64 {
    fn to_median_f32(self) -> f32 {
        self as f32
    }
}

impl ToMedianF32 for i32 {
    fn to_median_f32(self) -> f32 {
        self as f32
    }
}

impl ToMedianF32 for i64 {
    fn to_median_f32(self) -> f32 {
        self as f32
    }
}

impl ToMedianF32 for u32 {
    fn to_median_f32(self) -> f32 {
        self as f32
    }
}

impl ToMedianF32 for usize {
    fn to_median_f32(self) -> f32 {
        self as f32
    }
}

#[must_use]
pub fn compute_median<T>(x: &[T]) -> f32
where
    T: Copy + ToMedianF32 + PartialOrd,
{
    if x.is_empty() {
        return f32::NAN;
    }
    let copied = x
        .iter()
        .copied()
        .map(ToMedianF32::to_median_f32)
        .sorted_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .collect::<Vec<_>>();
    let len = copied.len();
    let midpoint = len / 2;
    if len & 1 != 0 {
        copied[midpoint]
    } else {
        (copied[midpoint] + copied[midpoint - 1]) * 0.5
    }
}

/// Compute percentile over an array.
/// Matches np.percentile with method = "linear"
///```
/// use paraphase::assembly::graph_stats::percentile_i32;
/// let a: &[i32] = &[10, 7, 4, 3, 2, 1];
/// assert!((percentile_i32(a, 50) - 3.5).abs() < 1e-10);
/// let a: &[i32] = &[
///     41, 20, 17, 48, 11, 16, 2, 21, 26, 40, 12, 7, 48, 14, 28, 44, 37, 23, 15, 41, 24, 18, 49,
///     39, 35, 18, 0, 22, 19, 34, 39, 29, 4, 10, 4, 11, 10, 30, 26, 33, 14, 17, 35, 9, 32, 47, 28,
///     11, 19, 30, 11, 35, 39, 42, 39, 37, 26, 47, 49, 1, 31, 35, 22, 11, 3, 25, 30, 9, 44, 24,
///     19, 14, 9, 35, 27, 36, 16, 10, 20, 44, 14, 6, 42, 13, 23, 36, 3, 43, 4, 8, 13, 31, 25, 10,
///     41, 20, 28, 47, 47, 9,
/// ];
/// assert!(
///     (percentile_i32(a, 80) - 39.0).abs() < 1e-10,
///     "{} vs expected 39.0",
///     percentile_i32(a, 80)
/// );
/// assert!(
///     (percentile_i32(a, 60) - 28.4).abs() < 1e-10,
///     "{} vs expected 28.4",
///     percentile_i32(a, 60)
/// );
/// ```
#[must_use]
pub fn percentile_i32(x: &[i32], percentile: i32) -> f64 {
    if x.is_empty() {
        return f64::NAN;
    }
    let copied = x.iter().copied().sorted().collect::<Vec<_>>();
    let len = copied.len();
    const ALPHA: i32 = 1;
    const BETA: i32 = 1;
    let percentile = f64::from(percentile) * 0.01;
    let n_mul = f64::from(len as i32 - ALPHA - BETA + 1);

    let virtual_index = percentile * n_mul + f64::from(ALPHA);
    let accessed_index = (virtual_index as usize) - 1;
    let fract = virtual_index.fract();
    let accessed_val = f64::from(copied[accessed_index]);
    if fract != 0.0 {
        accessed_val * (1. - fract) + f64::from(copied[accessed_index + 1]) * fract
    } else {
        accessed_val
    }
}

/// Find last index that is not 'x', or 0 if no such base exists.
/// ```
/// paraphase::toolkit::util::init_log(log::LevelFilter::Debug);
/// use paraphase::assembly::graph_stats::count_suffix_x;
/// assert_eq!(count_suffix_x(&[0, 1, b'x']), 1, "0, 1, x");
/// assert_eq!(count_suffix_x(&[0, 1]), 1, "0, 1");
/// assert_eq!(count_suffix_x(&[0, 1, b'x', 1]), 3, "0x1x");
/// assert_eq!(count_suffix_x(b"xxx122"), 5, "xxx122");
/// assert_eq!(count_suffix_x(b"xxxxxx"), 0, "xxxxxx");
/// assert_eq!(count_suffix_x(b"1xxxxx"), 0, "1xxxxx");
/// assert_eq!(count_suffix_x(b"1xxxxx1"), 6, "1xxxxx1");
/// ```
#[must_use]
pub fn count_suffix_x(x: &[u8]) -> usize {
    x.iter().rposition(|item| *item != b'x').unwrap_or(0)
}
