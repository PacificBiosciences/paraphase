/// Lightweight wrapper around `std::ops::Range` with serde support and
/// ordering/display helpers used across paraphase coordinate logic.
#[derive(PartialEq, Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Range<T: std::cmp::Eq + std::cmp::Ord + std::marker::Copy + std::fmt::Display> {
    pub inner: std::ops::Range<T>,
}

impl<
        T: std::cmp::Eq
            + std::cmp::Ord
            + std::marker::Copy
            + std::fmt::Display
            + std::default::Default,
    > std::default::Default for Range<T>
{
    fn default() -> Self {
        Self::new(T::default(), T::default())
    }
}

impl<T: std::cmp::Eq + std::cmp::Ord + std::marker::Copy + std::fmt::Display> std::ops::DerefMut
    for Range<T>
{
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.inner
    }
}

impl<
        T: std::cmp::Eq
            + std::cmp::Ord
            + std::marker::Copy
            + std::fmt::Display
            + std::fmt::Debug
            + std::ops::Sub<Output = T>
            + std::convert::TryInto<usize>,
    > Range<T>
{
    ///
    /// Length of range when used half-closed (bam/bed convention).
    /// Ideally, we would use `std::ops::Range` `len()`, but `std::iter::ExactSizeIterator`, where this is implemented,
    /// is not available for `i64`.
    ///
    /// `is_empty` is provided by `std::ops::Range`.
    ///
    ///```
    /// let r = paraphase::toolkit::range::Range::<i64>::new(0, 10);
    /// assert!(!r.is_empty());
    /// assert_eq!(r.len(), 10);
    /// let r = paraphase::toolkit::range::Range::<i64>::new(0, 0);
    /// assert!(r.is_empty());
    /// assert_eq!(r.len(), 0);
    /// ```
    ///
    #[must_use]
    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize {
        (self.end - self.start).try_into().unwrap_or_else(|_| {
            log::warn!(
                "Range::len conversion to usize failed for range {}; returning fallback value 0.",
                self
            );
            0
        })
    }
}

impl<T: std::cmp::Eq + std::cmp::Ord + std::marker::Copy + std::fmt::Display> std::ops::Deref
    for Range<T>
{
    type Target = std::ops::Range<T>;
    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl<T: std::cmp::Ord + std::marker::Copy + std::fmt::Display> std::cmp::Eq for Range<T> {}

impl<T: std::cmp::Eq + std::cmp::Ord + std::marker::Copy + std::fmt::Display> std::cmp::Ord
    for Range<T>
{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (self.inner.start, self.inner.end).cmp(&(other.inner.start, other.inner.end))
    }
}

impl<T: std::cmp::Eq + std::cmp::Ord + std::marker::Copy + std::fmt::Display> std::cmp::PartialOrd
    for Range<T>
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl<T: std::cmp::Eq + std::cmp::Ord + std::marker::Copy + std::fmt::Display> Range<T> {
    /// Construct a half-open interval `[start, end)`.
    pub fn new(start: T, end: T) -> Self {
        let inner = std::ops::Range::<T> { start, end };
        Self { inner }
    }
}

impl<T: std::cmp::Eq + std::cmp::Ord + std::marker::Copy + std::fmt::Display> std::fmt::Display
    for Range<T>
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "[{}-{})", self.start, self.end)
    }
}

impl<T: std::cmp::Eq + std::cmp::Ord + std::marker::Copy + std::fmt::Display>
    std::convert::From<(T, T)> for Range<T>
{
    fn from(x: (T, T)) -> Range<T> {
        Range::new(x.0, x.1)
    }
}

impl<T: std::cmp::Eq + std::cmp::Ord + std::marker::Copy + std::fmt::Display>
    std::convert::From<&(T, T)> for Range<T>
{
    fn from(x: &(T, T)) -> Range<T> {
        Range::new(x.0, x.1)
    }
}

pub type I64 = Range<i64>;
pub type I32 = Range<i32>;
