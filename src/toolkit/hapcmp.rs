use crate::toolkit::util::DError;
use anyhow::anyhow;
use vstr::VStr;

/// Result from comparing two haplotypes/fingerprint strings.
/// matches, mismatches, and extensions are counted.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd, Default)]
pub struct HapCompare {
    pub matches: i32,
    pub mismatches: i32,
    pub extra: i32,
}

/// Identify if a site label is a gap or a flank.
#[must_use]
pub fn is_gap_or_flank(x: impl Into<u8>) -> bool {
    let x = x.into();
    x == b'x'
}

impl HapCompare {
    /// Create a `HapCompare` object.
    #[must_use]
    pub fn new(matches: i32, mismatches: i32, extra: i32) -> Self {
        Self {
            matches,
            mismatches,
            extra,
        }
    }

    /// Core comparison function. Works on any type that can create a `VStr`, which includes
    /// &str, String, and `VStr`.
    /// To use with `VString`, call the .`vstr()` function.
    ///
    /// ```
    /// use paraphase::toolkit::hapcmp::HapCompare;
    /// use vstr::VStr;
    /// let res = HapCompare::from_haps(VStr::from("12x"), VStr::from("x22")).unwrap();
    /// assert_eq!(res.matches, 1);
    /// assert_eq!(res.mismatches, 0);
    /// assert_eq!(res.extra, 2);
    /// let res = HapCompare::from_haps(VStr::from("11x"), VStr::from("x22")).unwrap();
    /// assert_eq!(res.matches, 0);
    /// assert_eq!(res.mismatches, 1);
    /// assert_eq!(res.extra, 2);
    /// let res = HapCompare::from_haps(VStr::from("12xx"), VStr::from("xx22")).unwrap();
    /// assert_eq!(res.matches, 0);
    /// assert_eq!(res.mismatches, 0);
    /// assert_eq!(res.extra, 4);
    /// let res = HapCompare::from_haps(VStr::from("1234"), VStr::from("1234")).unwrap();
    /// assert_eq!(res.matches, 4);
    /// assert_eq!(res.mismatches, 0);
    /// assert_eq!(res.extra, 0);
    /// let res = HapCompare::from_haps(VStr::from("1134"), VStr::from("1234")).unwrap();
    /// assert_eq!(res.matches, 3);
    /// assert_eq!(res.mismatches, 1);
    /// assert_eq!(res.extra, 0);
    /// ```
    pub fn from_haps<'a, 'b>(
        lhs: impl Into<VStr<'a>>,
        rhs: impl Into<VStr<'b>>,
    ) -> Result<Self, DError> {
        let lhs = lhs.into();
        let rhs = rhs.into();
        if lhs.len() != rhs.len() {
            Err(anyhow!("Lengths don't match: {lhs} and {rhs}"))?;
        }
        Ok(lhs.iter().copied().zip(rhs.iter().copied()).fold(
            Self::default(),
            |mut acc, (lhb, rhb)| {
                match i32::from(is_gap_or_flank(lhb)) | i32::from(is_gap_or_flank(rhb)) {
                    0 => {
                        if lhb == rhb {
                            acc.matches += 1;
                        } else {
                            acc.mismatches += 1;
                        }
                    }
                    1 => {
                        acc.extra += 1;
                    }
                    _ => {}
                }
                acc
            },
        ))
    }
}
