use crate::toolkit::util::DError;

use vstr::VString;

use itertools::Itertools;

use std::collections::BTreeMap;
use std::fmt;

use super::Error;

/// Represents a read and an optional alignment id.
/// `align_id` lets us distinguish between multiple matches per read.
/// For instance, multiple repeat units in one read.
#[derive(
    Hash, PartialEq, Clone, PartialOrd, Ord, Eq, Default, serde::Deserialize, serde::Serialize,
)]
pub struct ReadAlignmentId {
    pub read_name: String,
}

impl fmt::Display for ReadAlignmentId {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.read_name)
    }
}

impl fmt::Debug for ReadAlignmentId {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{self}")
    }
}

impl ReadAlignmentId {
    /// Construct from a name, which implies that there is only one alignment for this read.
    pub fn from_name(read_name: impl Into<String>) -> Self {
        let read_name = read_name.into();
        Self {
            read_name,
            //align_id: None,
        }
    }
    /// Construct a `ReadAlignmentId` from a name and a count. This is used with multiple alignments per read.
    pub fn from_name_and_id(read_name: impl Into<String>, id: impl Into<i32>) -> Self {
        let read_name = read_name.into();
        let _align_id = Some(id.into());
        Self {
            read_name,
            //align_id,
        }
    }
    /// Return the canonical read key used in map lookups and JSON joins.
    #[must_use]
    pub fn unique_name(&self) -> String {
        self.read_name.to_string()
    }
}

impl<T: Into<String>> std::convert::From<T> for ReadAlignmentId {
    fn from(x: T) -> ReadAlignmentId {
        ReadAlignmentId::from_name(x.into())
    }
}

/// Struct for phasing site.
#[derive(Clone, Debug, Default)]
pub struct PhasingSite {
    pub pos: i64,
    pub ref_base: String,
    pub var_base: String,
    pub deletion: String,
}

impl fmt::Display for PhasingSite {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}_{}_{}", self.pos, self.ref_base, self.var_base)
    }
}

impl PhasingSite {
    /// Create `PhasingSite` from an input string.
    /// # Errors
    /// 1. If malformed. Expects "{int}_{ref}_{var}" or "{int}_del...".
    pub fn try_new(x: &str) -> Result<Self, DError> {
        if let Some((pos, ref_base, var_base)) = x.split_terminator('_').next_tuple() {
            let pos = pos.parse::<i64>()?;
            let ref_base = ref_base.to_owned();
            let var_base = var_base.to_owned();
            Ok(Self {
                pos,
                ref_base,
                var_base,
                deletion: String::new(),
            })
        } else if let Some((pos, var)) = x.split_terminator('_').next_tuple() {
            let pos = pos.parse::<i64>()?;
            if var.starts_with("del") {
                Ok(Self {
                    pos,
                    deletion: var.to_owned(),
                    ..Default::default()
                })
            } else {
                Err(Error::new(format!("var {var}")).into())
            }
        } else {
            Err(Error::new(format!("input {x}")).into())
        }
    }
}

pub type ReadFingerprintMap = BTreeMap<ReadAlignmentId, VString>;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::toolkit::util::DError;

    #[test]
    fn read_alignment_id_helpers_ok() {
        let id = ReadAlignmentId::from_name("read1");
        assert_eq!(id.read_name, "read1");
        assert_eq!(id.unique_name(), "read1");
        assert_eq!(id.to_string(), "read1");

        let id_from: ReadAlignmentId = String::from("read2").into();
        assert_eq!(id_from.read_name, "read2");

        // `from_name_and_id` currently preserves only read_name.
        let id_with_align = ReadAlignmentId::from_name_and_id("read3", 7);
        assert_eq!(id_with_align.read_name, "read3");
    }

    #[test]
    fn phasing_site_try_new_ok() -> Result<(), DError> {
        let snv = PhasingSite::try_new("100_A_T")?;
        assert_eq!(snv.pos, 100);
        assert_eq!(snv.ref_base, "A");
        assert_eq!(snv.var_base, "T");
        assert_eq!(snv.deletion, "");
        assert_eq!(snv.to_string(), "100_A_T");

        let del = PhasingSite::try_new("200_del6310")?;
        assert_eq!(del.pos, 200);
        assert_eq!(del.deletion, "del6310");
        Ok(())
    }

    #[test]
    fn phasing_site_try_new_rejects_malformed() {
        let err = PhasingSite::try_new("100_A")
            .expect_err("malformed phasing site should fail")
            .to_string();
        assert!(err.contains("var A"), "unexpected error: {err}");
    }
}
