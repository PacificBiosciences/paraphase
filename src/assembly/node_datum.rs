use std::fmt;
use vstr::VString;

/// Core datum for the graph representation.
/// Has the position, the sequence, and the bitflag associated.
/// Curently, it only includes `is_del`, but we could extend it in the future.
#[derive(Clone, Default, Hash, Eq, PartialEq, PartialOrd, Ord)]
pub struct NodeDatum {
    pub pos: u32,
    pub hap: VString,
    pub flag: u32,
}

impl std::fmt::Debug for NodeDatum {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "NodeDatum{{{}@{}|flag:{}}}",
            self.hap, self.pos, self.flag
        )
    }
}

#[repr(u32)]
/// Bit flags stored in `NodeDatum.flag`.
pub enum NodeFlag {
    IsDel = 1,
}

#[repr(u32)]
/// Bit offsets for `NodeFlag` values.
pub enum NodeFlagShift {
    IsDel = 0,
}

impl std::ops::Deref for NodeDatum {
    type Target = VString;
    fn deref(&self) -> &Self::Target {
        &self.hap
    }
}

impl NodeDatum {
    /// Create a NodeDatum, with a sequence and a position within site space.
    pub fn new(hap: impl Into<VString>, pos: u32) -> Self {
        let hap = hap.into();
        Self { hap, pos, flag: 0 }
    }
    /// True if the datum has been deleted from the graph.
    #[must_use]
    pub fn is_del(&self) -> bool {
        (self.flag & (NodeFlag::IsDel as u32)) != 0
    }
    /// Marks node as deleted.
    pub fn set_is_del(&mut self) {
        self.flag |= NodeFlag::IsDel as u32;
    }
}

static_assertions::const_assert_eq!(std::mem::size_of::<NodeDatum>(), 32);

#[derive(Debug)]
pub struct NodeDatumParserError {
    msg: String,
}

impl std::error::Error for NodeDatumParserError {}

impl std::fmt::Display for NodeDatumParserError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "NodeDatumParserError{{Failed to parse NodeDatum from str. message: \"{}\"}}",
            self.msg
        )
    }
}

impl std::convert::TryFrom<&str> for NodeDatum {
    type Error = NodeDatumParserError;
    fn try_from(x: &str) -> Result<Self, Self::Error> {
        use itertools::Itertools;
        if let Some((hap, pos)) = x.split_terminator('-').next_tuple() {
            match pos.parse::<u32>() {
                Ok(pos) => Ok(NodeDatum::new(hap, pos)),
                Err(e) => Err(NodeDatumParserError {
                    msg: format!("Failed to parse position from NodeDatum string {x}. Error {e:?}"),
                }),
            }
        } else {
            Err(NodeDatumParserError {
                msg: format!("Error: could not decompose string into NodeDatum. String: {x}"),
            })
        }
    }
}

impl<T> std::convert::From<(u8, T)> for NodeDatum
where
    T: std::convert::TryInto<u32> + std::fmt::Debug + std::marker::Copy,
{
    /// Build a single-base node datum from `(base, pos)` tuple.
    fn from(x: (u8, T)) -> Self {
        let (base, pos) = x;
        let hap = VString::from(vec![base]);
        let pos = pos.try_into().unwrap_or_else(|_| {
            log::warn!(
                "Failed to convert node position to u32: {:?}. Falling back to position 0.",
                x.1
            );
            0
        });
        Self { hap, pos, flag: 0 }
    }
}
