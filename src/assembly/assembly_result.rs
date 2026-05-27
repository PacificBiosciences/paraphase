use serde::{Deserialize, Serialize};
use vstr::VString;

use std::collections::BTreeSet;
use std::fmt;
use std::sync::OnceLock;

/// Results from running `variant_graph::Graph`
/// Has all possible haplotypes, the main set selected, and the highest copy number supported by the data.
/// `main_haps` corresponds to `ass_haps` in `Paraphase.paraphase.phaser` and `main_haps` in `Paraphase.phaser.haplotype_assembler`.
/// `final_haps` corresponds to `original_haps` in `Paraphase.paraphase.phaser` and `final_haps` in `Paraphase.phaser.haplotype_assembler`.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct AssemblyResult {
    pub main_haps: AssembledPaths,
    pub final_haps: AssembledPaths,
    pub highest_cn: usize,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct AssemblyResultForJson {
    pub main_haps: Vec<String>,
    pub final_haps: Vec<String>,
    pub highest_cn: usize,
}

impl AssemblyResultForJson {
    pub fn new(x: &AssemblyResult) -> Self {
        let main_haps = x
            .main_haps
            .iter()
            .map(std::string::ToString::to_string)
            .collect::<Vec<_>>();
        let final_haps = x
            .final_haps
            .iter()
            .map(std::string::ToString::to_string)
            .collect::<Vec<_>>();
        Self {
            main_haps,
            final_haps,
            highest_cn: x.highest_cn,
        }
    }
}

/// Container for storing assembled paths and associated data.
/// For now, the values are counts, but we may extend this.
/// The keys are `VString` objects, so they support random access.
#[derive(Clone, Default, Eq, PartialEq, Serialize, Deserialize)]
pub struct AssembledPaths {
    pub paths: BTreeSet<VString>,
}

impl fmt::Debug for AssembledPaths {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{{{}}}",
            itertools::intersperse(
                self.paths.iter().map(|x| format!("\"{x}\"")),
                String::from(",")
            )
            .collect::<String>()
        )
    }
}
impl AssembledPaths {
    /// Empty multiset
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Build `AssembledPaths` from a `BTreeSet`
    #[must_use]
    pub fn from_set(x: BTreeSet<VString>) -> Self {
        Self { paths: x }
    }

    /// Iterate over the keys
    pub fn path_iter(&self) -> impl Iterator<Item = &VString> {
        self.paths.iter()
    }

    /// Get a reference to inner `BTreeSet<VString>`
    #[must_use]
    pub fn inner(&self) -> &BTreeSet<VString> {
        &self.paths
    }

    /// Extract only the `BTreeSet` from `AssembledPaths`
    #[must_use]
    pub fn into_inner(self) -> BTreeSet<VString> {
        self.paths
    }

    /// Generate from a set of items which can be converted into `VString`.
    /// `Vec<u8>`, `&[u8]`, `&str`, `String`, `VStr`, and `VString` all support this.
    #[must_use]
    pub fn from_seqs<U>(iter: impl IntoIterator<Item = U>) -> Self
    where
        U: std::convert::Into<VString>,
    {
        let paths = iter
            .into_iter()
            .map(std::convert::Into::into)
            .collect::<BTreeSet<_>>();
        Self { paths }
    }
}
impl std::ops::Deref for AssembledPaths {
    type Target = BTreeSet<VString>;
    fn deref(&self) -> &Self::Target {
        &self.paths
    }
}

impl std::ops::DerefMut for AssembledPaths {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.paths
    }
}

impl std::ops::Index<usize> for AssembledPaths {
    type Output = VString;
    fn index(&self, x: usize) -> &Self::Output {
        self.iter().nth(x).unwrap_or_else(|| {
            log::warn!(
                "AssembledPaths index {} out of bounds (len {}); returning empty fallback path.",
                x,
                self.len()
            );
            empty_path_fallback()
        })
    }
}

impl std::ops::Index<&[u8]> for AssembledPaths {
    type Output = VString;
    fn index(&self, x: &[u8]) -> &Self::Output {
        self.iter().find(|seq| &seq[..] == x).unwrap_or_else(|| {
            log::warn!(
                "AssembledPaths key '{}' not found; returning empty fallback path.",
                String::from_utf8_lossy(x)
            );
            empty_path_fallback()
        })
    }
}

fn empty_path_fallback() -> &'static VString {
    static EMPTY_PATH: OnceLock<VString> = OnceLock::new();
    EMPTY_PATH.get_or_init(VString::default)
}
