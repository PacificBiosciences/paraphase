use crate::assembly::variant_graph::VGraph;
use crate::phaser::Exception;
use crate::toolkit::low_complexity::LowConfidenceSites;

use vstr::VString;

use indexmap::IndexSet;
use itertools::Itertools;
use petgraph::{stable_graph::NodeIndex, Direction::Outgoing};
use rust_htslib::{bam, htslib};
use std::{
    collections::BTreeMap,
    error,
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

pub type DError = crate::error::ParaphaseError;
pub type DResult = Result<(), DError>;

pub type MResult<T> = std::result::Result<T, Box<dyn error::Error>>;

pub type HashMap<K, V> = fnv::FnvHashMap<K, V>;
pub type HashSet<K> = fnv::FnvHashSet<K>;

/// `NotImplementedError`
/// Thrown when code that has not been implemented is executed.
#[derive(Debug, Clone, Default)]
pub struct NotImplementedError {
    pub msg: String,
}

impl std::error::Error for NotImplementedError {}

impl std::fmt::Display for NotImplementedError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "NotImplementedError{{msg: {}}}", self.msg)
    }
}

impl NotImplementedError {
    #[must_use]
    /// Create a `NotImplementedError` with a custom message.
    pub fn new(x: impl Into<String>) -> Self {
        Self { msg: x.into() }
    }
}

/// A port of `all_simple_paths` from `petgraph`, but skipping deleted nodes.
/// Necessary to support dynamic graph.
///
#[must_use]
pub fn all_simple_paths(
    graph: &VGraph,
    from: NodeIndex<u32>,
    to: NodeIndex<u32>,
    min_intermediate_nodes: usize,
    max_intermediate_nodes: Option<usize>,
) -> Vec<Vec<NodeIndex<u32>>> {
    type NodeId = NodeIndex<u32>;
    // how many nodes are allowed in simple path up to target node
    // it is min/max allowed path length minus one, because it is more appropriate when implementing lookahead
    // than constantly add 1 to length of current path
    let max_length = if let Some(l) = max_intermediate_nodes {
        l + 1
    } else {
        graph.node_count() - 1
    };

    let min_length = min_intermediate_nodes + 1;

    // list of visited nodes
    let mut visited: IndexSet<NodeId> = IndexSet::from_iter(Some(from));
    // list of childs of currently exploring path nodes,
    // last elem is list of childs of last visited node
    let mut stack = vec![graph.neighbors_directed(from, Outgoing)];

    std::iter::from_fn(move || {
        while let Some(children) = stack.last_mut() {
            // Here is where we skip the deleted nodes.
            let mut get_next = || -> Option<petgraph::stable_graph::NodeIndex<_>> {
                let mut ret = None;
                for res in children.by_ref() {
                    if let Some(weight) = graph.node_weight(res) {
                        if weight.is_del() {
                            continue;
                        }
                        ret = Some(res);
                        break;
                    }
                }
                ret
            };
            if let Some(child) = get_next() {
                if visited.len() < max_length {
                    if child == to {
                        if visited.len() >= min_length {
                            let path = visited.iter().copied().chain(Some(to)).collect::<Vec<_>>();
                            return Some(path);
                        }
                    } else if !visited.contains(&child) {
                        visited.insert(child);
                        stack.push(graph.neighbors_directed(child, Outgoing));
                    }
                } else {
                    if (child == to || children.any(|v| v == to)) && visited.len() >= min_length {
                        let path = visited.iter().copied().chain(Some(to)).collect::<Vec<_>>();
                        return Some(path);
                    }
                    stack.pop();
                    visited.pop();
                }
            } else {
                stack.pop();
                visited.pop();
            }
        }
        None
    })
    .collect::<Vec<_>>()
}

#[must_use]
/// Resolve a test-data path under `tests/data/`.
pub fn test_file(x: &str) -> std::path::PathBuf {
    let mut path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    for s in &["tests", "data", x] {
        path.push(s);
    }
    path
}

// from https://stackoverflow.com/questions/35045996/check-if-a-command-is-in-path-executable-as-process
#[must_use]
/// Find an executable by name on `PATH`.
pub fn get_program_from_path(program: &str) -> Option<std::path::PathBuf> {
    std::env::var("PATH").ok().and_then(|path| {
        path.split_terminator(':')
            .map(std::path::PathBuf::from)
            .map(|p| p.join(program))
            .find(|p| std::fs::metadata(p).is_ok())
    })
}

/// Initialize env-logger with millisecond timestamps at the requested level.
pub fn init_log(level: log::LevelFilter) {
    let mut builder = env_logger::builder();
    builder.format_timestamp_millis().filter_level(level);
    #[cfg(test)]
    {
        builder.is_test(true);
    }
    let _ = builder.try_init();
}

/// Read an indexed bam file. This can take a url or a local path.
/// If the index is not directly adjacent to the bam file, use `read_indexed_bam_with_index`.
/// # Errors
/// 1. `rust_htslib::errors::Error` if file not found or corrupted.
pub fn read_indexed_bam(
    input: impl Into<String>,
) -> Result<bam::IndexedReader, rust_htslib::errors::Error> {
    let input = input.into();
    if let Ok(url) = url::Url::parse(&input) {
        bam::IndexedReader::from_url(&url)
    } else {
        bam::IndexedReader::from_path(&input)
    }
}

/// Read an indexed BAM or CRAM file, attaching the reference when needed for CRAM decoding.
///
/// This can take a URL or a local path.
///
/// # Errors
/// 1. `rust_htslib::errors::Error` if file not found, corrupted, or CRAM reference setup fails.
pub fn read_indexed_bam_with_reference(
    input: impl Into<String>,
    reference: impl AsRef<Path>,
) -> Result<bam::IndexedReader, rust_htslib::errors::Error> {
    let mut reader = read_indexed_bam(input)?;
    let _ = reader.set_reference(reference);
    Ok(reader)
}

///
/// This comes from pysam's base-quality filtering in the pileups, which uses the next base position after a deletion.
/// in `rust-htslib`, it returns `Option<usize>` and `None` if `is_del`.
/// This function lets us get `qpos` from `bam_pileup1_t` even if the read is a deletion.
///```
/// use paraphase::toolkit::util::{raw_qpos, test_file};
/// use rust_htslib::bam::Read;
/// let mut reader =
///     rust_htslib::bam::Reader::from_path(test_file("bams/HG00733.smn1.bam")).unwrap();
/// for pileup in reader.pileup() {
///     let pileup = pileup.unwrap();
///     let qposes = pileup.alignments().map(|x| x.qpos()).collect::<Vec<_>>();
///     let qraws = pileup
///         .alignments()
///         .map(|x| raw_qpos(&x))
///         .collect::<Vec<_>>();
///     for (qp, qr) in qposes.iter().zip(qraws.iter()) {
///         assert!(
///             qp.is_none() || qp.unwrap() == *qr,
///             "qpos: {qposes:?} qr: {qraws:?}"
///         );
///     }
/// }
/// ```
#[must_use]
pub fn raw_qpos<'a>(x: &'a bam::pileup::Alignment<'a>) -> usize {
    static_assertions::const_assert_eq!(
        std::mem::size_of::<&bam::pileup::Alignment<'_>>(),
        std::mem::size_of::<&htslib::bam_pileup1_t>()
    );
    let ptr = (x as *const bam::pileup::Alignment<'_>).cast::<&htslib::bam_pileup1_t>();
    unsafe { *ptr }.qpos as usize
}

pub type GeneLevelConfig = BTreeMap<String, serde_yaml::Value>;

lazy_static::lazy_static! {
    pub static ref GIT_DESCRIBE: String =
        option_env!("VERGEN_GIT_DESCRIBE").unwrap_or("unknown").to_string();
    pub static ref FULL_VERSION: String = env!("CARGO_PKG_VERSION").to_string();
    pub static ref FULL_VERSION_PROGRAM: String =
        format!("{}-{}",
        env!("CARGO_PKG_NAME"),
        env!("CARGO_PKG_VERSION"));
    pub static ref CLI_COMMAND: String = std::env::args().collect::<Vec<_>>().join(" ");
}

#[must_use]
/// Build a BAM header with `@PG` metadata for this paraphase invocation.
pub fn output_bam_header(template: &bam::HeaderView) -> bam::Header {
    let mut header = bam::Header::from_template(template);
    let mut pg = bam::header::HeaderRecord::new(b"PG");
    pg.push_tag(b"ID", env!("CARGO_PKG_NAME"));
    pg.push_tag(b"PN", env!("CARGO_PKG_NAME"));
    pg.push_tag(b"VN", &FULL_VERSION[..]);
    pg.push_tag(b"CL", &CLI_COMMAND[..]);
    header.push_record(&pg);
    header
}

/// Read a bam file. This can take a url or a local path.
/// # Errors
/// 1. `rust_htslib::errors::Error` if file not found or corrupted.
pub fn read_bam(input: impl AsRef<Path>) -> Result<bam::Reader, String> {
    let input = input.as_ref().to_string_lossy().into_owned();
    if let Ok(url) = url::Url::parse(&input) {
        bam::Reader::from_url(&url)
    } else if input == "-" || input == "/dev/stdin" {
        bam::Reader::from_stdin()
    } else {
        bam::Reader::from_path(&input)
    }
    .map_err(|x| format!("Error: {x}"))
}

/// Read a BAM or CRAM file, attaching the reference when needed for CRAM decoding.
///
/// # Errors
/// 1. `rust_htslib::errors::Error` if file not found, corrupted, or CRAM reference setup fails.
pub fn read_bam_with_reference(
    input: impl AsRef<Path>,
    reference: impl AsRef<Path>,
) -> Result<bam::Reader, String> {
    let input = input.as_ref().to_string_lossy().into_owned();
    let mut reader = if let Ok(url) = url::Url::parse(&input) {
        bam::Reader::from_url(&url)
    } else if input == "-" || input == "/dev/stdin" {
        bam::Reader::from_stdin()
    } else {
        bam::Reader::from_path(&input)
    }
    .map_err(|x| format!("Error: {x}"))?;

    let _ = reader.set_reference(reference);
    Ok(reader)
}

/// Count records in a bam file.
///
/// Consumes the reader stream.
pub fn count_records(x: &mut impl bam::Read) -> usize {
    x.rc_records().count()
}

/// Accesses all sample names from a bam header.
/// Since a given sample can have several RG tags, we return a map
/// from sample to a vector of matching RG tags.
///
///```
/// use paraphase::toolkit::util::{sample_names, test_file};
/// use rust_htslib::bam::{Header, Read, Reader};
/// let bam = Reader::from_path(test_file("header_only.bam")).unwrap();
/// let expected_names = [("UnnamedSample", "default")]
///     .into_iter()
///     .map(|(key, val)| (key.to_string(), vec![val.to_string()]))
///     .collect::<std::collections::BTreeMap<_, _>>();
/// assert_eq!(
///     sample_names(&Header::from_template(bam.header())),
///     expected_names,
/// );
/// ```
#[must_use]
pub fn sample_names(view: &bam::Header) -> BTreeMap<String, Vec<String>> {
    let mut ret = BTreeMap::<String, Vec<_>>::new();
    if let Some(map) = view.to_hashmap().remove("RG") {
        for (sample, id) in map.iter().filter_map(|x| {
            x.get("SM")
                .and_then(|sample| x.get("ID").map(|id| (sample.clone(), id.clone())))
        }) {
            ret.entry(sample).or_default().push(id);
        }
    }
    ret
}

#[must_use]
/// Load sample names from a BAM or CRAM path by reading its header.
///
/// Returns `None` when the input cannot be opened.
pub fn sample_names_from_input(x: &PathBuf) -> Option<BTreeMap<String, Vec<String>>> {
    use bam::Read;
    let reader = read_bam(x).ok()?;
    Some(sample_names(&bam::Header::from_template(reader.header())))
}

/// Load sample names from a BAM or CRAM path by reading its header and attaching a reference for CRAM.
///
/// Returns `None` when the input cannot be opened.
pub fn sample_names_from_input_with_reference(
    x: &PathBuf,
    reference: impl AsRef<Path>,
) -> Option<BTreeMap<String, Vec<String>>> {
    use bam::Read;
    let reader = read_bam_with_reference(x, reference).ok()?;
    Some(sample_names(&bam::Header::from_template(reader.header())))
}

/*
 * BAM_CIGAR_TYPE  QUERY  REFERENCE
 * --------------------------------
 * BAM_CMATCH      1      1
 * BAM_CINS        1      0
 * BAM_CDEL        0      1
 * BAM_CREF_SKIP   0      1
 * BAM_CSOFT_CLIP  1      0
 * BAM_CHARD_CLIP  0      0
 * BAM_CPAD        0      0
 * BAM_CEQUAL      1      1
 * BAM_CDIFF       1      1
 * BAM_CBACK       0      0
 */
#[must_use]
#[inline]
pub fn consumes_ref(x: bam::record::Cigar) -> bool {
    use bam::record::Cigar::{Del, Diff, Equal, Match, RefSkip};
    matches!(x, Del(_) | RefSkip(_) | Match(_) | Diff(_) | Equal(_))
}

#[must_use]
#[inline]
pub fn consumes_qry(x: bam::record::Cigar) -> bool {
    use bam::record::Cigar::{Diff, Equal, Ins, Match, SoftClip};
    matches!(x, Ins(_) | SoftClip(_) | Match(_) | Diff(_) | Equal(_))
}

#[must_use]
/// Returns compile-time crate name (`CARGO_PKG_NAME` env).
pub fn crate_name() -> &'static str {
    env!("CARGO_PKG_NAME")
}

#[must_use]
/// Load all sequences from a faidx reader into owned uppercase byte strings.
pub fn load_all_seqs(index: &rust_htslib::faidx::Reader) -> Vec<vstr::VString> {
    load_all_seqs_view(index)
}

/// Read FASTA `(name, sequence)` pairs from a plain-text FASTA file.
///
/// When `uppercase` is true, sequence lines are uppercased before storage.
pub fn seq_name_pairs(
    path: &std::path::Path,
    uppercase: bool,
) -> Result<Vec<(VString, VString)>, DError> {
    use std::io::BufRead;

    let file = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(file);

    let mut ret = Vec::with_capacity(8);
    let mut current_name: Option<VString> = None;
    let mut current_seq = Vec::<u8>::new();

    for line in reader.lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            if let Some(name) = current_name.take() {
                ret.push((name, VString::from(current_seq.as_slice())));
                current_seq.clear();
            }
            let name = header.split_whitespace().next().unwrap_or_default();
            if name.is_empty() {
                return Err(Exception::new("FASTA header missing sequence name").into());
            }
            current_name = Some(VString::from(name));
            continue;
        }
        if line.starts_with(';') || line.trim().is_empty() {
            continue;
        }
        if current_name.is_none() {
            return Err(Exception::new("Invalid FASTA: sequence line before header").into());
        }
        current_seq.extend_from_slice(line.as_bytes());
    }

    if let Some(name) = current_name.take() {
        ret.push((name, VString::from(current_seq.as_slice())));
    }

    if ret.is_empty() {
        return Err(Exception::new("No FASTA records found").into());
    }

    if uppercase {
        for (_, seq) in &mut ret {
            seq.make_ascii_uppercase();
        }
    }

    Ok(ret)
}

#[must_use]
/// Return all contig names from an faidx reader.
pub fn faidx_names(index: &rust_htslib::faidx::Reader) -> Vec<vstr::VString> {
    (0..index.n_seqs() as i32)
        .filter_map(|v| index.seq_name(v).ok().map(VString::from))
        .collect::<Vec<_>>()
}

#[must_use]
/// Fetch all indexed sequences from an faidx reader.
pub fn load_all_seqs_view(index: &rust_htslib::faidx::Reader) -> Vec<VString> {
    (0..index.n_seqs())
        .filter_map(|x| index.seq_name(x as i32).ok())
        .filter_map(|name| index.fetch_seq(&name, 0, i64::MAX as usize).ok())
        .map(VString::from)
        .collect::<Vec<_>>()
}

#[must_use]
/// Parse homopolymer sites from a text file.
/// Format: "{pos}\t{characters}"
pub fn parse_homopolymers(path: &std::path::Path) -> LowConfidenceSites {
    let reader = File::open(path).ok().map(BufReader::new);
    LowConfidenceSites::from_map(
        reader
            .into_iter()
            .flat_map(|reader| reader.lines())
            .filter_map(Result::ok)
            .filter_map(|x| {
                let (k, v) = x.split_terminator('\t').next_tuple()?;
                let k = k.parse::<i64>().ok()?;
                let chars = v
                    .split_terminator(',')
                    .filter_map(|x| x.chars().next().map(|ch| ch as u8))
                    .collect::<std::collections::BTreeSet<_>>();
                Some((k, chars))
            })
            .collect::<BTreeMap<_, _>>(),
    )
}

pub trait DeletionInsensitiveCompare {
    /// Computes if two fingerprints are matches apart from gaps/deletions and soft-clips.
    /// Implemented for structures containing bytes.
    #[must_use]
    fn same_without_dels(&self, other: &Self) -> bool;
    #[must_use]
    fn is_del(x: u8) -> bool {
        matches!(x, b'x' | b'0')
    }
}

impl DeletionInsensitiveCompare for str {
    fn same_without_dels(&self, other: &Self) -> bool {
        self.as_bytes().same_without_dels(&other.as_bytes())
    }
}

impl DeletionInsensitiveCompare for &[u8] {
    /// This may need extension for leading xs followed by softclips,
    /// but let's hope it is rare e.g., " xxx0024".
    fn same_without_dels(&self, other: &Self) -> bool {
        let possible_match = self.len() == other.len()
            && self
                .iter()
                .zip(other.iter())
                .all(|(x, y)| x == y || Self::is_del(*x) || Self::is_del(*y));
        if !possible_match {
            return false;
        }
        let mut self_g = self
            .iter()
            .group_by(|x| *x)
            .into_iter()
            .map(|(k, v)| (*k, v.count() as i32))
            .collect::<Vec<_>>();
        let mut other_g = other
            .iter()
            .group_by(|x| *x)
            .into_iter()
            .map(|(k, v)| (*k, v.count() as i32))
            .collect::<Vec<_>>();
        let trim_ends = |x: &mut Vec<(u8, i32)>| {
            for idx in (0..x.len()).rev() {
                let reg = &x[idx];
                if !Self::is_del(reg.0) {
                    break;
                }
                x.pop();
            }
            let end_pos = x.iter().position(|x| !Self::is_del(x.0)).unwrap_or(x.len());
            x.drain(0..end_pos);
            /*
            while x.len() > 0 && Self::is_del(x.first().unwrap().0) {
                x.swap_remove(0);
            }
            */
        };
        trim_ends(&mut self_g);
        trim_ends(&mut other_g);
        let no_remaining_clips = |x: &[(u8, i32)]| -> bool {
            [x.last(), x.first()].iter().flatten().all(|x| x.0 != b'0')
        };
        [self_g, other_g].iter().all(|x| no_remaining_clips(x))
    }
}

impl<const N: usize> DeletionInsensitiveCompare for [u8; N] {
    fn same_without_dels(&self, other: &Self) -> bool {
        (&self[..]).same_without_dels(&&other[..])
    }
}

impl DeletionInsensitiveCompare for vstr::VStr<'_> {
    fn same_without_dels(&self, other: &Self) -> bool {
        (&self[..]).same_without_dels(&&other[..])
    }
}

impl DeletionInsensitiveCompare for vstr::VString {
    fn same_without_dels(&self, other: &Self) -> bool {
        self.vstr().same_without_dels(&other.vstr())
    }
}

#[cfg(test)]
mod tests {
    use crate::toolkit::util::{DResult, DeletionInsensitiveCompare};

    #[test]
    fn test_config_ok() {
        let res = crate::config::Region::load(None);
        assert_eq!(res["smn1"]["genes"].as_str(), Some("SMN1,SMN2"));
        assert_eq!(res["smn1"]["pivot_site"].as_i64(), Some(70_951_946));
        assert_eq!(res["CYP2D6"]["genes"].as_str().unwrap(), "CYP2D6");
    }

    #[test]
    fn test_gene_config_parser_ok() -> DResult {
        use std::collections::BTreeSet;
        let conf = crate::config::Gene::try_load(None)?;
        assert!(conf.genes_to_call.is_empty());
        assert_eq!(
            conf.no_vcf_genes,
            ["CFH", "CFHR3"]
                .into_iter()
                .map(String::from)
                .collect::<BTreeSet<_>>()
        );
        assert_eq!(
            conf.two_reference_regions_genes,
            ["smn1", "pms2", "strc", "ikbkg", "ncf1"]
                .into_iter()
                .map(String::from)
                .collect::<BTreeSet<_>>()
        );
        assert_eq!(
            conf.no_genome_depth_genes,
            ["pms2", "neb", "cfc1", "ikbkg", "opn1lw", "rccx"]
                .into_iter()
                .map(String::from)
                .collect::<BTreeSet<_>>()
        );
        Ok(())
    }

    #[test]
    fn same_without_dels_ok() {
        use DeletionInsensitiveCompare;
        // Simple case: matches
        assert!("12111".same_without_dels("12x11"));
        assert!(b"12111".same_without_dels(b"12x11"));
        let lhs = vstr::VString::from("121x1x11");
        let rhs = vstr::VString::from("121211x1");
        assert!(lhs.same_without_dels(&rhs));

        // Mismatched length: should be false
        let lhs = vstr::VString::from("121x1x11");
        let rhs = vstr::VString::from("121211x1a");
        assert!(!lhs.same_without_dels(&rhs));

        // Mismatched: mismatches
        //                                   *
        let lhs = vstr::VString::from("121x1x11");
        let rhs = vstr::VString::from("12121121");
        assert!(!lhs.same_without_dels(&rhs));

        // Allow softclip at end
        let lhs = vstr::VString::from("121x1x11");
        let rhs = vstr::VString::from("12121000");
        assert!(lhs.same_without_dels(&rhs));

        // Do not allow internal softclip.
        let lhs = vstr::VString::from("121x1x11");
        let rhs = vstr::VString::from("12121000");
        assert!(lhs.same_without_dels(&rhs));

        // Do allow internal gaps.
        let lhs = vstr::VString::from("121x1x11");
        let rhs = vstr::VString::from("12x21000");
        assert!(lhs.same_without_dels(&rhs));
    }
}
