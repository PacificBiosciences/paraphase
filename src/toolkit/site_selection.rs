use crate::config::Locus as LocusConfig;
use crate::phaser::Exception;
use crate::phaser::Phaser;
use crate::toolkit::deletion::Datum as DeletionDatum;
use crate::toolkit::range;
use crate::toolkit::util::{self, DError, HashMap};
use itertools::Itertools;
use rust_htslib::bam::{self, IndexedReader, Read as BamRead};
use std::fmt;
use vstr::{VStr, VString};

use std::collections::{BTreeMap, BTreeSet};
use std::hash::BuildHasherDefault;

/// `CandidateSite`
/// Holds position (`i64`), `ref_seq` and `var_seq` (both `VString`).
///
///
/// Use 0-based indexing internally but displays as 1-based to be consistent with SAM indexing.
///
/// For deletion entries, stored as `del` for `ref_seq` and `num_bases_deleted` as `var_seq`.
#[derive(Clone, Hash, PartialEq, PartialOrd, Ord, Eq, serde::Serialize, serde::Deserialize)]
pub struct CandidateSite {
    pub pos: i64,
    pub ref_seq: VString,
    pub var_seq: VString,
}

impl fmt::Display for CandidateSite {
    fn fmt(&self, format: &mut fmt::Formatter) -> fmt::Result {
        write!(format, "{}_{}_{}", self.pos + 1, self.ref_seq, self.var_seq)
    }
}
impl fmt::Debug for CandidateSite {
    fn fmt(&self, format: &mut fmt::Formatter) -> fmt::Result {
        write!(format, "{self}")
    }
}

impl CandidateSite {
    /// Create a `CandidateSite` from an integer, a reference sequence, and a variant sequence.
    #[must_use]
    pub fn new(
        pos: impl TryInto<i64>,
        ref_seq: impl Into<VString>,
        var_seq: impl Into<VString>,
    ) -> Self {
        let ref_seq = ref_seq.into();
        let var_seq = var_seq.into();
        let pos = pos.try_into().unwrap_or_else(|_| {
            log::warn!(
                "CandidateSite::new received a position that cannot be converted to i64; clamping to i64::MAX."
            );
            i64::MAX
        });
        Self {
            pos,
            ref_seq,
            var_seq,
        }
    }
    #[must_use]
    pub fn is_unit_length(&self) -> bool {
        self.ref_seq.len() == 1 && self.var_seq.len() == 1
    }

    #[must_use]
    pub fn reference_length(&self) -> usize {
        self.ref_seq.len()
    }
}

impl std::convert::From<&DeletionDatum> for CandidateSite {
    fn from(del: &DeletionDatum) -> Self {
        let pos = del.raw.start;
        let ref_seq = VString::from("del");
        let var_seq = VString::from(del.raw.len().to_string());
        CandidateSite {
            pos,
            ref_seq,
            var_seq,
        }
    }
}

impl std::str::FromStr for CandidateSite {
    type Err = Exception;
    fn from_str(x: &str) -> Result<Self, Self::Err> {
        let mut it = x.split_terminator('_');
        let (pos, ref_seq, var_seq) = it.next_tuple().ok_or_else(|| {
            Exception::new(String::from(
                "Malformatted CandidateSite string. Expected {pos}_{ref_seq}_{var_seq}",
            ))
        })?;
        let pos = pos.parse::<i64>().map_err(|e| {
            Exception::new(format!(
                "Failed to parse position (\"{pos}\") as i64. Parse error: {e:?}"
            ))
        })? - 1; // Correct for one-based indexing.
        Ok(Self {
            pos,
            ref_seq: ref_seq.into(),
            var_seq: var_seq.into(),
        })
    }
}

/// Settings for generating candidate sites for inclusion in the wfa graph.
#[derive(Clone, Debug)]
pub struct Settings {
    pub min_vaf: f64,
    pub min_read_support: i32,
    // pub regions_to_check: Vec<range::I64>,
    pub min_base_quality: u8, // Min base quality for use in fingerprinting.
    pub min_candidate_base_quality: u8, // Min base quality for use in raw pileups.
    pub min_mean_base_quality: u8, // Min mean base quality for windows for inclustion.
    pub max_indel_size: i32,
    pub max_candidate_seqs: i32,
    pub trusted_read_support: i32,
    pub permit_list: BTreeMap<i64, VString>, /* a list of positions which are permitted even if filters fail. */
    pub targeted: bool,
}

impl std::default::Default for Settings {
    fn default() -> Self {
        Self {
            min_vaf: 0.11,
            min_read_support: 5,
            // regions_to_check: vec![],
            min_base_quality: 25, // Match the default in paraphase.phaser.Phaser
            min_mean_base_quality: 25,
            min_candidate_base_quality: 13, // Match the default in pysam.Pileup
            max_indel_size: 24,
            max_candidate_seqs: 3,
            trusted_read_support: 20,
            permit_list: BTreeMap::new(),
            targeted: false,
        }
    }
}

impl Settings {
    /// Determines if an indel is small enough to use.
    #[must_use]
    pub fn indel_size_passes(&self, size: usize) -> bool {
        size as i32 <= self.max_indel_size
    }

    /// Update `site_selection::Settings` from `LocusConfig` from yaml.
    ///
    /// # Panics
    /// 1. `white_list` or `permit_list` field is not a mapping as expected.
    /// 2. key from permit list is not integral as expected.
    /// 3. value from permit list is not a string as expected.
    pub fn update_from_settings(&mut self, config: &LocusConfig) {
        if let Some(list) = config.get("white_list").or(config.get("permit_list")) {
            if let Some(list) = list.as_mapping() {
                self.permit_list = list
                    .iter()
                    .filter_map(|(k, v)| {
                        let (Some(key), Some(value)) = (k.as_i64(), v.as_str()) else {
                            return None;
                        };
                        Some((key - 1, VString::from(value))) // Subtract 1 to account for 1-based inputs.
                    })
                    .collect::<BTreeMap<_, _>>();
            }
        }
    }

    /// Build `site_selection::Settings` from `LocusConfig` from yaml.
    ///
    /// # Panics
    /// 1. `white_list` or `permit_list` field is not a mapping as expected.
    /// 2. key from permit list is not integral as expected.
    /// 3. value from permit list is not a string as expected.
    #[must_use]
    pub fn new_from_settings(config: &LocusConfig) -> Self {
        let mut ret = Self::default();
        ret.update_from_settings(config);
        ret
    }
}

#[inline]
#[must_use]
/// Return `x` with optional strand marking (reverse as lowercase).
pub fn maybe_strand_mark_char(x: u8, is_rev: bool, mark_strand: bool) -> u8 {
    if mark_strand {
        strand_mark_char(x, is_rev)
    } else {
        x.to_ascii_uppercase()
    }
}

#[inline]
#[must_use]
/// Strand-mark one base: uppercase for forward, lowercase for reverse.
pub fn strand_mark_char(x: u8, is_rev: bool) -> u8 {
    if is_rev {
        x.to_ascii_lowercase()
    } else {
        x.to_ascii_uppercase()
    }
}

/// Calculate sequences in a pileup, including indels.
/// Reproduces functionality in `pysam.libcalignedsegment.Pileup.get_query_sequences(add_indels=True)`.
///
/// Runtime improvement options: vectorize, remove array access checks.
///
///
/// # Arguments
///
/// * `x` - a Pileup to work with.
/// * `ref_seq` - reference sequence.
/// * `settings` - Settings, which specifies min base quality. paraphase's default is 13, but it is hidden inside pysam.
/// * `mark_strand` - whether or not to mark inserted sequences with strand. paraphase maps the data to upper-case.
///                   if false, all bases are upper-cased.
#[must_use]
pub fn query_seq_counter(
    x: &bam::pileup::Pileup,
    ref_seq: &[u8],
    settings: &Settings,
    offset: i64,
    aln2seq: &mut HashMap<String, VString>,
) -> (BTreeMap<VString, i32>, HashMap<VString, usize>) {
    let mut ret = BTreeMap::<VString, i32>::new();
    let mut first_seen = HashMap::<VString, usize>::default();
    let min_base_quality = settings.min_candidate_base_quality;
    let pos = x.pos();
    let offset = offset as usize;
    let mut seen_idx = 0usize;
    for aln in x.alignments() {
        let query_pos_raw = util::raw_qpos(&aln);
        let record = aln.record();
        let qname = VStr::from(record.qname());
        let is_reverse = record.is_reverse();
        let bq = record.qual().get(query_pos_raw).copied().unwrap_or(0);
        if bq < min_base_quality {
            //log::trace!("Skipping read {qname} at pos {pos} for base quality {bq} < threshold {min_base_quality}",);
            continue;
        }
        // Cache query seq, as this is can be shared across thousands of sites.
        let entry = aln2seq
            .entry(format!("{qname}:{}:{}", record.pos(), record.flags()))
            .or_insert_with(|| record.seq().as_bytes().into());
        let seq = &entry[..];
        let mut query_seq = VString::default();
        let base = if !aln.is_del() && !aln.is_refskip() {
            seq.get(query_pos_raw).copied().unwrap_or(b'N')
        } else if aln.is_refskip() {
            if is_reverse {
                b'<'
            } else {
                b'>'
            }
        } else {
            b'*'
        };
        // is_del
        /*
            if p.is_refskip:
                if bam_is_rev(p.b):
                    kputc(b'<', buf)
                else:
                    kputc(b'>', buf)
            else:
                kputc(b'*', buf)

        */
        query_seq.push(base);
        let pos = pos as usize;
        //log::trace!(
        //    "{qname}@pos={pos}. qpos: {query_pos_raw}. query_seq: {} bases, {query_seq} seq",
        //    query_seq.len(),
        //);
        match aln.indel() {
            bam::pileup::Indel::Ins(x) => {
                debug_assert!(x > 0);
                query_seq.push(b'+');
                query_seq.extend_from_slice(x.to_string().as_bytes());
                // TODO: speed this up by using slice operations to copy out faster.
                for j in 1..=(x as usize) {
                    query_seq.push(seq[j + query_pos_raw]);
                }
            }
            bam::pileup::Indel::Del(x) => {
                debug_assert!(x > 0);
                query_seq.push(b'-');
                query_seq.extend_from_slice(x.to_string().as_bytes());
                for j in 1..=(x as usize) {
                    query_seq.push(ref_seq[j + pos - offset]);
                }
            }
            bam::pileup::Indel::None => {}
        }
        query_seq.make_ascii_uppercase();
        first_seen.entry(query_seq.clone()).or_insert(seen_idx);
        seen_idx += 1;
        //log::trace!(
        //    "[query_seq_counter] Inserting query {query_seq} at pos {pos} for read name {qname} with flag {}",
        //    record.flags()
        //);
        *ret.entry(query_seq).or_default() += 1;
    }
    (ret, first_seen)
}

pub type RawVariantCounts = HashMap<i64, BTreeMap<VString, i32>>;
type RawVariantFirstSeen = HashMap<i64, HashMap<VString, usize>>;
pub type RawVariantCountsString = HashMap<i64, BTreeMap<String, i32>>;

#[must_use]
/// Convert raw variant-count map keys from `VString` to `String` for JSON/debug output.
pub fn raw_variants_to_string(
    input: &HashMap<i64, BTreeMap<VString, i32>>,
) -> HashMap<i64, BTreeMap<String, i32>> {
    input
        .iter()
        .map(|(k, v)| {
            (
                *k,
                (v.iter()
                    .map(|(hap, count)| (hap.to_string(), *count))
                    .collect::<BTreeMap<_, _>>()),
            )
        })
        .collect::<HashMap<_, _>>()
}

/// Generates a multiset of query sequences at each position in the provided range.
/// # Errors
/// 1. Failure to seek `chr:start-stop` in `IndexedReader`e
fn position_seq_counter(
    x: &mut IndexedReader,
    chr: &str,
    positions: &range::I64,
    ref_seq: &[u8],
    settings: &Settings,
    offset: i64,
    aln2seq: &mut HashMap<String, VString>,
) -> Result<(RawVariantCounts, RawVariantFirstSeen), DError> {
    let mut res: HashMap<i64, BTreeMap<VString, i32>> =
        HashMap::with_capacity_and_hasher(positions.len(), BuildHasherDefault::default());
    let mut first_seen: RawVariantFirstSeen =
        HashMap::with_capacity_and_hasher(positions.len(), BuildHasherDefault::default());
    log::trace!(
        "Fetching {chr}:{}-{}. chrnames: {:?}",
        positions.start,
        positions.end,
        x.header()
            .target_names()
            .into_iter()
            .map(VStr::from)
            .collect::<Vec<_>>()
    );
    let target_names = x
        .header()
        .target_names()
        .into_iter()
        .map(VString::from)
        .collect::<Vec<_>>();
    let chr_bytes = &chr.as_bytes();
    let chr_prefix = chr.split_terminator('_').next().unwrap_or(chr);
    let tid = target_names
        .iter()
        .position(|x| &&x[..] == chr_bytes)
        .or(target_names
            .iter()
            .position(|x| &x[..] == chr_prefix.as_bytes()))
        .ok_or_else(|| {
            DError::from(format!(
                "Failed to find chromosome for inputs {target_names:?} and query {chr}"
            ))
        })? as i32;
    let found = target_names
        .iter()
        .find(|x| &&x[..] == chr_bytes)
        .or(target_names
            .iter()
            .find(|x| &x[..] == chr_prefix.as_bytes()))
        .ok_or_else(|| {
            DError::from(format!(
                "Failed to find chromosome for inputs {target_names:?} and query {chr}"
            ))
        })?;
    log::trace!(
        "Piling up with tid = {tid}/{found} from chr {chr}:{}-{} and targets {target_names:?}",
        positions.start,
        positions.end
    );
    /*
    assert_eq!(
        x.header().target_names().len(),
        1,
        "Making sure there is only one allows us to just use tid = 0"
    );
    */
    let mut num_used = 0usize;
    let start = positions.start + 1;
    let end = positions.end + 1;
    // Paraphase uses 1-based coordinates but queries 0-based in pysam.
    // We add an offset (1) to use the same coordinates and get the same positions.
    x.fetch((tid, start, end))?;
    for x in x.pileup() {
        let x = x?;
        let pos = i64::from(x.pos());
        if pos < start {
            continue;
        }
        if pos >= end {
            break;
        }
        num_used += 1;
        let (seqs, seq_order) = query_seq_counter(&x, ref_seq, settings, offset, aln2seq);
        res.insert(pos, seqs);
        first_seen.insert(pos, seq_order);
    }
    log::debug!("Selected {num_used} candidate sites");
    Ok((res, first_seen))
}

pub type VariantMap = HashMap<i64, Vec<(VString, VString)>>;
pub type FilteredSitesForJson = (
    HashMap<i64, Vec<(String, String)>>,
    HashMap<i64, (String, String)>,
);

#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct FilteredSites {
    pub variants: VariantMap,
    pub variants_no_phasing: HashMap<i64, (VString, VString)>,
}

impl FilteredSites {
    #[must_use]
    pub fn to_json_friendly(&self) -> FilteredSitesForJson {
        let variants = self
            .variants
            .iter()
            .map(|(key, vec)| {
                (
                    *key,
                    vec.iter()
                        .map(|(x, y)| (x.to_string(), y.to_string()))
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<HashMap<_, _>>();
        let variants_no_phasing = self
            .variants_no_phasing
            .iter()
            .map(|(key, (x, y))| (*key, (x.to_string(), y.to_string())))
            .collect::<HashMap<_, _>>();
        (variants, variants_no_phasing)
    }
}

///```
/// use paraphase::toolkit::site_selection::raw_pile_to_string;
/// use vstr::VString;
/// let mut pile =
///     paraphase::toolkit::util::HashMap::<i64, std::collections::BTreeMap<VString, i32>>::default(
///     );
/// assert_eq!(raw_pile_to_string(&pile).unwrap(), "{}");
/// pile.insert(
///     0,
///     [(VString::from("ACGT"), 2)]
///         .into_iter()
///         .collect::<std::collections::BTreeMap<VString, i32>>(),
/// );
/// assert_eq!(raw_pile_to_string(&pile).unwrap(), "{0: {\"ACGT\": 2}\n}");
/// pile.entry(1).or_default().insert(VString::from("AGTG"), 1);
/// assert_eq!(
///     raw_pile_to_string(&pile).unwrap(),
///     "{0: {\"ACGT\": 2},\n1: {\"AGTG\": 1}\n}"
/// );
/// ```
pub fn raw_pile_to_string(
    x: &HashMap<i64, BTreeMap<VString, i32>>,
) -> Result<String, std::fmt::Error> {
    let x = x
        .iter()
        .map(|(k, v)| {
            (
                k,
                v.iter()
                    .map(|(vkey, vval)| (vkey.vstr(), vval))
                    .collect::<BTreeMap<_, _>>(),
            )
        })
        .collect::<BTreeMap<_, _>>();
    use std::fmt::Write;
    let mut ret = String::with_capacity(x.len() * 10);
    ret.push('{');
    let num_x = x.len();
    for (x_id, (pos, items)) in x.iter().enumerate() {
        write!(ret, "{pos}:")?;
        ret.push_str(" {");
        ret.push_str(
            &Itertools::intersperse(
                items
                    .iter()
                    .map(|(seq, count)| format!("\"{seq}\": {count}")),
                String::from(","),
            )
            .collect::<String>(),
        );
        ret.push_str(if x_id == num_x - 1 { "}\n" } else { "},\n" });
    }
    ret.push('}');
    Ok(ret)
}

///```
/// assert!(paraphase::toolkit::site_selection::seq_is_indel(b"BG-"));
/// assert!(!paraphase::toolkit::site_selection::seq_is_indel(b"BG"));
/// assert!(!paraphase::toolkit::site_selection::seq_is_indel(
///     b"ACGTG_1_3"
/// ));
/// ```
#[must_use]
pub fn seq_is_indel(x: &[u8]) -> bool {
    x.iter().any(|x| matches!(x, b'+' | b'-'))
}

///
/// Filter raw sites + select variant sites.
///
/// # Errors
/// 1. Faidx fetch failure.
/// 2. Failure to make local faidx.
/// 3. Error processing indel. (Should not happen.)
fn filtered_sites(
    settings: &Settings,
    phaser: &mut Phaser,
    raw_piles: &HashMap<i64, BTreeMap<VString, i32>>,
    raw_pile_first_seen: &RawVariantFirstSeen,
    regions_to_check: &[range::I64],
) -> Result<FilteredSites, DError> {
    let del_str = VString::from(&[b'*'][..]);
    let offset = phaser.offset();
    log::debug!(
        "Site-selection threshold: min_read_support={}",
        settings.min_read_support
    );
    //log::trace!("Initial raw pileup: {}", raw_pile_to_string(raw_piles)?);
    log::trace!(
        "Low-complexity masked sites at start of get_candidate_pos: {:?}",
        phaser.low_complexity_sites
    );
    let full_ref_seq = {
        log::debug!("Loading full locus reference sequence for candidate-site selection");
        let faidx = phaser.make_local_faidx()?;
        log::debug!(
            "Loading seq {:?} for faidx at {:?}",
            phaser.local_chr(),
            &phaser.local_reference()
        );
        for id in 0..faidx.n_seqs() {
            log::debug!(
                "Available faidx sequence[{id}] = {}",
                faidx.seq_name(id as i32)?
            );
        }
        let local_chr = phaser.local_chr().ok_or_else(|| {
            crate::phaser::Exception::new(format!(
                "No local chromosome parsed from realign region '{}'",
                phaser.realign_region
            ))
        })?;
        faidx
            .fetch_seq(local_chr, 0, i32::MAX as usize)?
            .to_ascii_uppercase()
    };
    log::debug!(
        "Loaded full locus reference sequence: length={}",
        full_ref_seq.len()
    );
    let mut variants: HashMap<i64, Vec<(VString, VString)>> =
        HashMap::with_capacity_and_hasher(8, BuildHasherDefault::default());
    let mut variants_no_phasing: HashMap<i64, (VString, VString)> =
        HashMap::with_capacity_and_hasher(8, BuildHasherDefault::default());
    let faidx = phaser.make_faidx()?; // Notice that this is global/genomic, not local.
                                      // Homozygous
                                      // May need to save more/less than 32768.
    let chr = phaser.chr().ok_or_else(|| {
        crate::phaser::Exception::new(format!(
            "No chromosome available for genomic faidx fetch in gene '{}'",
            phaser.gene_name()
        ))
    })?;
    log::trace!(
        "Fetching seq {chr} from faidx with names {:?}",
        util::faidx_names(&faidx)
    );
    let cached_faidx_seq = faidx
        .fetch_seq(chr, 0, i64::MAX as usize)?
        .to_ascii_uppercase();
    for (pos, pileup) in raw_piles {
        let pos = *pos;
        let depth = pileup.values().sum::<i32>();
        let del_count = pileup.get(&del_str).copied().unwrap_or(0);
        log::trace!(
            "Pileup depth summary at position {pos}: del_count={del_count}, total_depth={depth:?}"
        );
        assert!(
            depth >= del_count,
            "depth {depth} should be >= del_count {del_count}"
        );
        let depth_without_dels = depth - del_count;
        //log::trace!("Depth excluding deletions at position {pos}: {depth_without_dels}");
        let offset_pos = (pos - offset) as usize;
        let pass_for_dels = del_count < settings.min_read_support || phaser.allow_del_bases(pos);
        if depth < settings.min_read_support || !pass_for_dels {
            log::trace!(
                "variant filtered out due to low depth or high deletion count. depth {depth}, deletion count {del_count}, min_read_support {}",
                settings.min_read_support,
            );
            continue;
        }
        let ref_seq = VString::from(full_ref_seq[offset_pos..=offset_pos].to_ascii_uppercase());
        assert!(
            ref_seq.len() == 1,
            "ref seq {ref_seq:?} should be of size 1 because of the fetch."
        );
        // We are avoiding building a counter object.
        // Technically this could be faster using a heap, but we don't expect many variants at a site and so this is fair.
        let first_seen = raw_pile_first_seen.get(&pos);
        let most_common = pileup
            .iter()
            .sorted_by(|(item_a, count_a), (item_b, count_b)| {
                let by_count = count_b.cmp(count_a);
                if by_count != std::cmp::Ordering::Equal {
                    return by_count;
                }
                // Python Counter tie behavior: preserve first encounter order.
                let a_idx = first_seen
                    .and_then(|m| m.get(*item_a))
                    .copied()
                    .unwrap_or(usize::MAX);
                let b_idx = first_seen
                    .and_then(|m| m.get(*item_b))
                    .copied()
                    .unwrap_or(usize::MAX);
                a_idx.cmp(&b_idx)
            })
            .filter(|x| &x.0[..] != b"*")
            .map(|(x, c)| (x.vstr(), c))
            .take(settings.max_candidate_seqs as usize)
            .collect::<Vec<_>>();
        assert!(
            most_common.windows(2).all(|x| x[1].1 <= x[0].1),
            "Most common should be sorted in descending order: {most_common:?}"
        );
        let counter_len = pileup.len() - usize::from(pileup.contains_key(&del_str));
        let total_depth_above_minimum = f64::from(depth_without_dels - settings.min_read_support);
        // TODO: eventually replace with a cdf.
        // use `paraphase::util::site_probably_het`

        let depth_without_dels_threshold = f64::from(depth_without_dels) * 0.85;
        let most_common_threshold = if total_depth_above_minimum
            .partial_cmp(&depth_without_dels_threshold)
            .unwrap_or(std::cmp::Ordering::Equal)
            == std::cmp::Ordering::Greater
        {
            total_depth_above_minimum
        } else {
            depth_without_dels_threshold
        };
        log::trace!(
            "Most-common allele thresholding at position {pos}: threshold={most_common_threshold}, counts={most_common:?}, top_allele={:?}",
            most_common[0]
        );
        let is_homozygous = counter_len == 1
            || (counter_len >= 2 && f64::from(*most_common[0].1) > most_common_threshold);

        // Warning:
        // We make low complexity sites match exactly
        let hp_site_data = phaser.low_complexity_sites.get_0based(pos);
        let is_hp_site = hp_site_data.is_some();
        let allow_del = phaser.allow_del_bases(pos);
        if is_homozygous {
            //log::trace!(
            //    "Classified homozygous candidate site: offset_pos={offset_pos}, genomic_pos={pos}"
            //);
            let var_seq: VStr<'_> = most_common[0].0;
            if var_seq != ref_seq[..] {
                if !seq_is_indel(&var_seq) {
                    // Process SNV
                    // Homozygous sites in deletions are added as heterozygous sites.
                    log::trace!(
                        "Processing homozygous substitution candidate at position {pos}: var={var_seq}, ref={}, allow_del={}",
                        VStr::from(&ref_seq[..]),
                        allow_del
                    );
                    if allow_del && del_count >= settings.min_read_support && !is_hp_site {
                        log::trace!(
                            "is not a hp site. del count: {del_count}. Allow del bases: {}",
                            allow_del
                        );
                        variants
                            .entry(pos)
                            .or_default()
                            .push((ref_seq, VString::from(var_seq)));
                    } else if hp_site_data
                        .map_or(true, |x| !(var_seq.len() == 1 && x.contains(&var_seq[0])))
                    {
                        phaser
                            .hom_sites
                            .push(CandidateSite::new(pos, ref_seq, var_seq));
                    }
                } else if !is_hp_site {
                    // Process indel
                    log::trace!("Processing homozygous indel candidate at position {pos}");
                    let (processed_refseq, processed_varseq, indel_len) =
                        phaser.process_indel(pos, &ref_seq, &var_seq, &cached_faidx_seq)?;
                    if indel_len <= settings.max_indel_size {
                        phaser.hom_sites.push(CandidateSite::new(
                            pos,
                            processed_refseq,
                            processed_varseq,
                        ));
                    }
                }
            }
        } else if counter_len >= 2 {
            log::trace!("Classified potential heterozygous site: offset_pos={offset_pos}, genomic_pos={pos}");
            let found_ref = pileup
                .iter()
                .map(|(x, c)| (x.vstr(), c))
                .any(|x| x.0[..] == ref_seq[..]);
            //let found_ref = most_common.iter().any(|x| x.0[..] == ref_seq[..]);
            if found_ref || settings.permit_list.contains_key(&pos) {
                log::trace!(
                    "Reference allele observed among most-common pileup bases at position {pos}"
                );
                for (var_seq, _var_count) in most_common.iter().filter(|(seq, count)| {
                    let count = **count;
                    let is_refseq = seq[..] == ref_seq[..];
                    let sufficient_read_count = count >= settings.min_read_support;
                    let sufficient_vaf = f64::from(count) >= settings.min_vaf * f64::from(depth);
                    let count_above_trusted = count >= settings.trusted_read_support;
                    log::trace!(
                        "Evaluating heterozygous candidate at position {pos}: var={seq}, ref={ref_seq}, sufficient_count={sufficient_read_count}, sufficient_vaf={sufficient_vaf}, count={count}, trusted_read_support={}, min_read_support={}",
                        settings.trusted_read_support,
                        settings.min_read_support
                    );
                    !is_refseq && ((sufficient_read_count && sufficient_vaf) || (!settings.targeted && count_above_trusted))
                }) {
                    // Substitution
                    if !seq_is_indel(var_seq) {
                        log::trace!("SNV candidate classification at position {pos}: var_seq={var_seq}, is_hp_site={is_hp_site}");
                        debug_assert!(!var_seq.iter().any(|x| matches!(x, b'+' | b'-')));
                        if is_hp_site {
                            let var_seq_prohibited = var_seq.len() == 1
                                && hp_site_data
                                    .as_ref()
                                    .map_or(false, |forbid| forbid.contains(&var_seq[0]));
                            log::trace!(
                                "Homopolymer-context SNV check at position {pos}: var_seq={var_seq}, var_seq_prohibited={var_seq_prohibited}, forbidden_bases={:?}",
                                hp_site_data.as_ref().map_or(vec![], |x| x.iter().copied().collect::<Vec<_>>())
                            );
                            log::trace!(
                                "Neighbor homopolymer constraints at position {pos}: left={:?}, right={:?}",
                                if pos > 0 { phaser.low_complexity_sites.get(&(pos - 1)) } else { None },
                                phaser.low_complexity_sites.get(&(pos + 1))
                            );
                            if !var_seq_prohibited {
                                if hp_site_data.as_ref().map_or(false, |x| x.contains(&b'1')) {
                                    log::trace!("Homopolymer SNV accepted for phasing at position {pos}: var_seq={var_seq}");
                                    variants
                                        .entry(pos)
                                        .or_default()
                                        .push((ref_seq.clone(), (*var_seq).into()));
                                } else {
                                    log::trace!(
                                        "Homopolymer SNV routed to non-phasing bucket at position {pos}: var_seq={var_seq}, hp_site_data={hp_site_data:?}"
                                    );
                                    variants_no_phasing.entry(pos).or_insert_with(|| {
                                        ((&ref_seq[..]).into(), (&var_seq[..]).into())
                                    });
                                }
                            }
                        } else {
                            log::trace!("SNV candidate added to phasing variants at position {pos}: var_seq={var_seq}");
                            variants
                                .entry(pos)
                                .or_default()
                                .push(((&ref_seq[..]).into(), (&var_seq[..]).into()));
                        }
                    } else if !is_hp_site {
                        assert!(seq_is_indel(var_seq));
                        // Indel
                        let (ref_seq, var_seq, indel_len) =
                            phaser.process_indel(pos, &ref_seq, var_seq, &cached_faidx_seq)?;
                        if indel_len <= settings.max_indel_size {
                            variants_no_phasing
                                .entry(pos)
                                .or_insert_with(|| (ref_seq, var_seq));
                        }
                    } else {
                        log::trace!(
                            "Position {pos} is in a homopolymer context; skipping indel candidate handling."
                        );
                    }
                }
            }
        }
        // Now handle exclusion logic.
    }

    let excluded_variants = regions_to_check
        .iter()
        .flat_map(|x| {
            variants
                .keys()
                .filter(|pos| x.contains(*pos) && **pos > x.start)
        })
        .collect::<BTreeSet<_>>();
    log::trace!(
        "Excluding {}/{} variant sites. Variants: {variants:?}",
        excluded_variants.len(),
        variants.len()
    );

    for (pos, variants) in variants.iter().filter(|x| !excluded_variants.contains(x.0)) {
        if variants.len() == 1 {
            log::trace!(
                "Position {pos} has one passing variant and is retained as a candidate site: {variants:?}"
            );
            let var = &variants[0];
            phaser
                .candidate_sites
                .insert(CandidateSite::new(*pos, &var.0, &var.1));
        } else if let Some(permitted_variant) = settings.permit_list.get(pos) {
            log::trace!("Permitted variant override encountered at position {pos}; validating candidate list");
            if let Some((rseq, vseq)) = variants
                .iter()
                .find(|(_rseq, vseq)| vseq != permitted_variant)
            {
                phaser
                    .candidate_sites
                    .insert(CandidateSite::new(*pos, rseq, vseq));
            }
        } else {
            log::trace!(
                "Multiple competing variants remain at position {pos}; dropping this position from selected candidates: variants={variants:?}"
            );
        }
    }

    let excluded_variants = regions_to_check
        .iter()
        .flat_map(|x| {
            variants_no_phasing
                .keys()
                .filter(|pos| x.contains(*pos) && **pos > x.start)
        })
        .sorted()
        .unique()
        .collect::<BTreeSet<_>>();

    log::debug!(
        "Excluding {}/{} no-phasing variant sites.",
        excluded_variants.len(),
        variants_no_phasing.len()
    );

    for site in variants_no_phasing
        .iter()
        .filter(|x| !excluded_variants.contains(&x.0))
        .map(|(pos, variant)| CandidateSite::new(*pos, &variant.0, &variant.1))
    {
        phaser.het_sites_no_phasing.push(site);
    }

    log::debug!(
        "{} het sites no phasing from {} initial variants_no_phasing positions",
        phaser.het_sites_no_phasing.len(),
        variants_no_phasing.len()
    );
    log::debug!(
        "Homozygous site count retained for downstream use: {}",
        phaser.hom_sites.len()
    );
    phaser.het_sites = phaser
        .candidate_sites
        .iter()
        .sorted()
        .cloned()
        .collect::<Vec<_>>();
    log::debug!(
        "{} candidate sites. Het sites for phasing {:?}",
        phaser.candidate_sites.len(),
        phaser.candidate_sites
    );
    Ok(FilteredSites {
        variants,
        variants_no_phasing,
    })
}

/// Performs a pileup.
/// WARNING: positions are 0-based here to match htslib.
/// We are consistently using 0-based, which does not match paraphase. We are doing book-keeping to keep it in sync outside of this function.
/// # Errors
/// 1. Failure to fetch positions in `position_seq_counter`
/// 2. Errors in `filtered_sites`.
///   a. Faidx fetch failure.
///   b. Failure to make local faidx.
///   c. Error processing indel. (Should not happen.)
/// Run candidate-site discovery from read pileups for one phaser instance.
///
/// Returns phased candidates and non-phased variant candidates together with
/// updated phaser state (hom/het site sets).
pub fn pileup(
    x: &mut IndexedReader,
    positions: &range::I64,
    ref_seq: &[u8],
    phaser: &mut Phaser,
    settings: Option<&Settings>,
    regions_to_check: &[range::I64],
    chr: &str,
) -> Result<(FilteredSites, RawVariantCounts), DError> {
    let default = Settings::default();
    let settings = settings.unwrap_or(&default);
    let offset = phaser.offset();
    let mut aln2seq =
        HashMap::<String, VString>::with_capacity_and_hasher(256, BuildHasherDefault::default());
    let (raw_piles, raw_pile_first_seen) =
        position_seq_counter(x, chr, positions, ref_seq, settings, offset, &mut aln2seq)?;
    log::debug!(
        "Pileup context: region_chr={chr}, local_chr={:?}, genome_chr={:?}, aligned_locus_bam={:?}",
        phaser.local_chr(),
        phaser.chr(),
        phaser.realigned_bam_path(),
    );
    //log::trace!("Raw piles: {raw_piles:?}");
    let f = filtered_sites(
        settings,
        phaser,
        &raw_piles,
        &raw_pile_first_seen,
        regions_to_check,
    )?;
    Ok((f, raw_piles))
}

#[cfg(test)]
mod tests {
    use super::strand_mark_char;
    #[test]
    fn strand_mark_char_ok() {
        assert_eq!(b'A', strand_mark_char(b'A', false));
        assert_eq!(b'a', strand_mark_char(b'a', true));
        assert_eq!(b'A', strand_mark_char(b'A', false));
        assert_eq!(b'a', strand_mark_char(b'a', true));
        assert_eq!(b'G', strand_mark_char(b'G', false));
        assert_eq!(b'g', strand_mark_char(b'g', true));
        assert_eq!(b'G', strand_mark_char(b'G', false));
        assert_eq!(b'g', strand_mark_char(b'g', true));
        for base in b"1234567890QRS".iter().copied() {
            assert_eq!(base.to_ascii_lowercase(), strand_mark_char(base, true));
            assert_eq!(base.to_ascii_uppercase(), strand_mark_char(base, false));
        }
    }
}
