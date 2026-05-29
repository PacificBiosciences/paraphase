use crate::assembly::assembly_result::{AssembledPaths, AssemblyResult};
use crate::assembly::graph_stats::count_suffix_x;
use crate::assembly::variant_graph::{Graph, SegmentationClass};
use crate::toolkit::range::Range;
use crate::toolkit::util::{DError, DResult};

use anyhow::anyhow;
use itertools::Itertools;
use std::collections::{BTreeMap, BTreeSet};
use vstr::{VStr, VString};

type SubregionAssemblyResult = (
    bool,
    Vec<VString>,
    BTreeMap<VString, BTreeSet<VString>>,
    BTreeMap<VString, BTreeSet<VString>>,
);

impl Graph {
    /// Merge blocks together in order of decreasing size.
    pub fn merge_blocks_by_size(&mut self) -> Result<(), DError> {
        let mut iternum = 0usize;
        log::debug!(
            "[merge_blocks_by_size] self.pos_edge_map {:?}",
            self.pos_edge_map
        );
        while !self.pos_edge_map.is_empty() {
            iternum += 1;
            let order = self.order_blocks_by_size();
            log::trace!("Order: {order:?} at iteration {iternum} for merge_blocks_by_size");
            if order.len() == 1 {
                log::trace!(
                    "Only one block remains. Break!. Order: {order:?}. pem: {:?}",
                    self.formatted_edges()
                );
                break;
            }
            let mut i = 0;
            let mut total_success = false;
            for region in &order {
                let pos1 = region[0];
                let pos2 = region[1];
                i += 1;
                log::trace!("About to assemble subregion {pos1}..={pos2}");
                let (success, haps, _x, _y) = self.subregion_assembly(pos1, pos2, false)?;
                log::trace!("Assembled subregion {pos1}..={pos2}");
                if success {
                    log::trace!(
                        "Success for {pos1}..={pos2} yielded haps {haps:?} at iternum {iternum}. Current state before rm_add_edges {pos1}..={pos2}: {:?}", self.formatted_edges()
                    );
                    self.rm_add_edges(pos1, pos2, haps.iter())?;
                    log::trace!(
                        "Success for {pos1}..={pos2} yielded haps {haps:?} at iternum {iternum}. Current state after rm_add_edges {pos1}..={pos2}: {:?}", self.formatted_edges()
                    );
                    total_success = true;
                    break;
                }
                log::trace!("At iternum {iternum}, subregion assembly for {pos1}..={pos2} did not have success.");
            }
            if i == order.len() && !total_success {
                break;
            }
            log::trace!(
                "Success: {total_success}. Remaining in pos edge map: {:?} at iternum {iternum}",
                self.pos_edge_map.len()
            );
        }
        Ok(())
    }

    /// Sort blocks by size.
    /// Returns sorted vector of (node1, node2, minlen, maxlen)
    /// ordered by increasing (minlen, maxlen)
    #[must_use]
    pub(crate) fn order_blocks_by_size(&self) -> Vec<[u32; 4]> {
        let lengths = self
            .node_iter()
            .map(|x| (x.pos, x.hap.len() as u32))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        lengths
            .windows(2)
            .map(|window| {
                let len1 = window[0].1;
                let len2 = window[1].1;
                let node1 = window[0].0;
                let node2 = window[1].0;
                let minlen = std::cmp::min(len1, len2);
                let maxlen = std::cmp::max(len1, len2);
                [node1, node2, minlen, maxlen]
            })
            .sorted_by_key(|x| (-(x[2] as i32), -(x[3] as i32)))
            .collect::<Vec<_>>()
    }

    /// Break the graph up into easier/harder segments to assemble separately.
    pub fn segment_graph(&self) -> Result<BTreeMap<Range<u32>, SegmentationClass>, DError> {
        use SegmentationClass::{Ten, Two};
        log::trace!("dnhap before segmentation: {:?}", self.dnhap);
        let mut dnhap = self.dnhap.clone();
        let keys = dnhap.keys().copied().collect::<Vec<u32>>();
        let values = dnhap.values().copied().collect::<Vec<_>>();
        assert_eq!(values.len(), keys.len());
        values
            .windows(3)
            .zip(keys.iter().copied().skip(1))
            .enumerate()
            .for_each(|(window_index, (v, k))| {
                if v == [Ten, Two, Ten] {
                    log::trace!(
                        "Detected 10-2-10 pattern at window index {window_index} (centroid {}). Updating range {window_index}..={} to Ten. k: {k}",
                        window_index + 1,
                        window_index + 2
                    );
                    dnhap.insert(k, Ten);
                }
            });
        assert!(keys.windows(2).all(|x| x[1] >= x[0]));
        log::trace!("Updated local segmentation classes: {dnhap:?}");

        let mut segment_pos = 0u32;
        Ok(dnhap
            .iter()
            .group_by(|x| *x.1)
            .into_iter()
            .map(|(key, group)| {
                let group_len = group.count() as u32;
                let segment_end = segment_pos + group_len;
                let range = Range::new(segment_pos, segment_end - 1);
                segment_pos = segment_end;
                (range, key)
            })
            .collect::<BTreeMap<_, _>>())
    }

    /// Given strings to the left and right and connections between sub-haplotypes,
    /// generate additional possible complete haplotypes.
    #[must_use]
    pub(crate) fn rescue_missing<'a>(
        left: &[&'a VString],
        right: &[&'a VString],
        dnext: &BTreeMap<VString, BTreeSet<VString>>,
        dbefore: &BTreeMap<VString, BTreeSet<VString>>,
    ) -> BTreeSet<VString> {
        let sub_rescue = |side: &[&VString],
                          other_side: &[&VString],
                          match_dict: &BTreeMap<VString, BTreeSet<VString>>,
                          is_left: bool|
         -> BTreeSet<VString> {
            let merge_slices_base = |x: &VString, y: &VString| -> VString {
                let mut ret = x.clone();
                ret.extend_from_slice(&y[..]);
                ret
            };
            let merge_slices = |x: &VString, y: &VString, is_left: bool| -> VString {
                if is_left {
                    merge_slices_base(x, y)
                } else {
                    merge_slices_base(y, x)
                }
            };
            side.iter()
                .filter_map(|x| match_dict.get(*x).map(|matches| (x, matches)))
                .filter(|x| bytecount::count(x.0, b'x') <= 1)
                .filter_map(|(side_missing, next_hits)| {
                    let next_hit = next_hits.iter().next();
                    if next_hits.len() == 1
                        && (other_side.is_empty() || other_side[..] == [next_hit?])
                    {
                        Some(merge_slices(side_missing, next_hit?, is_left))
                    } else {
                        let cand_missing = other_side
                            .iter()
                            .filter(|x| next_hits.contains(**x))
                            .collect::<Vec<_>>();
                        if cand_missing.len() == 1
                            && other_side.len() == 1
                            && &other_side[0] == cand_missing[0]
                            && side.len() == 1
                        {
                            log::trace!(
                        "Merging {side_missing:?}. {next_hits:?} and {cand_missing:?} from {}", if is_left {"left"} else {"right"}
                    );
                            Some(merge_slices(side_missing, cand_missing[0], is_left))
                        } else {
                            None
                        }
                }
                }).collect::<BTreeSet<_>>()
        };
        let mut ret = sub_rescue(left, right, dnext, true);
        for item in sub_rescue(right, left, dbefore, false) {
            ret.insert(item);
        }
        log::trace!("After rescuing: {ret:?}");
        ret
    }

    /// Find potential partial haplotypes that match haps at left/right
    /// to use for assembling.
    #[must_use]
    pub(crate) fn get_missing<'a>(
        sub_haps_assembled: &BTreeSet<VString>,
        left: &'a [VString],
        right: &'a [VString],
    ) -> (Vec<&'a VString>, Vec<&'a VString>) {
        let left = left
            .iter()
            .filter(|x| {
                let len = x.len();
                !sub_haps_assembled
                    .iter()
                    .any(move |sub| sub.len() >= len && sub[..len] == x[..])
            })
            .collect::<Vec<&VString>>();
        let right = right
            .iter()
            .filter(|x| {
                let len = x.len();
                !sub_haps_assembled
                    .iter()
                    .any(move |sub| sub.len() >= len && sub[sub.len() - len..] == x[..])
            })
            .collect::<Vec<&VString>>();
        log::trace!(
            "Finding unmatched left/right haps using assembled sub-haps {sub_haps_assembled:?}. Unmatched left: {left:?}. Unmatched right: {right:?}"
        );
        (left, right)
    }

    /// Assemble region in pos1..=pos2
    pub fn subregion_assembly(
        &mut self,
        pos1: u32,
        pos2: u32,
        allow_x: bool,
    ) -> Result<SubregionAssemblyResult, DError> {
        let mut success = true;
        let mut curr_pos1 = pos1;
        let mut pos1_haps = self.get_owned_haps_by_pos(pos1);
        log::trace!("[subregion_assembly] Initial pos1_haps: {pos1_haps:?} for {pos1}..={pos2}. Pos edge map: {:?}", self.formatted_edges());
        if pos1_haps.is_empty() {
            return Err(anyhow::anyhow!(
                "pos1 haps is empty at {pos1}..={pos2}. Allow x: {allow_x}"
            )
            .into());
        }
        let mut dnext = BTreeMap::new();
        let mut dbefore = BTreeMap::new();
        let mut iternum = 1usize;
        while curr_pos1 < pos2 {
            let next_pos = self.get_next_pos(curr_pos1).ok_or_else(|| {
                anyhow::anyhow!(
                    "Missing next position: {curr_pos1} query, initial pos {pos1} and pos2 {pos2}"
                )
            })?;
            let next_haps = self.get_owned_haps_by_pos(next_pos);
            if ![&pos1_haps, &next_haps].iter().all(|x| !x.is_empty()) {
                return Err(anyhow::anyhow!(
                    "pos1 haps or next haps are empty. {pos1_haps:?}, {next_haps:?}"
                )
                .into());
            }
            let mut subregion_haps_assembled = BTreeSet::new();
            if pos1_haps.len() == 1 || next_haps.len() == 1 {
                for (hap1, hap2) in itertools::iproduct!(pos1_haps.iter(), next_haps.iter()) {
                    let mut new_hap = VString::from(hap1);
                    new_hap.extend_from_slice(hap2);
                    subregion_haps_assembled.insert(new_hap);
                }
                dnext.clear();
                dbefore.clear();
            } else {
                let new_nodes = [(pos1, &pos1_haps), (next_pos, &next_haps)]
                    .into_iter()
                    .map(|(k, value)| (k, value.iter().map(|x| x.vstr()).collect::<BTreeSet<_>>()))
                    .collect::<BTreeMap<_, _>>();
                let (sub_hap, dnext_new, dbefore_new) = self.merge_two_pos(&new_nodes)?;
                dnext = dnext_new
                    .into_iter()
                    .map(|(k, v)| {
                        (
                            VString::from(k),
                            v.into_iter().map(VString::from).collect::<BTreeSet<_>>(),
                        )
                    })
                    .collect::<BTreeMap<_, _>>();
                dbefore = dbefore_new
                    .into_iter()
                    .map(|(k, v)| {
                        (
                            VString::from(k),
                            v.into_iter().map(VString::from).collect::<BTreeSet<_>>(),
                        )
                    })
                    .collect::<BTreeMap<_, _>>();
                let thres = Self::get_len_threshold(pos1_haps[0].len(), next_haps[0].len());
                for (hap, _reads) in sub_hap.iter().filter(|(_hap, reads)| reads.len() >= thres) {
                    subregion_haps_assembled.insert(hap.clone());
                }
            }
            let (mut left, mut right) =
                Self::get_missing(&subregion_haps_assembled, &pos1_haps, &next_haps);
            if !left.is_empty() || !right.is_empty() {
                if !allow_x {
                    success = false;
                    return Ok((success, pos1_haps, dnext, dbefore));
                }
                let hap_len1 = pos1_haps[0].len();
                let hap_len2 = next_haps[0].len();
                if std::cmp::min(hap_len1, hap_len2) >= 3
                    && matches!(
                        self.dnhap
                            .get(&(next_pos - 1))
                            .map(SegmentationClass::is_assigned),
                        Some(true)
                    )
                    && self
                        .edge_info
                        .get(next_pos as usize - 1)
                        .and_then(|x| x.iter().copied().max())
                        .unwrap_or(0)
                        >= 2
                {
                    let mut missing = Self::rescue_missing(&left, &right, &dnext, &dbefore);
                    subregion_haps_assembled.append(&mut missing);
                    (left, right) =
                        Self::get_missing(&subregion_haps_assembled, &pos1_haps, &next_haps);
                }

                for missing in &right {
                    let mut hap = vec![b'x'; hap_len1];
                    hap.extend_from_slice(&missing[..]);
                    subregion_haps_assembled.insert(VString::from(hap));
                }
                let x_missing = vec![b'x'; hap_len2];
                for missing in &left {
                    let mut hap = (*missing).clone();
                    hap.extend_from_slice(&x_missing);
                    subregion_haps_assembled.insert(hap);
                }
            }
            pos1_haps = subregion_haps_assembled
                .into_iter()
                .collect::<Vec<VString>>();
            curr_pos1 = next_pos;
            iternum += 1;
            let _ = iternum;
        }
        Ok((success, pos1_haps, dnext, dbefore))
    }

    #[must_use]
    /// Support threshold for accepting merged sub-haplotype pairs.
    ///
    /// Short blocks require stronger support to avoid over-merging noisy links.
    fn get_len_threshold(x: usize, y: usize) -> usize {
        let maxlen = std::cmp::max(x, y);
        let minlen = std::cmp::min(x, y);
        if minlen > 2 || (minlen == 2 && maxlen >= 10) {
            1
        } else {
            2
        }
    }

    /// Extend a partially unknown haplotype (`x` placeholders) using adjacent
    /// assembled-block evidence from both sides when confident.
    fn extended_hap(&mut self, hap: VStr<'_>) -> Result<VString, DError> {
        if !hap.contains(&b'x') {
            return Ok(VString::from(hap));
        }
        let prefix_x_count = hap.iter().position(|x| *x != b'x').unwrap_or(hap.len());
        log::trace!("Hap: {hap} has length {}", hap.len());
        let suffix_x_count = count_suffix_x(&hap);
        log::trace!("pref x count: {prefix_x_count}. suffix: {suffix_x_count} for hap {hap}");
        let mut identical_bases_right = VString::default();
        let mut identical_bases_left = VString::default();
        if suffix_x_count < (self.nvar - 2) {
            log::trace!(
                "x count {suffix_x_count} < nvar - 2 ({}) for {hap}. dnhap: {:?}",
                self.nvar - 2,
                self.dnhap
            );
            let next_pos = (suffix_x_count + 1) as u32;
            let dnhap_is_present_and_assigned = matches!(
                self.dnhap
                    .get(&(suffix_x_count as u32))
                    .map(SegmentationClass::is_assigned),
                Some(true)
            );
            let min_support = self
                .edge_info
                .get(suffix_x_count)
                .and_then(|x| x.iter().copied().min())
                .unwrap_or(0);
            log::trace!(
                "dnhap_is_present_and_assigned on right for {hap}. pos edge map has minimum {min_support}. Pos edge map at {suffix_x_count}",
            );
            if self.pos_edge_map.contains_key(&suffix_x_count) {
                log::trace!(
                    "[extended_hap] self.pos_edge_map[&suffix_x_count] {:?}",
                    self.pos_edge_map[&suffix_x_count]
                );
            }
            if dnhap_is_present_and_assigned && min_support >= 2 {
                log::trace!("dnhap had the suffix, it was assigned, and pos_edge_map had at least min count of 2 for hap {hap}.");
                let previous_pos = self.get_previous_pos(next_pos).ok_or_else(|| anyhow::anyhow!("Failed to get previous pos for {next_pos} with suffix x count {suffix_x_count}"))?;
                log::trace!("Assembling subregion from {previous_pos}..={next_pos}");
                let (_, _, dnext, _dbefore) =
                    self.subregion_assembly(previous_pos, next_pos, false)?;
                log::trace!("dnext {dnext:?} and dbefore {_dbefore:?}");
                let block = VString::from(&hap[(previous_pos as usize)..=suffix_x_count]);
                if let Some(candidates) = dnext.get(&block).map(|candidates| {
                    log::trace!("candidates from dnext: {candidates:?}");
                    if candidates.len() <= 1 {
                        candidates.iter().cloned().collect::<Vec<_>>()
                    } else {
                        log::trace!("Filtering extension candidates by read support.");
                        self.filter_candidates(
                            hap,
                            candidates,
                            self.nvar - 1 - suffix_x_count,
                            false,
                        )
                    }
                }) {
                    log::trace!("filtered candidates : {candidates:?}");
                    if let Some(next_hap_len) = candidates.first().map(|x| x.len()) {
                        for idx in 0..next_hap_len {
                            let bases_per_candidate =
                                candidates.iter().map(|x| x[idx]).collect::<BTreeSet<_>>();
                            if bases_per_candidate.len() != 1 {
                                break;
                            }
                            if let Some(base) = bases_per_candidate.first() {
                                identical_bases_right.push(*base);
                            }
                        }
                    }
                }
            } else {
                log::trace!(
                    "Did not find the dnhap match for {suffix_x_count} as key and dnhap {:?}",
                    self.dnhap
                );
            }
        } else {
            log::trace!("No right-side extension candidates for hap {hap}");
        }
        if prefix_x_count >= 2 {
            let next_pos = prefix_x_count;
            let Some(nend) = self.get_nstart_nend().map(|x| x.1) else {
                return Ok(VString::from(hap));
            };
            log::trace!("next pos: {next_pos}, nend = {nend}");
            let hap_block = VString::from(if next_pos < nend as usize {
                let Some(next_next_pos) = self.get_next_pos(next_pos as u32) else {
                    return Ok(VString::from(hap));
                };
                &hap[prefix_x_count..next_next_pos as usize]
            } else {
                &hap[prefix_x_count..]
            });
            let dnhap_is_assigned = self
                .dnhap
                .get(&((prefix_x_count - 1) as u32))
                .map(SegmentationClass::is_assigned)
                == Some(true);
            let passing_edge_support = self
                .edge_info
                .get(prefix_x_count - 1)
                .and_then(|x| x.iter().copied().min())
                .unwrap_or(0);
            let extends_left = dnhap_is_assigned && passing_edge_support >= 2;
            log::trace!("Extending hap {hap:?} left? {extends_left}. dnhap: {dnhap_is_assigned}. Edge support: {passing_edge_support} when querying index {}. Edges: {:?}", (prefix_x_count - 1) as u32, self.edge_info.get(prefix_x_count - 1));
            if extends_left {
                let previous_pos = self
                    .get_previous_pos(next_pos as u32)
                    .ok_or_else(|| anyhow::anyhow!("prev pos missing in extend_left {next_pos}"))?;
                let dbefore = self
                    .subregion_assembly(previous_pos, next_pos as u32, false)?
                    .3;
                if let Some(candidates) = dbefore.get(&hap_block).map(|candidates| {
                    if candidates.len() <= 1 {
                        candidates.iter().cloned().collect::<Vec<_>>()
                    } else {
                        self.filter_candidates(hap, candidates, prefix_x_count, true)
                    }
                }) {
                    if !candidates.is_empty() {
                        if let Some(previous_hap_len) = candidates.first().map(|x| x.len()) {
                            for idx in (0..previous_hap_len).rev() {
                                let bases_per_candidate =
                                    candidates.iter().map(|x| x[idx]).collect::<BTreeSet<_>>();
                                if bases_per_candidate.len() != 1 {
                                    break;
                                }
                                if let Some(base) = bases_per_candidate.first() {
                                    identical_bases_left.push(*base);
                                }
                            }
                            identical_bases_left.reverse();
                        }
                    }
                }
            }
        } else {
            log::trace!("No left-side extension candidates.");
        }
        log::trace!(
            "Identity check for hap {hap}: left={identical_bases_left}, right={identical_bases_right}, prefix_x={prefix_x_count}, suffix_x={suffix_x_count}"
        );
        let mut ret = VString::default();
        ret.append(&mut vec![b'x'; prefix_x_count - identical_bases_left.len()]);
        ret.extend_from_slice(&identical_bases_left);
        ret.extend_from_slice(&hap[prefix_x_count..=suffix_x_count]);
        ret.extend_from_slice(&identical_bases_right);
        ret.resize(self.nvar, b'x');
        log::trace!(
            "Computed extension score for hap {hap}: left={identical_bases_left}, right={identical_bases_right}, final_score={ret}"
        );
        Ok(ret)
    }

    /// Extend all pivot-derived haplotypes using [`Self::extended_hap`].
    ///
    /// Any per-haplotype extension failure is downgraded to the original hap with
    /// a warning, so assembly can continue.
    #[must_use]
    pub(crate) fn extend_pivot_blocks<'a, T>(
        &mut self,
        x: impl IntoIterator<Item = T>,
    ) -> AssembledPaths
    where
        T: std::convert::Into<VStr<'a>> + 'a,
    {
        AssembledPaths::from_seqs(x.into_iter().map(|x| {
            let hap = x.into();
            let res = self.extended_hap(hap).unwrap_or_else(|e| {
                log::warn!("Failed to extend hap {hap}; keeping original hap. Error: {e:?}");
                VString::from(hap)
            });
            log::trace!("Original hap {hap:?} was extended to {res:?}");
            res
        }))
    }

    /// Performs assembly.
    /// Assumes that `init()` has been run successfully.
    pub fn run_asm(&mut self) -> Result<AssemblyResult, DError> {
        let init_haps = self.assemble_haps()?;
        log::debug!("Initially-assembled haps: {init_haps:?}");
        let hap_len = init_haps.iter().next().map_or(0, |x| x.len());
        if hap_len != self.nvar {
            Err(anyhow!(
                "Assembly failed because haplotypes have length {0} but expected {1}",
                hap_len,
                self.nvar
            ))?;
        }
        let final_haps = self.filter_low_support_haps(&init_haps)?;
        log::debug!("Filtered and assembled haps: {final_haps:?}");
        self.display_state("State after filtering");
        let (highest_cn, main_haps_candidates) = self.get_highest_cn(final_haps.iter())?;
        log::debug!(
            "Highest copy number possible: {highest_cn}. {} haplotype set candidates, {main_haps_candidates:?}",
            main_haps_candidates.len()
        );
        let extended_main_hap_candidates = main_haps_candidates
            .iter()
            .map(|x| {
                let extended = self.extend_pivot_blocks(x.iter());
                log::debug!("Extended block {x:?} to {extended:?}");
                extended
            })
            .collect::<Vec<_>>();
        log::debug!("Extended candidates: {extended_main_hap_candidates:?}");
        let main_haps = self.select_main_haps(&final_haps, extended_main_hap_candidates);
        log::debug!("Main haplotypes after candidate selection: {main_haps:?}");

        Ok(AssemblyResult {
            main_haps,
            final_haps,
            highest_cn,
        })
    }

    /// Merge the simple edges for segments that don't need additional assembly.
    pub fn merge_simple_edges_from_segments(
        &mut self,
        segments: &std::collections::BTreeMap<Range<u32>, SegmentationClass>,
    ) -> DResult {
        let num_segments = segments.len();
        for (segment_id, (segment, _height)) in segments
            .iter()
            .enumerate()
            .filter(|(_id, (_segment, height))| **height == SegmentationClass::Two)
        {
            let mut segment_start = segment.start;
            let mut segment_end = segment.end;
            if segment_id > 0 && segment_id != (num_segments - 1) {
                segment_start += 1;
            }
            if segment_end == (self.nvar - 2) as u32 {
                segment_end += 1;
            }
            if segment_end > segment_start {
                self.merge_edges_simple(segment_start, segment_end)?;
            }
        }
        Ok(())
    }

    pub fn merge_complex_regions(
        &mut self,
        segments: std::collections::BTreeMap<Range<u32>, SegmentationClass>,
    ) -> DResult {
        for (segment, _height) in segments.into_iter().filter(|(segment, height)| {
            *height == SegmentationClass::Ten && segment.end >= segment.start
        }) {
            self.merge_edges(segment.start, segment.end + 1)?;
        }
        Ok(())
    }

    pub fn assemble_haps(&mut self) -> Result<AssembledPaths, DError> {
        self.summarize_edges()?;
        let segments = self.segment_graph()?;
        log::trace!(
            "dnhap after segmentation: {:?}. Segments: {segments:?}. Pos edge map: {:?}. Formatted nodes: {:?}",
            self.dnhap,
            self.pos_edge_map,
            self.formatted_nodes()
        );
        assert!(
            !segments.is_empty(),
            "edges: {:?}. inputs {:?}",
            self.pos_edge_map,
            self.reads
        );
        if segments.len() == 1 && matches!(segments.values().next(), Some(SegmentationClass::Two)) {
            if let Some(ivl) = segments.keys().next() {
                return self.path(Range::<u32>::new(ivl.start, ivl.end + 1));
            }
        }

        self.merge_simple_edges_from_segments(&segments)?;
        log::trace!(
            "pos edge map after simple edge merging: {:?}. Formatted nodes {:?}. Formatted edges: {:?}",
            self.pos_edge_map,
            self.formatted_nodes(),
            self.formatted_edges(),
        );
        self.merge_complex_regions(segments)?;
        log::trace!(
            "pos edge map after complex edge merging: {:?}. Formatted nodes {:?}. Formatted pos_edge_map: {:?}",
                self.pos_edge_map, self.formatted_nodes(), self.formatted_edges(),
        );

        let Some((nstart, nend)) = self.get_nstart_nend() else {
            return Ok(AssembledPaths::new());
        };
        if nstart == nend {
            log::trace!("nstart == nend, no need to merge, {nstart}, {nend}");
            let candidates_to_return = AssembledPaths::from_seqs(self.node_iter().map(|x| &x.hap));
            if !candidates_to_return.is_empty()
                && matches!(candidates_to_return.first(), Some(x) if self.nvar == x.len())
            {
                return Ok(candidates_to_return);
            } else {
                return Ok(AssembledPaths::new());
            }
        }
        log::trace!("nstart = {nstart}, nend = {nend}");
        self.merge_where_possible(nstart, nend)?;

        log::trace!(
            "after merge_where_possible: {:?}. Formatted nodes {:?}. Formatted pos_edge_map: {:?}",
            self.pos_edge_map,
            self.formatted_nodes(),
            self.formatted_edges(),
        );

        self.merge_blocks_by_size()?;
        log::trace!(
            "after merge_blocks_by_size: {:?}. Formatted nodes {:?}. Formatted pos_edge_map: {:?}",
            self.pos_edge_map,
            self.formatted_nodes(),
            self.formatted_edges(),
        );

        let Some((nstart, nend)) = self.get_nstart_nend() else {
            return Ok(AssembledPaths::new());
        };
        Ok(if nstart == nend {
            log::debug!(
            "final subregion_assembly but nstart == nend: {:?}. Formatted nodes {:?}. Formatted pos_edge_map: {:?}. haps: {:?}",
                self.pos_edge_map, self.formatted_nodes(), self.formatted_edges(), self.node_iter().map(|x| &x.hap[..]).collect::<Vec<_>>());
            AssembledPaths::from_seqs(self.node_iter().map(|x| x.hap.clone()))
        } else {
            log::debug!(
            "before final subregion_assembly: {:?}. Formatted nodes {:?}. Formatted pos_edge_map: {:?}.",
                self.pos_edge_map, self.formatted_nodes(), self.formatted_edges());
            let (_success, haps, _x, _y) = self.subregion_assembly(nstart, nend, true)?;
            log::debug!(
            "after final subregion_assembly: {:?}. Formatted nodes {:?}. Formatted pos_edge_map: {:?}. haps: {haps:?}",
                self.pos_edge_map, self.formatted_nodes(), self.formatted_edges(),
        );
            AssembledPaths::from_seqs(haps)
        })
    }
}
