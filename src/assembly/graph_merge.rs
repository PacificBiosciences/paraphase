use crate::assembly::node_datum::NodeDatum;
use crate::assembly::variant_graph::{Graph, MatchingHapSet, MergePosResult, WeightedEdge};
use crate::io::json::ReadAlignmentId;
use crate::toolkit::hapcmp::HapCompare;
use crate::toolkit::range::Range;
use crate::toolkit::util::{DError, DResult, HashMap};

use petgraph::stable_graph::NodeIndex;
use vstr::{VStr, VString};

use itertools::iproduct;

use std::collections::{BTreeMap, BTreeSet};

impl Graph {
    /// Find the best matches for nodes given a query fingerprint.
    pub fn get_matching_hap<'a, 'b, C, IntoVstr: std::convert::Into<VStr<'a>>>(
        nodes: &'a BTreeMap<u32, C>,
        pos1: u32,
        hap: VStr<'b>,
    ) -> Result<(VStr<'b>, Option<MatchingHapSet<'a>>), DError>
    where
        &'a C: std::iter::IntoIterator<Item = IntoVstr>,
    {
        let candidates = nodes
            .get(&pos1)
            .ok_or_else(|| anyhow::anyhow!("Missing position"))?
            .into_iter()
            .map(std::convert::Into::into)
            .collect::<Vec<_>>();
        let pos1 = pos1 as usize;
        let pos2 = candidates[0].len();
        if pos1 + pos2 > hap.len() {
            return Err(anyhow::anyhow!("{pos1} and {pos2} for hap {hap} is out of bounds").into());
        }
        let hap1 = VStr::from(&hap[pos1..(pos1 + pos2)]);
        log::trace!(
            "Matching query hap segment {hap1} against candidates {candidates:?} (pos1={pos1}, segment_len={pos2})"
        );

        let mut matches: BTreeMap<_, Vec<_>> = BTreeMap::new();
        for cand in candidates {
            let cmp = HapCompare::from_haps(&cand[..], hap1)?;
            if cmp.mismatches > 0 || cmp.matches <= 0 {
                continue;
            }
            matches.entry(cmp.matches).or_default().push(cand);
        }

        let best_match = matches
            .into_iter()
            .max_by_key(|x| x.0)
            .map(|(key, value)| (key, value.clone()));
        Ok((hap1, best_match))
    }

    /// Progressively merge adjacent blocks in `[nstart, nend]` when subregion
    /// assembly can extend haplotype spans.
    pub fn merge_where_possible(&mut self, mut nstart: u32, nend: u32) -> DResult {
        log::trace!(
            "merge_where_possible start: pos_edge_map={:?}, nodes={:?}, formatted_edges={:?}",
            self.pos_edge_map,
            self.formatted_nodes(),
            self.formatted_edges(),
        );
        let mut inum = 1usize;
        while nstart < nend - 1 {
            let istr = format!("[merge_where_possible@{inum}");
            log::trace!(
                "{istr} Pre-assembly graph snapshot: edges={:?}, nodes={:?}",
                self.formatted_edges(),
                self.formatted_nodes()
            );
            let (success, final_haps, _dnext, _dbefore) =
                self.subregion_assembly(nstart, nend, false)?;
            log::trace!(
                "{istr} Subregion assembly result for {nstart}..={nend}: success={success}, final_haps={final_haps:?}"
            );
            let current_haps = self.get_hap_by_pos(nstart).collect::<Vec<_>>();
            log::trace!(
                "{istr} Current position haps for {nstart}..={nend}: {current_haps:?}; dnext={_dnext:?}, dbefore={_dbefore:?}"
            );
            let final_len = final_haps[0].len();
            let pos = nstart + (final_len as u32);
            let Some(block_ending_pos) = self.get_previous_pos(pos) else {
                break;
            };
            let current_len = current_haps[0].len();
            if current_len != final_len {
                log::trace!(
                    "{istr} Rebuilding edges between {nstart} and {block_ending_pos} after hap-length change."
                );
                self.rm_add_edges(nstart, block_ending_pos, final_haps.into_iter())?;
            } else {
                log::trace!(
                    "{istr} Haplotype lengths unchanged (current={current_len}, final={final_len}); skipping edge rebuild."
                );
            }
            nstart = pos;
            if success {
                log::trace!("{istr} Completed region assembly in {inum} iterations.");
                break;
            }
            inum += 1;
        }
        Ok(())
    }

    /// Given a pair of (position, hap) pairs in a dictionary, merge the relevant positions.
    /// Returns `sub_hap_support`, dnext, dbefore.
    pub fn merge_two_pos<'a>(
        &'a self,
        new_nodes: &'a BTreeMap<u32, BTreeSet<VStr<'_>>>,
    ) -> Result<MergePosResult<'a>, DError> {
        let mut sub_hap_support: BTreeMap<VString, Vec<&ReadAlignmentId>> = BTreeMap::new();
        let mut dnext: BTreeMap<VStr, BTreeSet<VStr>> = BTreeMap::new();
        let mut dbefore: BTreeMap<VStr, BTreeSet<VStr>> = BTreeMap::new();
        log::trace!("Merging two positions with node sets: {new_nodes:?}");
        for (read, fingerprint) in &self.reads {
            let hits = new_nodes
                .keys()
                .map(|x| {
                    let (hap, matches) =
                        Self::get_matching_hap(new_nodes, *x, fingerprint.vstr()).ok()?;
                    matches.map(|matches| (hap, matches))
                })
                .collect::<Vec<_>>();
            if hits.iter().any(Option::is_none) {
                continue;
            }
            let hits = hits
                .into_iter()
                .map(std::option::Option::unwrap)
                .collect::<Vec<_>>();

            let mut num_matches = hits.iter().map(|(_hap, matches)| matches.1.len());
            if let (Some(1), Some(1)) = (num_matches.next(), num_matches.next()) {
                let mut joined = VString::from(&hits[0].1 .1[0][..]);
                joined.extend_from_slice(&hits[1].1 .1[0][..]);
                sub_hap_support.entry(joined).or_default().push(read);
            }
            log::trace!(
                "Per-read match candidates: left={:?}, right={:?}",
                hits[0].1 .1,
                hits[1].1 .1
            );
            for (a, b) in iproduct!(hits[0].1 .1.iter(), hits[1].1 .1.iter()) {
                dnext.entry(*a).or_default().insert(*b);
                dbefore.entry(*b).or_default().insert(*a);
            }
        }
        Ok((sub_hap_support, dnext, dbefore))
    }

    /// Attempt to merge a specific region by assembly, then rewrite local graph
    /// edges/nodes if the merge succeeds.
    pub fn merge_edges(&mut self, first: u32, last: u32) -> DResult {
        let (success, new_haps, _, _) =
            self.subregion_assembly(first, last, /* allow_w= */ false)?;
        if success {
            log::trace!(
                "Subregion assembly succeeded for {first}..={last}; rebuilding edges with new_haps={new_haps:?}"
            );
            self.rm_add_edges(first, last, new_haps.into_iter())?;
        }
        Ok(())
    }

    /// Handle merging in simple portions of the graph.
    pub(crate) fn merge_edges_simple(&mut self, start: u32, end: u32) -> DResult {
        log::trace!(
            "Merging simple edge block from {start} to {end}. pos_edge_map before merge: {:?}",
            self.pos_edge_map
        );
        let path = self.path(Range::<u32>::new(start, end))?;
        self.rm_add_edges(start, end, path.iter().map(|x| x.vstr()))?;
        log::trace!(
            "Completed simple edge merge from {start} to {end}. pos_edge_map after merge: {:?}",
            self.pos_edge_map
        );
        Ok(())
    }

    /// Mark a node as deleted and keep id-map consistency across graph caches.
    ///
    /// The old datum is replaced with its deleted form in `node_id_map` while
    /// preserving the numeric node id.
    pub(crate) fn fully_remove_node(&mut self, node: NodeIndex) {
        log::trace!("Marking node {} as deleted", node.index());
        self.node_positions.remove(&(node.index() as u32));
        let mut weight = if let Some(weight) = self.node_weight_mut(node) {
            let weight_copy = weight.clone();
            assert!(!weight.is_del());
            weight.set_is_del();
            weight_copy
        } else {
            log::warn!("Node weight was unexpectedly missing while deleting node {node:?}");
            return;
        };
        let Some(id) = self.node_id_map.remove(&weight) else {
            log::warn!("Node id was unexpectedly missing from node_id_map for node {node:?}");
            return;
        };
        weight.set_is_del();
        self.node_id_map.insert(weight, id);
    }

    /// Rebuild the local graph for `[pos1, pos2]` using assembled path strings.
    ///
    /// This removes overlapping intermediate edges/nodes, inserts merged nodes,
    /// reconnects flanking edges, and refreshes edge support metadata.
    pub(crate) fn rm_add_edges<T: std::convert::Into<VString>>(
        &mut self,
        pos1: u32,
        pos2: u32,
        path: impl IntoIterator<Item = T>,
    ) -> DResult {
        log::trace!(
            "pos_edge_map before edge rebuild for range {pos1}..={pos2}: {:?}",
            self.pos_edge_map
        );
        let old_edges = self.pos_edge_map.clone();
        let prev = self.get_previous_pos(pos1);
        let start = prev.unwrap_or(0);
        let pos_edge_map_len = self.pos_edge_map.len();
        for idx in (start as usize..=pos2 as usize).filter(|x| *x < pos_edge_map_len) {
            if self.pos_edge_map.contains_key(&idx) {
                let edges_to_remove = self.pos_edge_map[&idx].clone();
                for edge in edges_to_remove {
                    self.rm_edge(edge)?;
                }
                self.pos_edge_map.retain(|k, _| *k != idx);
            }
        }
        if let Some(pos) = prev {
            self.pos_edge_map.retain(|k, _| *k != pos as usize);
        }
        let contain_range = pos1..=pos2;
        let nodes_to_remove = self
            .node_iter()
            .filter(|x| contain_range.contains(&x.pos))
            .filter_map(|x| self.id_node_const(x))
            .map(NodeIndex::from)
            .collect::<Vec<_>>();
        for node in nodes_to_remove {
            self.fully_remove_node(node);
        }
        let mut id_lookup = self
            .get_id_lookup()
            .into_iter()
            .map(|(x, y)| (x, y.clone()))
            .collect::<HashMap<_, _>>();

        let path = path
            .into_iter()
            .map(std::convert::Into::into)
            .collect::<Vec<_>>();
        let npaths = path.len();
        for (i, hap) in path.iter().enumerate() {
            let i = i as u32;
            let node = NodeDatum::new(hap.clone(), pos1);
            let (node_id, was_new) = self.fully_add_node(node.clone());
            log::trace!("Inserted rebuilt node with id={node_id}; was_new={was_new}");
            id_lookup.insert(
                node_id,
                self.node_weight(node_id.into())
                    .cloned()
                    .unwrap_or(node.clone()),
            );
            assert!(was_new);
            self.node_positions.insert(
                node_id,
                (
                    pos1 * 1000,
                    if npaths > 2 { i * 5 + 10 } else { (i + 1) * 10 },
                ),
            );
            let maybe_add_edge = |pos_edge_map: &mut BTreeMap<usize, Vec<WeightedEdge>>,
                                  pos: usize,
                                  node1: u32,
                                  node2: u32| {
                pos_edge_map.entry(pos).or_default();
                if let Some(pos_edges) = pos_edge_map.get_mut(&pos) {
                    if !pos_edges.iter().any(|x| (x.0, x.1) == (node1, node2)) {
                        pos_edges.push((node1, node2, 0));
                    }
                }
            };
            if let Some(old) = old_edges.get(&(pos2 as usize)) {
                for (from, to, _) in old {
                    let Some(from_node) = id_lookup.get(from) else {
                        continue;
                    };
                    let Some(last_base) = hap.last() else {
                        continue;
                    };
                    if *last_base != from_node.hap[0] {
                        continue;
                    }
                    let to = *to;
                    self.add_edge(node_id.into(), to.into(), vec![]);
                    maybe_add_edge(&mut self.pos_edge_map, pos1 as usize, node_id, to);
                }
            }
            if let Some((prev_pos, edges)) = prev
                .map(|x| x as usize)
                .and_then(|x| old_edges.get(&x).map(|edges| (x, edges)))
            {
                edges.iter().copied().for_each(|(node1, node2, _count)| {
                    let Some(node2_node) = id_lookup.get(&node2) else {
                        return;
                    };
                    if hap[0] == node2_node[0] {
                        self.add_edge(node1.into(), node_id.into(), vec![]);
                        maybe_add_edge(&mut self.pos_edge_map, prev_pos, node1, node_id);
                    }
                });
            }
        }
        Ok(())
    }
}
