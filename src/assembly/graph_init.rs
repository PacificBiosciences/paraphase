use crate::assembly::node_datum::NodeDatum;
use crate::assembly::variant_graph::{
    Graph, NodeMapType, ReadHapPair, ReadHapSupport, SegmentationClass, WeightedEdge,
};
use crate::io::json::{ReadAlignmentId, ReadFingerprintMap};
use crate::toolkit::hapcmp::HapCompare;
use crate::toolkit::util::{DResult, HashMap, HashSet};

use petgraph::stable_graph::NodeIndex;
use petgraph::visit::{EdgeRef, IntoEdgeReferences};

use std::collections::BTreeSet;
use std::hash::BuildHasherDefault;
use vstr::{VStr, VString};
impl Graph {
    /// Checks if the node is made entirely of non-1 or 2.
    /// This means it's not ref or non-ref, so it's a deletion or other special operation.
    #[inline]
    pub(crate) fn node_is_sv(node_id: &NodeDatum) -> bool {
        node_id.hap.iter().all(|x| !matches!(x, b'1' | b'2'))
    }

    /// Compute a map from node ids to node references.
    #[must_use]
    pub(crate) fn get_id_lookup(&self) -> HashMap<u32, &NodeDatum> {
        self.node_id_map
            .iter()
            .map(|(x, y)| (*y, x))
            .collect::<HashMap<u32, &NodeDatum>>()
    }

    /// Remove an edge from the graph using a `WeightedEdge`, which is numeric ids.
    pub(crate) fn rm_edge(&mut self, edge: WeightedEdge) -> DResult {
        let index = self
            .graph
            .find_edge(
                NodeIndex::new(edge.0 as usize),
                NodeIndex::new(edge.1 as usize),
            )
            .ok_or_else(|| {
                anyhow::anyhow!("Cannot remove edge that does not exist. [graph.find_edge()]")
            })?;
        self.graph.remove_edge(index).ok_or_else(|| {
            anyhow::anyhow!("Cannot remove edge that does not exist [graph.rm_edge(index)]")
        })?;
        Ok(())
    }

    /// Find reads that match given haplotypes and vice versa.
    pub fn match_reads_and_haps<'a, 'b, Iter, T>(
        reads: &'a ReadFingerprintMap,
        haps: Iter,
        min_match: Option<i32>,
    ) -> Result<ReadHapSupport<'a>, Box<dyn std::error::Error>>
    where
        vstr::VStr<'b>: From<T>,
        Iter: IntoIterator<Item = T>,
    {
        let haps = haps
            .into_iter()
            .map(VStr::from)
            .map(VString::from)
            .collect::<Vec<_>>();
        let min_match = min_match.unwrap_or(1);
        let mut support_reads_nonunique: HashMap<u32, Vec<ReadHapPair>> =
            HashMap::with_capacity_and_hasher(haps.len(), BuildHasherDefault::default());
        let mut support_reads: HashMap<u32, Vec<ReadHapPair>> =
            HashMap::with_capacity_and_hasher(haps.len(), BuildHasherDefault::default());
        for idx in 0..haps.len() {
            support_reads_nonunique.insert(idx as u32, vec![]);
            support_reads.insert(idx as u32, vec![]);
        }
        let mut supporting_haps_per_read: HashMap<ReadHapPair, Vec<u32>> =
            HashMap::with_capacity_and_hasher(reads.len() / 2, Default::default());
        for (fprint, info) in reads.iter().map(|(k, v)| {
            (
                v,
                ReadHapPair {
                    name: &k.read_name[..],
                    path: v,
                },
            )
        }) {
            let mut matching_hapids: Vec<u32> = Vec::new();
            for (hap_id, hap_to_extend) in haps.iter().map(|x| -> VStr<'_> { x.into() }).enumerate()
            {
                let result = HapCompare::from_haps(hap_to_extend, fprint.vstr())?;
                if result.mismatches == 0 && result.matches >= min_match && result.extra >= 0 {
                    let hap_id = hap_id as u32;
                    matching_hapids.push(hap_id);
                    support_reads_nonunique
                        .entry(hap_id)
                        .or_default()
                        .push(info);
                }
            }
            assert!(
                supporting_haps_per_read
                    .insert(info, matching_hapids.clone())
                    .is_none(),
                "Duplicate read id in map... is this possible?"
            );
            if matching_hapids.len() == 1 {
                support_reads
                    .entry(matching_hapids[0])
                    .or_default()
                    .push(info);
                continue;
            }

            for lh_hap_id in matching_hapids.iter().copied() {
                let has_overlap = matching_hapids.iter().copied().any(|rh_hap_id| {
                    rh_hap_id != lh_hap_id && {
                        let compare = HapCompare::from_haps(
                            &haps[lh_hap_id as usize],
                            &haps[rh_hap_id as usize],
                        );
                        let Ok(compare) = compare else {
                            return false;
                        };
                        std::cmp::max(compare.matches, compare.mismatches) > 0
                    }
                });
                if !has_overlap {
                    support_reads.entry(lh_hap_id).or_default().push(info);
                }
            }
        }
        Ok(ReadHapSupport {
            support_reads,
            support_reads_nonunique,
            supporting_haps_per_read,
        })
    }

    /// Add the node to all structures - the graph structure, `node_id_map`, and related.
    /// This is the only add node function you should be using.
    /// Returns (`node_id`, `was_new`) when adding.
    /// If it was not new, the original id is returned.
    #[must_use]
    pub fn fully_add_node(&mut self, datum: NodeDatum) -> (u32, bool) {
        let mut total_node_count = self.total_node_count;
        let (idx, was_new) = self.fetch_add_node(&datum, Some(&mut total_node_count));
        self.total_node_count = total_node_count;
        log::trace!("Registering node {datum:?}: id={idx}, was_new={was_new}");
        if was_new {
            log::trace!(
                "Inserted new node {datum:?} with id={idx}. graph_nodes={}, node_id_map_entries={}",
                self.graph.node_count(),
                self.node_id_map.len()
            );
            let graph_idx = self.add_node(datum).index() as u32;
            assert_eq!(graph_idx, idx);
        }
        (idx, was_new)
    }

    /// Adds the node to the map if not present. Also manages the number of nodes, outside of the structure for mutability.
    /// This function is expected to only be used by maintainers.
    #[must_use]
    pub fn fetch_add_node(
        &mut self,
        datum: &NodeDatum,
        num_nodes: Option<&mut usize>,
    ) -> (u32, bool) {
        if let Some(num_nodes) = num_nodes {
            Self::id_node_by_map(&mut self.node_id_map, datum, num_nodes)
        } else {
            let mut num_nodes = self.total_node_count;
            let res = Self::id_node_by_map(&mut self.node_id_map, datum, &mut num_nodes);
            self.total_node_count = num_nodes;
            res
        }
    }

    /// Id a node if present, but do not add.
    #[must_use]
    pub(crate) fn id_node_const(&self, datum: &NodeDatum) -> Option<u32> {
        self.node_id_map.get(datum).copied()
    }

    /// Id a node using a non-member map.
    /// This lets us have another mutable borrow on the Graph object.
    #[must_use]
    pub(crate) fn id_node_by_map(
        map: &mut NodeMapType,
        datum: &NodeDatum,
        num_nodes: &mut usize,
    ) -> (u32, bool) {
        match map.get(datum) {
            Some(x) => {
                // log::trace!("Found node: {datum:?} at id {x}");
                (*x, false)
            }
            None => {
                let id = *num_nodes as u32;
                *num_nodes += 1;
                log::trace!(
                    "Node {datum:?} was not present; assigning new id={id}. total_nodes now={}",
                    *num_nodes
                );
                map.insert(datum.clone(), id);
                (id, true)
            }
        }
    }

    /// Convert a `read_aln_id` into a numeric id.
    #[must_use]
    pub fn read_aln_id(&self, name: &ReadAlignmentId) -> Option<u32> {
        self.read_aln_id_map.get(name).copied()
    }

    #[inline]
    /// Return whether a byte is a concrete allele token (not `'x'` unknown).
    pub(crate) fn is_valid(x: &u8) -> bool {
        *x != b'x'
    }

    /// Count non-`'x'` tokens in a haplotype slice.
    #[must_use]
    fn count_valid(hap: &[u8]) -> usize {
        hap.iter().copied().filter(Self::is_valid).count()
    }

    /// Initialize graph nodes/edges from read haplotype fingerprints.
    ///
    /// This populates `graph`, `node_id_map`, and edge weights based on adjacent
    /// valid symbols in each read path.
    pub fn init_edges(&mut self) {
        let mut node_id_map = NodeMapType::new();
        std::mem::swap(&mut node_id_map, &mut self.node_id_map);
        log::trace!("Initializing graph edges from {} reads", self.reads.len());
        let mut total_node_count = self.total_node_count;
        for (read_aln_id, hap) in self.reads.iter().filter(|(_name, hap)| {
            let num_valid = Self::count_valid(&hap[..]);
            assert_eq!(hap.len(), self.nvar);
            num_valid > 0
        }) {
            let read_name = &read_aln_id.read_name[..];
            let Some(read_id) = self.read_aln_id(read_aln_id) else {
                continue;
            };
            log::trace!("Processing read {read_name} (id={read_id}) with hap={hap}");
            for (window_index, window) in hap
                .windows(2)
                .enumerate()
                .filter(|x| x.1.iter().all(Self::is_valid))
            {
                let (Some(from), Some(to)) = (window.first().copied(), window.get(1).copied())
                else {
                    continue;
                };
                let from_pos = NodeDatum::from((from, window_index));
                let to_pos = NodeDatum::from((to, window_index + 1));
                let (from_id, from_was_new) =
                    Self::id_node_by_map(&mut node_id_map, &from_pos, &mut total_node_count);
                if from_was_new {
                    let graph_id = self.graph.add_node(from_pos).index() as u32;
                    assert_eq!(from_id, graph_id);
                }
                let (to_id, to_was_new) =
                    Self::id_node_by_map(&mut node_id_map, &to_pos, &mut total_node_count);
                if to_was_new {
                    let graph_id = self.graph.add_node(to_pos);
                    assert_eq!(to_id as usize, graph_id.index());
                }
                let edge_idx = if let Some(e) = self.graph.find_edge(from_id.into(), to_id.into()) {
                    if let Some(weight) = self.graph.edge_weight_mut(e) {
                        weight.push(read_id);
                    }
                    e
                } else {
                    self.graph
                        .add_edge(from_id.into(), to_id.into(), vec![read_id])
                };
                log::trace!(
                    "Assigned edge {edge_idx:?} from {}@{from_id} to {}@{to_id}",
                    from as char,
                    to as char
                );
            }
        }
        std::mem::swap(&mut node_id_map, &mut self.node_id_map);
        std::mem::swap(&mut total_node_count, &mut self.total_node_count);
        let id_to_node = self.get_id_lookup();
        let mut pos_edge_map: std::collections::BTreeMap<usize, Vec<(u32, u32, u32)>> =
            std::collections::BTreeMap::new();
        for i in 0..(self.nvar - 1) {
            pos_edge_map.insert(i, vec![]);
        }
        let mut edge_info: Vec<Vec<u32>> = vec![vec![]; self.nvar - 1];
        for edge in self.graph.edge_references() {
            let source = edge.source();
            let target = edge.target();
            let weight = edge.weight().len() as u32;
            let v = (source.index() as u32, target.index() as u32, weight);
            let position_index = id_to_node
                .get(&(source.index() as u32))
                .map(|x| x.pos as usize);
            let Some(position_index) = position_index else {
                continue;
            };
            pos_edge_map.entry(position_index).or_default().push(v);
            edge_info[position_index].push(weight);
        }
        self.pos_edge_map = pos_edge_map;
        self.edge_info = edge_info;
    }

    #[must_use]
    /// Build per-position edge-support vectors from `pos_edge_map`.
    ///
    /// Each entry corresponds to one genomic transition position and contains
    /// support counts for all edges at that transition.
    pub fn make_edge_info(&self) -> Vec<Vec<u32>> {
        let id_to_node = self
            .get_id_lookup()
            .into_iter()
            .map(|(x, y)| (x, y.to_owned()))
            .collect::<HashMap<_, _>>();
        let mut edge_info: Vec<Vec<u32>> = vec![vec![]; self.nvar - 1];
        for edge in self.graph.edge_references() {
            let source = edge.source().index() as u32;
            let weight = edge.weight().len() as u32;
            let Some(position_index) = id_to_node.get(&source).map(|x| x.pos as usize) else {
                continue;
            };
            edge_info[position_index].push(weight);
        }
        edge_info
    }

    /// Check the number of different edges starting at each position.
    /// Summarize local edge complexity into `dnhap` segmentation labels.
    ///
    /// Uses max/min edge support at each position plus edge-count heuristics to
    /// classify regions as simple (`Two`) or complex (`Ten`/`Multiple` variants).
    pub fn summarize_edges(&mut self) -> DResult {
        let mut site1 = HashSet::default();
        let mut site2 = HashSet::default();
        let mut expected_site1 = HashSet::default();
        let mut expected_site2 = HashSet::default();
        for (idx, edges) in &self.pos_edge_map {
            let idx = *idx as u32;
            let nhap: SegmentationClass = if edges.iter().map(|x| x.2).max().unwrap_or(1) > 2
                && edges.iter().map(|x| x.2).min().unwrap_or(1) >= 2
            {
                let mut dlinks: HashMap<u32, BTreeSet<u32>> = HashMap::default();
                for x in [
                    &mut site1,
                    &mut site2,
                    &mut expected_site1,
                    &mut expected_site2,
                ] {
                    x.clear();
                }
                for node_index in self.graph.node_indices() {
                    let node: &NodeDatum = &self.graph[node_index];
                    let base = &node.hap[..];
                    let pos = node.pos;
                    match (pos as i32) - (idx as i32) {
                        0 => {
                            expected_site1.insert(base);
                        }
                        1 => {
                            expected_site2.insert(base);
                        }
                        _ => {}
                    }
                }
                for (from_id, to_id, _count) in edges {
                    let from = &self.graph[NodeIndex::new(*from_id as usize)];
                    let to = &self.graph[NodeIndex::new(*to_id as usize)];
                    site1.insert(&from.hap[..]);
                    site2.insert(&to.hap[..]);
                    dlinks.entry(*from_id).or_default().insert(*to_id);
                }
                if (&site1, &site2) == (&expected_site1, &expected_site2) {
                    if dlinks.values().any(|x| x.len() > 1) {
                        SegmentationClass::Ten
                    } else {
                        SegmentationClass::Two
                    }
                } else {
                    SegmentationClass::Unassigned
                }
            } else {
                SegmentationClass::Unassigned
            };
            self.dnhap.insert(idx, nhap);
        }
        if self.settings.pivot_index > 0 {
            if let Some(height) = self
                .dnhap
                .get_mut(&((self.settings.pivot_index - 1) as u32))
            {
                if *height == SegmentationClass::Two {
                    *height = SegmentationClass::Ten;
                }
            }
        }
        Ok(())
    }

    /// Computes the minimum edge support for given read counts.
    /// Higher coverage samples need slightly higher thresholds.
    #[must_use]
    /// Compute the dynamic minimum edge-support threshold for pruning.
    ///
    /// Threshold is scaled by read depth and constrained by graph settings.
    pub fn edge_support_threshold<T: TryInto<usize> + Copy>(num_reads: &[T]) -> u32 {
        let converted = num_reads
            .iter()
            .filter_map(|x| match (*x).try_into() {
                Ok(v) => Some(v),
                Err(_) => {
                    log::warn!("edge_support_threshold received a non-convertible read count; skipping it.");
                    None
                }
            })
            .collect::<Vec<_>>();
        let min_count = num_reads
            .iter()
            .filter_map(|x| (*x).try_into().ok())
            .min()
            .unwrap_or(0);
        let total_min = converted
            .iter()
            .copied()
            .filter(|x| *x != min_count)
            .min()
            .unwrap_or(0);
        if total_min > 10 {
            (total_min / 25).clamp(2, 5) as u32
        } else {
            u32::from(total_min > 5)
        }
    }

    /// Remove edges which have low support.
    /// Remove low-support edges and update per-position edge caches.
    ///
    /// This preserves high-confidence transitions while dropping weak links that
    /// are unlikely to contribute to consistent haplotype paths.
    pub fn prune_edges(&mut self) -> DResult {
        log::debug!("Pruning low-support graph edges");
        let id_to_node = self
            .get_id_lookup()
            .into_iter()
            .map(|(x, y)| (x, y.to_owned()))
            .collect::<HashMap<_, _>>();
        assert!(!self.pos_edge_map.is_empty());
        let mut tmp_pos_edge = std::collections::BTreeMap::<usize, Vec<WeightedEdge>>::default();
        std::mem::swap(&mut tmp_pos_edge, &mut self.pos_edge_map);
        let mut edges_to_remove: std::collections::BTreeMap<usize, Vec<(u32, u32, u32)>> =
            std::collections::BTreeMap::new();
        for (pos, edges) in &tmp_pos_edge {
            let num_edges = edges.iter().map(|x| x.2).collect::<Vec<_>>();
            let has_sv = edges
                .iter()
                .flat_map(|x| [x.0, x.1].into_iter())
                .any(|x| Self::node_is_sv(&id_to_node[&x]));

            let threshold = Self::edge_support_threshold(&num_edges[..]);
            log::trace!(
                "Edge-pruning context at pos {pos}: has_sv={has_sv}, threshold={threshold}, edge_supports={num_edges:?}"
            );
            for edge in edges {
                let (from, to, count) = edge;
                let outgoing_from = self
                    .graph
                    .edges_directed(
                        NodeIndex::new(*from as usize),
                        petgraph::Direction::Outgoing,
                    )
                    .count();
                let incoming_to = self
                    .graph
                    .edges_directed(NodeIndex::new(*to as usize), petgraph::Direction::Incoming)
                    .count();
                let nonunique_edges = (outgoing_from > 1) && (incoming_to > 1);
                let fails = nonunique_edges && ((*count <= threshold) || (has_sv && *count == 1));
                if fails {
                    edges_to_remove.entry(*pos).or_default().push(*edge);
                }
            }
        }
        std::mem::swap(&mut tmp_pos_edge, &mut self.pos_edge_map);
        log::trace!("Edges selected for removal: {edges_to_remove:?}");
        for (pos, edges) in &edges_to_remove {
            for edge in edges {
                self.rm_edge(*edge)?;
                if let Some(val) = self.pos_edge_map.get_mut(pos) {
                    val.retain(|x| x != edge);
                }
            }
        }

        let mut read_ids_to_delete = Vec::new();
        for (read_name, read_seq) in &mut self.reads {
            let mut new_seq = read_seq.clone();
            let mut num_found = 0usize;
            let mut total = 0usize;
            for read_idx in (0..(read_seq.len() - 1))
                .filter(|x| read_seq[*x..(*x + 2)].iter().all(Self::is_valid))
            {
                total += 1;
                let from = read_seq[read_idx];
                let to = read_seq[read_idx + 1];
                let from_pos = NodeDatum::from((from, read_idx));
                let to_pos = NodeDatum::from((to, read_idx + 1));
                let from_pos = self.node_id_map.get(&from_pos).copied().unwrap_or(u32::MAX);
                let to_pos = self.node_id_map.get(&to_pos).copied().unwrap_or(u32::MAX);
                let found = self.pos_edge_map[&read_idx]
                    .iter()
                    .any(|x| (x.0, x.1) == (from_pos, to_pos));
                if found {
                    num_found += 1;
                } else {
                    log::trace!(
                        "Removed edge support between {from}/{from_pos} and {to}/{to_pos} at site {read_idx}; previous read segment={}",
                        VStr::from(&new_seq[read_idx..read_idx + 2])
                    );
                    new_seq[read_idx] = b'x';
                    new_seq[read_idx + 1] = b'x';
                }
            }
            log::trace!(
                "Read {read_name}: retained {num_found}/{total} traversed edges after pruning"
            );
            if bytecount::count(&new_seq, b'x') <= self.nvar - 2 {
                *read_seq = new_seq;
            } else {
                read_ids_to_delete.push(read_name.clone());
            }
        }

        for read_aln_id in read_ids_to_delete {
            assert!(
                self.reads.remove(&read_aln_id).is_some(),
                "Failed to remove read_aln_id: {read_aln_id}"
            );
        }

        self.edge_info = self.make_edge_info();
        log::trace!(
            "Pruned edges; edge info after construction: {:?}",
            self.edge_info
        );
        self.display_state("Pruned edges.");
        Ok(())
    }

    /// Check for possible haps made from missing for later extension.
    /// Compute read-to-haplotype matching support for all current graph paths.
    ///
    /// Stores both unique and non-unique support in `read_support`.
    fn match_all_reads_and_haps(&mut self) -> DResult {
        const ALPHABET: &[u8] = b"01234xACGT";
        for (i, pos, node) in itertools::iproduct!(ALPHABET.iter(), 0..self.nvar)
            .map(|(i, pos)| (i, pos, NodeDatum::from((*i, pos))))
        {
            if self.node_id_map.contains_key(&node) {
                continue;
            }
            let mut hap = VString::from(vec![b'x'; pos]);
            hap.reserve_exact(self.nvar);
            hap.push(*i);
            hap.resize(self.nvar, b'x');
            let read_support = Self::match_reads_and_haps(
                &self.reads_original,
                &[hap],
                Some(self.settings.min_match_count),
            )?;
            if read_support
                .support_reads
                .get(&0)
                .map_or(0, std::vec::Vec::len)
                > 1
            {
                let (node_id, node_was_new) =
                    Self::id_node_by_map(&mut self.node_id_map, &node, &mut self.total_node_count);
                debug_assert!(node_was_new, "Node should have been new at construction.");
                let graph_id = self.graph.add_node(node);
                debug_assert_eq!(graph_id.index() as u32, node_id);
            }
        }

        Ok(())
    }

    /// For debug purposes, print nodes/edges.
    /// Emit a compact debug snapshot of graph state.
    ///
    /// Intended for trace-level diagnostics during graph construction/merging.
    pub(crate) fn display_state(&self, msg: impl Into<String>) {
        let msg = msg.into();
        log::debug!(
            "{msg}. Nodes: {:?}. Edges: {:?}",
            self.formatted_nodes(),
            self.formatted_edges()
        );
    }

    /// Initialization:
    /// Set up graph, prune edges, and match reads and haplotypes.
    /// Initialize graph internals end-to-end before assembly.
    ///
    /// Runs edge construction, summarization, support computation, and pruning.
    pub fn init(&mut self) -> DResult {
        self.init_edges();
        self.display_state("[init] Edges initialized.");

        self.prune_edges()?;
        self.display_state("[init] Edges pruned.");

        self.match_all_reads_and_haps()?;
        self.display_state("[init] reads matched to haplotypes.");

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::Graph;

    #[test]
    fn threshold_count_ok() {
        assert_eq!(
            Graph::edge_support_threshold(&[1usize, 10usize, 10usize]),
            1
        );
        assert_eq!(Graph::edge_support_threshold(&[24, 9]), 2);
    }
}
