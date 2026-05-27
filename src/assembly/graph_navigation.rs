use crate::assembly::assembly_result::AssembledPaths;
use crate::assembly::node_datum::NodeDatum;
use crate::assembly::variant_graph::Graph;
use crate::toolkit::range::Range;
use crate::toolkit::util::{all_simple_paths, DError};
use vstr::{VStr, VString};

use itertools::intersperse;
use std::collections::BTreeSet;

impl Graph {
    /// Iterate over all nodes.
    /// Unlike the node iteration provided by petgraph, this is aware of deleted nodes.
    pub fn node_iter(&self) -> impl Iterator<Item = &NodeDatum> {
        self.node_id_map.keys().filter(|x| !x.is_del())
    }

    /// Generate potential paths between nodes at positions.
    pub fn path(&self, ivl: Range<u32>) -> Result<AssembledPaths, DError> {
        let mut ret = AssembledPaths::new();
        let sources = self.pos_edge_map[&(ivl.start as usize)]
            .iter()
            .filter(|x| self.find_edge(x.0.into(), x.1.into()).is_some())
            .map(|x| x.0)
            .collect::<BTreeSet<_>>();
        let targets = self.pos_edge_map[&((ivl.end - 1) as usize)]
            .iter()
            .filter(|x| self.find_edge(x.0.into(), x.1.into()).is_some())
            .map(|x| x.1)
            .collect::<BTreeSet<_>>();
        for (source, target) in
            itertools::iproduct!(sources.iter().copied(), targets.iter().copied())
        {
            all_simple_paths(
                &self.graph,
                source.into(),
                target.into(),
                0,
                Some(ivl.len()),
            )
            .into_iter()
            .map(|path| {
                intersperse(
                    path.into_iter()
                        .filter_map(|x| self.graph.node_weight(x).map(|n| n.hap.to_string())),
                    String::default(),
                )
                .collect::<VString>()
            })
            .for_each(|hap| {
                ret.insert(hap);
            });
        }
        Ok(ret)
    }

    /// Find the highest/lowest positions corresponding to current nodes.
    #[must_use]
    pub fn get_nstart_nend(&self) -> Option<(u32, u32)> {
        let mut ret = (u32::MAX, 0);
        for node in self.node_iter() {
            ret.0 = std::cmp::min(ret.0, node.pos);
            ret.1 = std::cmp::max(ret.1, node.pos);
        }
        if ret.0 == u32::MAX {
            None
        } else {
            Some(ret)
        }
    }

    /// Yield token strings at position.
    /// Could be faster - we just enumerate all nodes and filter for pos currently.
    /// We should be able to use a `BTree` query to do this in linearithmic.
    pub fn get_hap_by_pos(&self, pos: u32) -> impl Iterator<Item = VStr<'_>> {
        self.node_id_map
            .keys()
            .filter(move |x| !x.is_del() && x.pos == pos)
            .map(|x| x.hap.vstr())
    }

    /// Find the next position after `<pos>`
    /// Skip deleted nodes and find the first node that is relevant.
    /// We query the btree for the position in question and walk to the next one to account for arbitary strings in tie-breaking.
    #[must_use]
    pub fn get_next_pos(&self, pos: u32) -> Option<u32> {
        let lower_bound = NodeDatum::new(VString::from("\x7f\x7f\x7f\x7f\x7f\x7f\x7f\x7f"), pos);
        let positions = self.node_id_map.range((
            std::ops::Bound::Excluded(lower_bound),
            std::ops::Bound::Unbounded,
        ));
        positions
            .filter(|x| !x.0.is_del())
            .find(|x| x.0.pos > pos)
            .map(|x| x.0.pos)
    }

    /// Find the first position before `<pos>`
    /// Skip deleted nodes and find the first node that is relevant.
    #[must_use]
    pub fn get_previous_pos(&self, pos: u32) -> Option<u32> {
        let positions = self.node_id_map.range((
            std::ops::Bound::Unbounded,
            std::ops::Bound::Excluded(NodeDatum::from((0u8, pos))),
        ));
        positions.rev().find(|x| !x.0.is_del()).map(|x| x.0.pos)
    }

    /// Retrieve haplotypes at position, converted to owned representations.
    /// For views alone, use `Graph::get_hap_by_pos`.
    #[must_use]
    pub fn get_owned_haps_by_pos(&self, pos: u32) -> Vec<VString> {
        self.get_hap_by_pos(pos)
            .map(VString::from)
            .collect::<Vec<_>>()
    }
}
