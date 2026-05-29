use crate::assembly::node_datum::NodeDatum;
use crate::assembly::variant_graph::{Graph, IxType};
use crate::toolkit::util::HashMap;

use itertools::Itertools;
use petgraph::graph::EdgeIndex;

impl Graph {
    /// `paraphase_format_node`: format the paraphase format (stored in graph) to the string representation used in paraphase.
    /// This is mostly for readability and testing.
    #[must_use]
    pub fn paraphase_format_node(
        &self,
        node_index: impl TryInto<IxType>,
        node_id_lookup: Option<&HashMap<u32, &NodeDatum>>,
    ) -> Option<String> {
        let node_index = node_index.try_into().ok()?;
        let node_id_local = if node_id_lookup.is_none() {
            Some(self.get_id_lookup())
        } else {
            None
        };
        let node_id_lookup = if let Some(node_id_lookup) = node_id_lookup {
            node_id_lookup
        } else {
            node_id_local.as_ref()?
        };
        node_id_lookup
            .get(&node_index)
            .map(|node| format!("{}-{}", node.hap, node.pos))
    }

    /// `paraphase_format_edge`: format the paraphase format (stored in graph) to the string representation used in paraphase.
    /// This is mostly for readability and testing.
    #[must_use]
    pub fn paraphase_format_edge(
        &self,
        edge: EdgeIndex<IxType>,
        node_id_lookup: Option<&HashMap<u32, &NodeDatum>>,
    ) -> Option<(String, String)> {
        self.graph
            .edge_endpoints(edge)
            .map(|(source, target)| {
                (
                    self.paraphase_format_node(source.index(), node_id_lookup),
                    self.paraphase_format_node(target.index(), node_id_lookup),
                )
            })
            .and_then(|x| match x {
                (Some(y), Some(z)) => Some((y, z)),
                _ => None,
            })
    }

    /// Display all nodes from the graph as a vector of strings.
    #[must_use]
    pub fn formatted_nodes(&self) -> Vec<String> {
        self.node_iter()
            .filter_map(|x| self.node_id_map.get(x).copied())
            .filter_map(|id| self.paraphase_format_node(id, None))
            .sorted()
            .collect::<Vec<_>>()
    }

    /// Display an edge as a string.
    #[must_use]
    pub fn format_edge(&self, x: &(u32, u32, u32)) -> Option<(String, String)> {
        let edge_index = self.find_edge(x.0.into(), x.1.into());
        edge_index.and_then(|edge_index| self.paraphase_format_edge(edge_index, None))
    }

    /// Display all nodes from the graph as a vector of string pairs in (from, to) format.
    #[must_use]
    pub fn formatted_edges(&self) -> Vec<Vec<(String, String)>> {
        self.pos_edge_map
            .values()
            .map(|edges| {
                edges
                    .iter()
                    .filter_map(|edge| self.format_edge(edge))
                    .sorted()
                    .collect::<Vec<(String, String)>>()
            })
            .collect::<Vec<Vec<(String, String)>>>()
    }
}
