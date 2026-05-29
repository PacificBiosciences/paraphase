use itertools::Itertools;
use paraphase::{
    assembly::{
        assembly_result::AssembledPaths,
        variant_graph::{Graph as VariantGraph, SegmentationClass},
    },
    io::json::{ParaphaseOutput, ParsedParaphaseOutputJSON},
    toolkit::{
        range::Range,
        util::{test_file, DError, DResult},
    },
};
use std::collections::*;

#[test]
fn graph_asm_test1() -> Result<(), DError> {
    let file = test_file("jsons/cfc1_graph_test1.json.xz");
    let res = std::process::Command::new("xz")
        .arg("-dc")
        .arg(file)
        .stdout(std::process::Stdio::piped())
        .spawn()?;
    let reader = std::io::BufReader::new(res.stdout.unwrap());
    let json: serde_json::Value = serde_json::from_reader(reader)?;
    let _asm_objects = json
        .as_object()
        .ok_or_else(|| anyhow::anyhow!("Not a dictionary"))?
        .iter()
        .map(|(k, v)| ParaphaseOutput::from_json(v, k.clone(), None))
        .collect::<Vec<_>>();

    let all_inputs = ParsedParaphaseOutputJSON::from_json(json, None)?;
    let gene = "cfc1";
    let gene_input = all_inputs
        .gene_data
        .get(gene)
        .ok_or_else(|| anyhow::anyhow!("cfc1 not found"))?;
    log::debug!("gene input: {gene_input:?}");
    let mut gene_vg = VariantGraph::try_from(gene_input)?;
    let paths = gene_vg.construct()?;
    let paths = paths
        .final_haps
        .into_inner()
        .into_iter()
        .map(|x| x.to_string())
        .sorted()
        .collect::<Vec<_>>();
    assert_eq!(
        paths,
        vec![
            String::from("111211"),
            String::from("212211"),
            String::from("xxx122")
        ]
    );
    Ok(())
}

#[test]
fn graph_asm_test2() -> DResult {
    use crate::SegmentationClass::*;
    let file = test_file("jsons/cfc1_graph_test2.json.xz");
    let all_inputs = ParsedParaphaseOutputJSON::from_path(file, None)?;
    let gene = "cfc1";
    let gene_input = all_inputs
        .gene_data
        .get(gene)
        .ok_or_else(|| anyhow::anyhow!("cfc1 not found"))?;
    // let pivot_pos = gene_input.pivot_index(None);
    // log::info!("Pivot pos: {pivot_pos:?}");
    let mut graph = VariantGraph::try_from(gene_input)?;
    graph.init()?;
    let found_edges = graph.formatted_edges();
    let expected: Vec<Vec<(String, String)>> = [
        vec![("1-0", "1-1"), ("2-0", "1-1")],
        vec![("1-1", "2-2"), ("1-1", "1-2")],
        vec![("2-2", "1-3"), ("1-2", "1-3"), ("1-2", "2-3")],
        vec![("2-3", "1-4"), ("1-3", "2-4"), ("1-3", "1-4")],
    ]
    .into_iter()
    .map(|subvec| {
        subvec
            .into_iter()
            .map(|(from, to)| (String::from(from), String::from(to)))
            .sorted()
            .collect::<Vec<_>>()
    })
    .collect::<Vec<_>>();
    assert_eq!(found_edges, expected);
    assert_eq!(
        graph.formatted_nodes(),
        vec!["1-0", "1-1", "2-0", "2-2", "1-2", "1-3", "2-3", "1-4", "2-4"]
            .into_iter()
            .sorted()
            .collect::<Vec<_>>()
    );
    // init graph worked.
    // Now summarize edges
    graph.summarize_edges()?;
    let expected_dnhap = [(0, Ten), (1, Ten), (2, Ten), (3, Ten)]
        .into_iter()
        .collect::<BTreeMap<u32, SegmentationClass>>();
    // Check nodes/hap after summarizing edges
    assert_eq!(expected_dnhap, graph.dnhap);
    assert_eq!(
        graph.formatted_nodes(),
        vec!["1-0", "1-1", "2-0", "2-2", "1-2", "1-3", "2-3", "1-4", "2-4"]
            .into_iter()
            .sorted()
            .collect::<Vec<_>>()
    );
    let segments = graph.segment_graph()?;
    assert_eq!(
        segments,
        [(Range::<u32>::new(0, 3), Ten)]
            .into_iter()
            .collect::<BTreeMap<_, _>>()
    );
    assert_eq!(
        graph.formatted_nodes(),
        vec!["1-0", "1-1", "2-0", "2-2", "1-2", "1-3", "2-3", "1-4", "2-4"]
            .into_iter()
            .sorted()
            .collect::<Vec<_>>()
    );

    let found_edges = graph.formatted_edges();
    let format_expected_edges = |x: &[Vec<(&str, &str, i32)>]| -> Vec<Vec<(String, String)>> {
        x.iter()
            .map(|subvec| {
                subvec
                    .iter()
                    .copied()
                    .sorted()
                    .unique()
                    .map(|(from, to, _count)| (String::from(from), String::from(to)))
                    .sorted()
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    };
    let expected_edges = format_expected_edges(&[
        vec![
            ("1-0", "1-1", 37),
            ("2-0", "1-1", 13),
            ("1-0", "1-1", 37),
            ("2-0", "1-1", 13),
        ],
        vec![
            ("1-1", "2-2", 2),
            ("1-1", "1-2", 10),
            ("1-1", "2-2", 2),
            ("1-1", "1-2", 10),
        ],
        vec![
            ("2-2", "1-3", 10),
            ("1-2", "1-3", 23),
            ("1-2", "2-3", 11),
            ("2-2", "1-3", 10),
            ("1-2", "1-3", 23),
            ("1-2", "2-3", 11),
        ],
        vec![
            ("2-3", "1-4", 8),
            ("1-3", "2-4", 16),
            ("1-3", "1-4", 7),
            ("2-3", "1-4", 8),
            ("1-3", "2-4", 16),
            ("1-3", "1-4", 7),
        ],
    ]);
    assert_eq!(found_edges, expected_edges);
    graph.merge_simple_edges_from_segments(&segments)?;
    graph.merge_complex_regions(segments)?;

    let expected_pos_edge_map = [
        vec![
            ("1-0", "1-1", 37),
            ("2-0", "1-1", 13),
            ("1-0", "1-1", 37),
            ("2-0", "1-1", 13),
        ],
        vec![
            ("1-1", "2-2", 2),
            ("1-1", "1-2", 10),
            ("1-1", "2-2", 2),
            ("1-1", "1-2", 10),
        ],
        vec![
            ("2-2", "1-3", 10),
            ("1-2", "1-3", 23),
            ("1-2", "2-3", 11),
            ("2-2", "1-3", 10),
            ("1-2", "1-3", 23),
            ("1-2", "2-3", 11),
        ],
        vec![
            ("2-3", "1-4", 8),
            ("1-3", "2-4", 16),
            ("1-3", "1-4", 7),
            ("2-3", "1-4", 8),
            ("1-3", "2-4", 16),
            ("1-3", "1-4", 7),
        ],
    ];
    log::debug!(
        "pos_edge_map after merging simple + complex edges: {:?}",
        graph.pos_edge_map
    );
    let expected_pos_edge_map = format_expected_edges(&expected_pos_edge_map);
    let found_edges = graph.formatted_edges();
    assert_eq!(found_edges, expected_pos_edge_map);
    // Now we have done init, summarize edges, segment_graph, and merge_simple_edges_from_segments and merge_complex_regiong
    // Next, we merge where possible.
    let (nstart, nend) = graph.get_nstart_nend().expect("No nodes in graph?");
    assert!(nstart != nend);
    graph.merge_where_possible(nstart, nend)?;

    let expected_pos_edge_map = format_expected_edges(&[vec![
        ("1121-0", "2-4", 0),
        ("1121-0", "1-4", 0),
        ("2111-0", "2-4", 0),
        ("2111-0", "1-4", 0),
        ("1111-0", "2-4", 0),
        ("1111-0", "1-4", 0),
        ("1112-0", "1-4", 0),
    ]]);

    let found_edges = graph.formatted_edges();
    assert_eq!(found_edges, expected_pos_edge_map);

    // And now let's finish assemble_haps
    graph.merge_blocks_by_size()?;
    let (nstart, nend) = graph.get_nstart_nend().expect("Nodes expected in graph");
    assert_eq!((nstart, nend), (0, 4));
    let haps = {
        let (_success, haps, _x, _y) = graph.subregion_assembly(nstart, nend, true)?;
        AssembledPaths::from_seqs(haps)
    };
    assert_eq!(
        haps,
        AssembledPaths::from_seqs(["1111x", "21112", "11121", "11211"].into_iter())
    );

    Ok(())
}
