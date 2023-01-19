import pytest
from paraphase.haplotype_assembler import VariantGraph


class TestVariantGraph(object):
    def test_get_thres(self):
        reads = {"r1": "11", "r2": "22"}
        graph = VariantGraph(reads, 0)
        num_reads = [1, 10, 10]
        assert graph.get_thres(num_reads) == 1

    def test_filter_low_support_haps(self):
        reads = {
            "r1": "11x",
            "r2": "22x",
            "r3": "x12",
            "r4": "x22",
            "r5": "x2x",
            "r6": "x22",
            "r7": "x11",
            "r8": "xx1",
            "r9": "xx1",
        }
        graph = VariantGraph(reads, 0)
        haps = ["112", "111", "222"]
        filtered_haps = graph.filter_low_support_haps(haps, min_count=3)
        assert sorted(filtered_haps) == ["111", "222"]
        filtered_haps = graph.filter_low_support_haps(haps, min_count=1)
        assert sorted(filtered_haps) == ["111", "112", "222"]

    def test_initialize_graph(self):
        reads = {
            "r1": "11x",
            "r2": "11x",
            "r3": "22x",
            "r4": "22x",
            "r5": "x12",
            "r6": "x22",
            "r7": "x2x",
            "r8": "x22",
            "r9": "x11",
            "r10": "x11",
            "r11": "x22",
            "r12": "x22",
            "r13": "x22",
            "r14": "x22",
            "r15": "x11",
            "r16": "x11",
            "r17": "x11",
            "r18": "x11",
        }
        graph = VariantGraph(reads, 0)
        graph.initialize_graph()
        assert sorted(graph.nodes) == ["1-0", "1-1", "1-2", "2-0", "2-1", "2-2"]
        assert graph.edges == [
            ["1-0", "1-1"],
            ["2-0", "2-1"],
            ["2-1", "2-2"],
            ["1-1", "1-2"],
        ]
        assert graph.edge_per_position == {
            0: [("1-0", "1-1", 2), ("2-0", "2-1", 2)],
            1: [("2-1", "2-2", 6), ("1-1", "1-2", 6)],
        }
        assert "r5" not in graph.reads

        reads = {
            "r1": "11x",
            "r2": "11x",
            "r3": "22x",
            "r4": "22x",
            "r5": "x12",
            "r6": "x22",
            "r7": "x2x",
            "r8": "x22",
            "r9": "x11",
            "r10": "x11",
            "r11": "x22",
            "r12": "x22",
            "r13": "x22",
            "r14": "x22",
            "r15": "x11",
            "r16": "x11",
            "r17": "x11",
        }
        graph = VariantGraph(reads, 0)
        graph.initialize_graph()
        assert sorted(graph.nodes) == ["1-0", "1-1", "1-2", "2-0", "2-1", "2-2"]
        assert graph.edges == [
            ["1-0", "1-1"],
            ["2-0", "2-1"],
            ["1-1", "2-2"],
            ["2-1", "2-2"],
            ["1-1", "1-2"],
        ]
        assert graph.edge_per_position == {
            0: [("1-0", "1-1", 2), ("2-0", "2-1", 2)],
            1: [("1-1", "2-2", 1), ("2-1", "2-2", 6), ("1-1", "1-2", 5)],
        }
        assert "r5" in graph.reads

    # test one variant, two variants
    def test_run(self):
        reads = {"r1": "11", "r2": "22", "r3": "11", "r4": "22", "r5": "11", "r6": "22"}
        graph = VariantGraph(reads, pivot_pos=1)
        graph.initialize_graph()
        haps = graph.assemble_haplotypes(debug=False, make_plot=False)
        assert sorted(haps) == ["11", "22"]

    def test_get_highest_cn(self):
        reads = {"r1": "11111", "r2": "22222"}
        graph = VariantGraph(reads, 0)
        highest_cn, _ = graph.get_highest_cn(["11111", "22222"])
        assert highest_cn == 2
        highest_cn, _ = graph.get_highest_cn(["121xxx", "111111", "xxx111"])
        assert highest_cn == 2
        highest_cn, _ = graph.get_highest_cn(["1211xx", "111111", "xxx111"])
        assert highest_cn == 3

    def test_get_next_pos(self):
        reads = {"r1": "11", "r2": "22"}
        graph = VariantGraph(reads, 0)
        graph.nodes = [
            "121-0",
            "111-0",
            "222-3",
            "111-3",
        ]
        next_pos, next_haps = graph.get_next_pos(1)
        assert next_pos == 3
        assert next_haps == ["222", "111"]

        next_pos, next_haps = graph.get_next_pos(0)
        assert next_pos == 3
        assert next_haps == ["222", "111"]

        next_pos, next_haps = graph.get_next_pos(3)
        assert next_pos is None
        assert next_haps is None

    def test_get_previous_pos(self):
        reads = {"r1": "11", "r2": "22"}
        graph = VariantGraph(reads, 0)
        graph.nodes = [
            "121-0",
            "111-0",
            "222-3",
            "111-3",
        ]
        next_pos, next_haps = graph.get_previous_pos(1)
        assert next_pos == 0
        assert next_haps == ["121", "111"]

        next_pos, next_haps = graph.get_previous_pos(3)
        assert next_pos == 0
        assert next_haps == ["121", "111"]

        next_pos, next_haps = graph.get_previous_pos(0)
        assert next_pos is None
        assert next_haps is None

    def test_path(self):
        reads = {"r1": "11", "r2": "22"}
        graph = VariantGraph(reads, 0)
        graph.initialize_graph()
        paths = graph.path(0, 1)
        assert sorted(paths) == ["11", "22"]

    def test_summerize_edges(self):
        reads = {"r1": "11", "r2": "22"}
        graph = VariantGraph(reads, 0)
        graph.nodes = [
            "1-0",
            "2-0",
            "1-1",
            "2-1",
            "1-2",
            "2-2",
            "1-3",
            "2-3",
            "1-4",
            "2-4",
        ]
        graph.edge_per_position = {
            0: [("1-0", "1-1", 2), ("2-0", "2-1", 3)],
            1: [("1-1", "2-2", 2), ("2-1", "2-2", 6), ("1-1", "1-2", 5)],
            2: [("1-2", "2-3", 2), ("2-2", "2-3", 6), ("1-2", "1-3", 5)],
            3: [("1-3", "2-4", 2), ("2-3", "2-4", 6), ("1-3", "1-4", 5)],
        }
        graph.edge_info = {
            0: [2, 3],
            1: [2, 3],
            2: [2, 3],
            3: [2, 3],
        }
        graph.summerize_edges()
        assert graph.dnhap == {
            0: 2,
            1: 10,
            2: 10,
            3: 10,
        }

    def test_get_segments(self):
        reads = {"r1": "11", "r2": "22"}
        graph = VariantGraph(reads, 0)
        graph.dnhap = {
            0: 2,
            1: 10,
            2: 10,
            3: 10,
        }
        assert graph.get_segments() == {
            (0, 0): 2,
            (1, 3): 10,
        }

        graph.dnhap = {
            0: 2,
            1: 2,
            2: 10,
            3: 10,
        }
        assert graph.get_segments() == {
            (0, 1): 2,
            (2, 3): 10,
        }

        graph.dnhap = {
            0: 10,
            1: 2,
            2: 2,
            3: 2,
        }
        assert graph.get_segments() == {
            (0, 3): 2,
        }

        graph.dnhap = {
            0: 2,
            1: 2,
            2: 2,
            3: 10,
        }
        assert graph.get_segments() == {
            (0, 3): 2,
        }

        graph.dnhap = {
            0: 10,
            1: 2,
            2: 2,
            3: 10,
        }
        assert graph.get_segments() == {
            (0, 0): 10,
            (1, 2): 2,
            (3, 3): 10,
        }

    def test_rescue_missing(self):
        left_missing = ["11"]
        right_missing = ["22"]
        dnext = {"11": ["22"]}
        dbefore = {"22": ["11"]}
        rescued_haps = VariantGraph.rescue_missing(
            left_missing, right_missing, dnext, dbefore
        )
        assert rescued_haps == {"1122"}

        left_missing = ["111"]
        right_missing = ["122"]
        dnext = {"121": {"111", "122"}, "111": {"111", "122"}}
        dbefore = {"111": {"111", "121"}, "122": {"111", "121"}}
        rescued_haps = VariantGraph.rescue_missing(
            left_missing, right_missing, dnext, dbefore
        )
        assert rescued_haps == {"111122"}

    def test_get_missing(self):
        reads = {"r1": "11", "r2": "22"}
        graph = VariantGraph(reads, 0)
        haps_left = ["11", "22"]
        haps_right = ["12", "21"]
        test_haps = ["1112"]
        left_missing, right_missing = graph.get_missing(
            test_haps, haps_left, haps_right
        )
        assert left_missing == ["22"]
        assert right_missing == ["21"]

    def test_compare_two_haps(self):
        match, mismatch, extend = VariantGraph.compare_two_haps("12x", "x22")
        assert match == 1
        assert mismatch == 0
        assert extend == 1

        match, mismatch, extend = VariantGraph.compare_two_haps("12xx", "xx22")
        assert match == 0

    def test_get_matching_hap(self):
        nodes = {0: ["11", "22"], 2: ["11", "22"]}
        hap1, best_match, max_match = VariantGraph.get_matching_hap(nodes, 0, "x11x")
        assert hap1 == "x1"
        assert best_match == ["11"]
        assert max_match == 1
        hap1, best_match, max_match = VariantGraph.get_matching_hap(nodes, 2, "x11x")
        assert hap1 == "1x"
        assert best_match == ["11"]
        assert max_match == 1

    def test_merge_two_pos(self):
        reads = {
            "r1": "x11x",
            "r2": "11xx",
            "r3": "x11x",
            "r4": "x22x",
            "r5": "xx22",
        }
        graph = VariantGraph(reads, 0)
        nodes = {0: ["11", "22"], 2: ["11", "22"]}
        sub_hap_support_reads, dnext, dbefore = graph.merge_two_pos(nodes)
        assert sub_hap_support_reads == {"1111": ["r1", "r3"], "2222": ["r4"]}
        assert dnext == {"11": {"11"}, "22": {"22"}}
        assert dbefore == {"11": {"11"}, "22": {"22"}}

        reads = {
            "r1": "111x",
            "r2": "211x",
            "r3": "x111",
            "r4": "x112",
            "r5": "xx12",
            "r6": "1111",
        }
        graph = VariantGraph(reads, 0)
        nodes = {0: ["11", "21"], 2: ["11", "12"]}
        sub_hap_support_reads, dnext, dbefore = graph.merge_two_pos(nodes)
        assert sub_hap_support_reads == {"1111": ["r6"]}
        assert dnext == {"11": {"11", "12"}, "21": {"11", "12"}}
        assert dbefore == {"11": {"11", "21"}, "12": {"11", "21"}}

    def test_subregion_assembly(self):
        reads = {
            "r1": "x211xx",
            "r2": "x112xx",
            "r3": "xxx22x",
            "r4": "121xxx",
            "r5": "xx1222",
            "r6": "1211xx",
            "r7": "xxx22x",
        }
        graph = VariantGraph(reads, 0)
        graph.nodes = [
            "121-0",
            "111-0",
            "222-3",
            "111-3",
        ]
        success, ass_haps, dnext, dbefore = graph.subregion_assembly(0, 3)
        assert success is True
        assert sorted(ass_haps) == ["111222", "121111"]

        # more ambiguity
        reads = {
            "r1": "x2111x",
            "r2": "x11xxx",
            "r3": "xxx22x",
            "r4": "x2111x",
            "r5": "xxx122",
            "r6": "1211xx",
            "r7": "xx1122",
            "r8": "xx11xx",
        }
        graph = VariantGraph(reads, 0)
        graph.nodes = [
            "121-0",
            "111-0",
            "122-3",
            "111-3",
        ]
        success, ass_haps, dnext, dbefore = graph.subregion_assembly(0, 3)
        assert success is False
        assert dnext == {"121": {"111", "122"}, "111": {"111", "122"}}
        assert dbefore == {"111": {"111", "121"}, "122": {"111", "121"}}
        success, ass_haps, dnext, dbefore = graph.subregion_assembly(0, 3, allow_x=True)
        assert sorted(ass_haps) == ["111xxx", "121111", "xxx122"]
        # enable rescue missing haps
        graph.edge_info.setdefault(2, [3, 4, 5])
        graph.dnhap.setdefault(2, 10)
        success, ass_haps, dnext, dbefore = graph.subregion_assembly(0, 3, allow_x=True)
        assert success is True
        assert sorted(ass_haps) == ["111122", "121111"]

        # only one possible node at one position
        reads = {
            "r1": "x21x",
            "r2": "x11x",
            "r3": "xxx1",
            "r4": "121x",
            "r5": "xx1x",
            "r6": "121x",
        }
        graph = VariantGraph(reads, 0)
        graph.nodes = [
            "121-0",
            "111-0",
            "1-3",
        ]
        success, ass_haps, dnext, dbefore = graph.subregion_assembly(0, 3)
        assert success is True
        assert sorted(ass_haps) == ["1111", "1211"]

    def test_rm_add_edges(self):
        reads = {
            "r1": "x11x",
            "r2": "11xx",
            "r3": "x11x",
            "r4": "x22x",
            "r5": "xx22",
        }
        graph = VariantGraph(reads, 0)
        graph.nodes = ["1-0", "1-1", "1-2", "1-3", "2-0", "2-1", "2-2", "2-3"]
        graph.edges = [
            ["1-0", "1-1"],
            ["2-0", "2-1"],
            ["1-1", "1-2"],
            ["2-1", "2-2"],
            ["1-2", "1-3"],
            ["2-2", "2-3"],
        ]
        graph.edge_per_position = {
            0: [("1-0", "1-1", 2), ("2-0", "2-1", 2)],
            1: [("1-1", "1-2", 2), ("2-1", "2-2", 2)],
            2: [("1-2", "1-3", 2), ("2-2", "2-3", 2)],
        }
        graph.rm_add_edges(1, 2, ["11", "22"])
        assert graph.nodes == ["1-0", "1-3", "2-0", "2-3", "11-1", "22-1"]
        assert sorted(graph.edges) == [
            ["1-0", "11-1"],
            ["11-1", "1-3"],
            ["2-0", "22-1"],
            ["22-1", "2-3"],
        ]
        assert graph.edge_per_position == {
            0: [("1-0", "11-1", None), ("2-0", "22-1", None)],
            1: [("11-1", "1-3", None), ("22-1", "2-3", None)],
        }

        # at the edge
        graph.nodes = ["1-0", "1-1", "1-2", "1-3", "2-0", "2-1", "2-2", "2-3"]
        graph.edges = [
            ["1-0", "1-1"],
            ["2-0", "2-1"],
            ["1-1", "1-2"],
            ["2-1", "2-2"],
            ["1-2", "1-3"],
            ["2-2", "2-3"],
        ]
        graph.edge_per_position = {
            0: [("1-0", "1-1", 2), ("2-0", "2-1", 2)],
            1: [("1-1", "1-2", 2), ("2-1", "2-2", 2)],
            2: [("1-2", "1-3", 2), ("2-2", "2-3", 2)],
        }
        graph.rm_add_edges(0, 1, ["11", "22"])
        assert sorted(graph.nodes) == ["1-2", "1-3", "11-0", "2-2", "2-3", "22-0"]
        assert sorted(graph.edges) == [
            ["1-2", "1-3"],
            ["11-0", "1-2"],
            ["2-2", "2-3"],
            ["22-0", "2-2"],
        ]
        assert graph.edge_per_position == {
            0: [("11-0", "1-2", None), ("22-0", "2-2", None)],
            2: [("1-2", "1-3", 2), ("2-2", "2-3", 2)],
        }

    ## need some changes here
    def test_extend_pivot_blocks(self):
        reads = {
            "r1": "x2111x",
            "r2": "x11xxx",
            "r3": "xxx22x",
            "r4": "x2111x",
            "r5": "xxx122",
            "r6": "1211xx",
            "r7": "xx1122",
            "r8": "xx11xx",
            "r9": "x111x1",
            "r0": "x1x12x",
        }
        graph = VariantGraph(reads, 0)
        graph.nodes = [
            "121-0",
            "111-0",
            "122-3",
            "111-3",
        ]
        graph.dnhap = {2: 10}
        graph.edge_info = {2: [2, 3]}
        haps = ["111xxx"]
        extended_haps = graph.extend_pivot_blocks(haps)
        assert extended_haps == ["1111xx"]

    def test_order_blocks_by_size(self):
        reads = {
            "r1": "111x",
            "r2": "211x",
            "r3": "x111",
            "r4": "x112",
            "r5": "xx12",
            "r6": "1111",
        }
        graph = VariantGraph(reads, 0)
        graph.nodes = ["11-0", "22-0", "11-2", "22-2", "111-4", "222-4"]
        block_order = graph.order_blocks_by_size()
        assert block_order == [(2, 4, 2, 3), (0, 2, 2, 2)]

    def test_match_reads_and_haplotypes(self):
        hap_list = ["21xx", "xx12"]
        haplotype_per_read = {"r1": "x11x", "r2": "x11x"}
        read_support = VariantGraph.match_reads_and_haplotypes(
            haplotype_per_read, hap_list
        )
        assert read_support.unique == {
            "21xx": ["x11x", "x11x"],
            "xx12": ["x11x", "x11x"],
        }
