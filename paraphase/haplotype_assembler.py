# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


import networkx as nx
import matplotlib.pyplot as plt
import copy
from pprint import pprint
import math
import numpy as np
from collections import namedtuple


class VariantGraph:
    HaplotypesReads = namedtuple("HaplotypesReads", "unique nonunique by_read")

    def __init__(self, reads, pivot_pos, figure_id="fig"):
        self.nvar = len(list(reads.values())[0])
        if self.nvar <= 1:
            raise Exception("Too few variants")
        self.reads = reads
        self.reads_original = copy.deepcopy(reads)
        self.pivot_pos = pivot_pos
        self.figure_id = figure_id
        # these three get updated throughout
        self.edge_per_position = {}
        self.nodes = []
        self.edges = []
        # only used when first constructing edges
        self.next_per_node = {}
        self.previous_per_node = {}
        # positions to plot nodes
        self.node_positions = {}
        self.edge_info = {}
        self.dnhap = {}

    def add_edge(self, a, b):
        """Add edges between two nodes"""
        self.edges.append([a, b])

    def rm_edge(self, a, b):
        """Remove edges between two nodes"""
        if [a, b] in self.edges:
            self.edges.remove([a, b])

    def get_positions(self):
        """Get a dictionary of the positions of each node in the graph plot"""
        for i in [1, 2]:
            for j in range(self.nvar):
                self.node_positions.setdefault(f"{i}-{j}", (j * 1000, i * 10))
        for n in self.nodes:
            i, j = n.split("-")
            if len(i) == 1 and f"{i}-{j}" not in self.node_positions:
                i = int(i)
                j = int(j)
                self.node_positions.setdefault(f"{i}-{j}", (j * 1000, i * 10))

    def label_pos_offset(self):
        """Add offset to labels in the graph plot"""
        pos_with_offset = {}
        y_off = 0.5  # offset on the y axis
        for k in self.node_positions.keys():
            v = self.node_positions[k]
            if int(k.split("-")[0]) == 1:
                pos_with_offset[k] = (v[0], v[1] - y_off)
            else:
                pos_with_offset[k] = (v[0], v[1] + y_off)
        return pos_with_offset

    def visualize(self, ax):
        """Make a graph plot"""
        G = nx.DiGraph()
        for n in self.nodes:
            G.add_node(n)
        G.add_edges_from(self.edges)
        nx.draw_networkx(
            G,
            pos=self.node_positions,
            with_labels=False,
            node_size=40,
            node_color="red",
            ax=ax,
            arrows=False,
        )
        label_pos = self.label_pos_offset()
        text = nx.draw_networkx_labels(
            G, pos=label_pos, font_size=25, ax=ax, clip_on=False
        )
        for _, t in text.items():
            t.set_rotation(25)
        plt.plot()

    def path(self, pos1, pos2):
        """Get simple path between two nodes"""
        G = nx.DiGraph()
        G.add_edges_from(self.edges)
        sources = set([a[0] for a in self.edge_per_position[pos1]])
        targets = set([a[1] for a in self.edge_per_position[pos2 - 1]])
        paths = []
        for source in sources:
            for target in targets:
                kk = list(
                    nx.all_simple_edge_paths(G, source, target, cutoff=pos2 - pos1)
                )
                for path in kk:
                    hap = path[0][0].split("-")[0]
                    for site in path:
                        hap += site[1].split("-")[0]
                    paths.append("".join(hap))
        return list(set(paths))

    def plotfig(self):
        """Create figure"""
        figure = plt.gcf()
        figure.set_size_inches(60, 30)
        plt.savefig(f"{self.figure_id}.png")
        plt.close(figure)

    def get_hap_by_pos(self, pos):
        """Given a position, return all the haplotypes starting there"""
        haps = []
        for node in self.nodes:
            if int(node.split("-")[1]) == pos:
                haps.append(node.split("-")[0])
        return haps

    def get_nstart_nend(self):
        """Get the smallest and biggest position numbers in the graph"""
        pos_left = sorted(list(set([int(a.split("-")[1]) for a in self.nodes])))
        if pos_left == []:
            return 0, 0
        nstart = min(pos_left)
        nend = max(pos_left)
        return nstart, nend

    def get_next_pos(self, pos):
        """Get next position and its haplotype blocks in the graph"""
        positions = sorted(list(set([int(a.split("-")[1]) for a in self.nodes])))
        if pos >= positions[-1]:
            return None, None
        for position in positions:
            if position > pos:
                next_pos = position
                break
        next_haps = self.get_hap_by_pos(next_pos)
        return next_pos, next_haps

    def get_previous_pos(self, pos):
        """Get previous position and its haplotype blocks in the graph"""
        positions = sorted(list(set([int(a.split("-")[1]) for a in self.nodes])))
        if pos <= positions[0]:
            return None, None
        for position in positions[::-1]:
            if position < pos:
                prev_pos = position
                break
        prev_haps = self.get_hap_by_pos(prev_pos)
        return prev_pos, prev_haps

    def get_highest_cn(self, final_haps):
        """Return the highest possible copy number given the assembled haplotypes"""
        main_haps_candidates = []
        for i in reversed(range(self.nvar)):
            lhaps = [a for a in final_haps if a[i] != "x"]
            if lhaps not in main_haps_candidates:
                main_haps_candidates.append(lhaps)
        main_haps_candidates = sorted(
            main_haps_candidates, key=lambda x: len(x), reverse=True
        )
        main_haps = main_haps_candidates[0]
        highest_cn = len(main_haps)
        return (
            highest_cn,
            main_haps_candidates,
        )

    def filter_low_support_haps(self, haps_to_assess, min_count=3, debug=False):
        """Remove haplotypes that have a low number of supporting reads"""
        reads = self.reads_original
        read_support = self.match_reads_and_haplotypes(reads, haps_to_assess)
        good_reads = read_support.unique
        filtered_ass_haps = []
        for test_hap in good_reads:
            nread = len(good_reads[test_hap])
            if nread > 0:
                filtered_ass_haps.append(test_hap)
            elif debug:
                print(f"removing low support hap {test_hap} {nread}")
        haps_to_assess = filtered_ass_haps
        while True:
            read_support = self.match_reads_and_haplotypes(reads, haps_to_assess)
            good_reads = read_support.unique
            filtered_ass_haps = []
            for test_hap in good_reads:
                nread = len(good_reads[test_hap])
                if nread >= min_count:
                    filtered_ass_haps.append(test_hap)
                elif debug:
                    print(f"removing low support hap {test_hap} {nread}")
            if len(filtered_ass_haps) == len(haps_to_assess):
                break
            else:
                haps_to_assess = filtered_ass_haps
        if debug:
            print(
                "All removed haplotypes are ",
                set(haps_to_assess) - set(filtered_ass_haps),
            )
        return filtered_ass_haps

    @staticmethod
    def get_thres(num_reads):
        """Get thresholds for pruning some noisy links"""
        min_num = min(num_reads) if num_reads != [] else 0
        remain = [a for a in num_reads if a != min_num]
        if remain != []:
            if min(remain) > 10:
                return min(5, max(2, math.floor(min(remain) / 25)))
            if min(remain) > 5:
                return 1
        return 0

    # update the edges
    def get_edges(self, debug=False):
        """Create and filter edges"""
        dedge = {}
        reads = self.reads
        for read in reads:
            hap = reads[read]
            nvar = self.nvar
            if hap.count("x") < nvar - 1:
                for i in range(nvar - 1):
                    base1 = hap[i]
                    base2 = hap[i + 1]
                    if "x" not in [base1, base2]:
                        node1 = base1 + "-" + str(i)
                        node2 = base2 + "-" + str(i + 1)
                        dedge.setdefault((node1, node2), []).append(read)
        for node1, node2 in dedge:
            pos = int(node1[2:])
            self.edge_per_position.setdefault(pos, []).append(
                (node1, node2, len(dedge[(node1, node2)]))
            )
            self.next_per_node.setdefault(node1, []).append(node2)
            self.previous_per_node.setdefault(node2, []).append(node1)
        new_edge_per_position = {}
        for i in range(self.nvar - 1):
            new_edge_per_position.setdefault(i, [])

        # only remove those edges that have low support and won't cause an
        # orphan node if removed
        edges_removed = []
        for pos in self.edge_per_position:
            num_reads = [at[2] for at in self.edge_per_position[pos]]
            # when variants are inside a deletion region
            has_sv = True in [
                at[0].split("-")[0] not in ["1", "2"]
                or at[1].split("-")[0] not in ["1", "2"]
                for at in self.edge_per_position[pos]
            ]
            thres = self.get_thres(num_reads)
            for node1, node2, n in self.edge_per_position[pos]:
                if (
                    len(self.next_per_node[node1]) > 1
                    and len(self.previous_per_node[node2]) > 1
                    and n <= thres
                ) or (
                    has_sv
                    and n == 1
                    and len(self.next_per_node[node1]) > 1
                    and len(self.previous_per_node[node2]) > 1
                ):
                    edges_removed.append((node1, node2, n))
                else:
                    new_edge_per_position.setdefault(pos, []).append((node1, node2, n))
        self.edge_per_position = new_edge_per_position
        for node1, node2, n in edges_removed:
            if debug:
                print("removing", node1, node2, n)
            self.next_per_node[node1].remove(node2)
            self.previous_per_node[node2].remove(node1)
        for pos in self.edge_per_position:
            num_reads = [at[2] for at in self.edge_per_position[pos]]
            self.edge_info.setdefault(pos, num_reads)

    def summerize_edges(self):
        """Check the number of different edges starting at each position"""
        edge_per_position = self.edge_per_position
        for pos in edge_per_position:
            nhap = -1
            if (
                pos in self.edge_info
                and self.edge_info[pos] != []
                and min(self.edge_info[pos]) >= 2
                and max(self.edge_info[pos]) > 2
            ):
                dlinks = {}
                site1 = set()
                site2 = set()
                expected_site1 = set()
                expected_site2 = set()
                for node in self.nodes:
                    base, base_pos = node.split("-")
                    base_pos = int(base_pos)
                    if base_pos == pos:
                        expected_site1.add(base)
                    elif base_pos == pos + 1:
                        expected_site2.add(base)
                for node1, node2, _ in edge_per_position[pos]:
                    site1.add(node1.split("-")[0])
                    site2.add(node2.split("-")[0])
                    dlinks.setdefault(node1, set()).add(node2)
                # add some test here
                if sorted(site1) == sorted(expected_site1) and sorted(site2) == sorted(
                    expected_site2
                ):
                    if [a for a in dlinks if len(dlinks[a]) > 1] == []:
                        nhap = 2
                    else:
                        nhap = 10
            self.dnhap.setdefault(pos, nhap)
        # break the graph at the pivot pos, can make it optional
        if self.pivot_pos - 1 in self.dnhap:
            if self.dnhap[self.pivot_pos - 1] == 2:
                self.dnhap[self.pivot_pos - 1] = 10

    def get_segments(self):
        """Break graph into segments"""
        segments = {}
        dnhap = dict(sorted(self.dnhap.items()))
        dnhap_bk = copy.deepcopy(self.dnhap)
        list_keys = list(dnhap_bk.keys())
        list_values = list(dnhap_bk.values())
        if list_values[:3] == [10, 2, 2] and list_values != [10, 2, 2, 10]:
            dnhap[list_keys[0]] = 2
        if list_values[-3:] == [2, 2, 10] and list_values != [10, 2, 2, 10]:
            dnhap[list_keys[-1]] = 2
        for i, val in enumerate(list_values):
            if i > 0 and i < len(list_values) - 1:
                if list_values[i - 1] == 10 and val == 2 and list_values[i + 1] == 10:
                    dnhap[list_keys[i]] = 10

        lpos = [0]
        last_nhap = dnhap[0]
        for pos in dnhap:
            if pos > 0:
                nhap = dnhap[pos]
                if nhap != last_nhap:
                    segments.setdefault((min(lpos), max(lpos)), last_nhap)
                    lpos = [pos]
                    last_nhap = nhap
                else:
                    lpos.append(pos)
        segments.setdefault((min(lpos), max(lpos)), last_nhap)
        return segments

    def update_reads(self):
        """
        After pruning edges, update the corresponding sites in reads to "x"
        """
        reads = self.reads
        edge_per_position = self.edge_per_position
        new_reads = {}
        for read in reads:
            hap = reads[read]
            new_hap = list(hap)
            nvar = len(hap)
            for i in range(nvar - 1):
                base1 = hap[i]
                base2 = hap[i + 1]
                if "x" not in [base1, base2]:
                    found = False
                    for node1, node2, n in edge_per_position[i]:
                        if node1 == base1 + "-" + str(i) and node2 == base2 + "-" + str(
                            i + 1
                        ):
                            found = True
                    if found is False:
                        new_hap[i] = "x"
                        new_hap[i + 1] = "x"
            if new_hap.count("x") <= nvar - 2:
                new_reads.setdefault(read, "".join(new_hap))
        self.reads = new_reads

    def merge_two_pos(self, new_nodes, debug=False):
        """
        Parameters:
            new_nodes (dict): two positions and their associated partial haplotypes
        Returns:
            sub_hap_support_reads (dict): dict of (merged haps):(supporting reads)
            dnext (dict): dict of (hap1):(set of potential linking hap2s downstream)
            dbefore (dict): dict of (hap1):(set of potential linking hap2s upstream)
        """
        reads = self.reads
        sub_hap_support_reads = {}
        dnext = {}
        dbefore = {}
        for read in reads:
            hap = reads[read]
            best_matches = []
            nmatches = []
            for pos in new_nodes:
                _, best_match, max_match = self.get_matching_hap(new_nodes, pos, hap)
                nmatches.append(max_match)
                best_matches.append(best_match)
            if None not in best_matches:
                # if debug:
                #    print(read, hap, best_matches, nmatches)
                num_matches = [len(a) for a in best_matches]
                if num_matches == [1, 1]:
                    joined_haplotype = best_matches[0][0] + "" + best_matches[1][0]
                    sub_hap_support_reads.setdefault(joined_haplotype, []).append(read)
                for a in best_matches[0]:
                    for b in best_matches[1]:
                        dnext.setdefault(a, set()).add(b)
                        dbefore.setdefault(b, set()).add(a)
        return sub_hap_support_reads, dnext, dbefore

    @staticmethod
    def get_missing(subregion_haps_assembled, left_original, right_original):
        """
        find those that cannot span
        Parameters:
            subregion_haps_assembled (lst): assembled haps
            left_original (lst): left candidates
            right_original (lst): right candidates
        Returns:
            left_missing (lst): left candidates missing
            right_missing (lst): right candidates missing
        """
        left_missing = []
        right_missing = []
        for test_hap in left_original:
            found = False
            for ass_hap in subregion_haps_assembled:
                if test_hap == ass_hap[: len(test_hap)]:
                    found = True
            if found is False:
                left_missing.append(test_hap)
        for test_hap in right_original:
            found = False
            for ass_hap in subregion_haps_assembled:
                if test_hap == ass_hap[0 - len(test_hap) :]:
                    found = True
            if found is False:
                right_missing.append(test_hap)
        return left_missing, right_missing

    # get haplotypes between pos1 and pos2
    def subregion_assembly(self, pos1, pos2, allow_x=False, debug=False):
        """
        Assemble the region between pos1 and pos2
        Parameters:
            allow_x (bool): if true, missing candidates will be output with x's appended
        Returns:
            bool: assembly success
            pos1_haps (lst of str): assembled haps
            dnext (dict):
            dbefore (dict):
        """
        pos1_next = pos1
        pos1_haps = self.get_hap_by_pos(pos1)

        # extend positions one by one
        while pos1_next < pos2:
            next_pos, next_haps = self.get_next_pos(pos1_next)
            assert pos1_haps != [] and next_haps != []
            subregion_haps_assembled = set()

            if len(pos1_haps) == 1 or len(next_haps) == 1:
                for a in pos1_haps:
                    for b in next_haps:
                        subregion_haps_assembled.add(a + b)
                        dnext = {}
                        dbefore = {}
            else:
                new_nodes = {pos1: pos1_haps, next_pos: next_haps}
                if debug:
                    print("positions", pos1, next_pos)
                    print(new_nodes)
                sub_hap_support_reads, dnext, dbefore = self.merge_two_pos(
                    new_nodes, debug=debug
                )

                # do some filtering
                thres = 2  # 3
                hap_len1 = len(pos1_haps[0])
                hap_len2 = len(next_haps[0])
                hap_lens = [hap_len1, hap_len2]
                if min(hap_lens) > 2 or (min(hap_lens) == 2 and max(hap_lens) >= 10):
                    thres = 1  # 2
                for sub_hap in sub_hap_support_reads:
                    if debug:
                        print(sub_hap, len(sub_hap_support_reads[sub_hap]))
                    if len(sub_hap_support_reads[sub_hap]) >= thres:
                        subregion_haps_assembled.add(sub_hap)

            left_missing, right_missing = self.get_missing(
                subregion_haps_assembled, pos1_haps, next_haps
            )

            if right_missing != [] or left_missing != []:
                if allow_x is False:
                    return False, pos1_haps, dnext, dbefore

                if debug:
                    pprint(dnext)
                # rescue missing haps
                if (
                    hap_len1 >= 3
                    and hap_len2 >= 3
                    and (next_pos - 1) in self.dnhap
                    and self.dnhap[next_pos - 1] != -1
                    and (next_pos - 1) in self.edge_info
                    and self.edge_info[next_pos - 1] != []
                    and min(self.edge_info[next_pos - 1]) >= 2
                ):
                    if debug:
                        print("missing", left_missing, right_missing)
                    rescued_haps = self.rescue_missing(
                        left_missing, right_missing, dnext, dbefore, debug=debug
                    )
                    subregion_haps_assembled = subregion_haps_assembled.union(
                        rescued_haps
                    )

                    left_missing, right_missing = self.get_missing(
                        subregion_haps_assembled, pos1_haps, next_haps
                    )
                    if debug:
                        print("missing2", left_missing, right_missing)

                # append x's to missing haps
                if right_missing != []:
                    len_x = len(pos1_haps[0])
                    for test_hap in right_missing:
                        subregion_haps_assembled.add("x" * len_x + test_hap)
                if left_missing != []:
                    len_x = len(next_haps[0])
                    for test_hap in left_missing:
                        subregion_haps_assembled.add(test_hap + "x" * len_x)

            if debug:
                print("after merging ", subregion_haps_assembled)
            pos1_haps = list(subregion_haps_assembled)
            pos1_next = next_pos

        return True, pos1_haps, dnext, dbefore

    @staticmethod
    def rescue_missing(left_missing, right_missing, dnext, dbefore, debug=False):
        """
        Two scenarios:
        1. there is only one link and that one link is the only one in the other
        missing set, or the other missing set is empty
        2. there can be multiple links but only one of them is in the other missing
        set, and that one link is the only one in the other missing set
        """
        rescued_haps = set()
        for missing_hap in left_missing:
            if missing_hap in dnext and missing_hap.count("x") <= 1:
                possible_links = list(dnext[missing_hap])
                if len(possible_links) == 1 and right_missing in [
                    [],
                    [possible_links[0]],
                ]:
                    if debug:
                        print("merging", missing_hap, possible_links)
                    rescued_haps.add(missing_hap + possible_links[0])
                else:
                    possible_links_also_missing = [
                        a for a in possible_links if a in right_missing
                    ]
                    if (
                        len(possible_links_also_missing) == 1
                        and right_missing == possible_links_also_missing
                        and len(left_missing) == 1
                    ):
                        if debug:
                            print(
                                "merging",
                                missing_hap,
                                possible_links,
                                possible_links_also_missing,
                            )
                        rescued_haps.add(missing_hap + possible_links_also_missing[0])
        for missing_hap in right_missing:
            if missing_hap in dbefore and missing_hap.count("x") <= 1:
                possible_links = list(dbefore[missing_hap])
                if len(possible_links) == 1 and left_missing in [
                    [],
                    [possible_links[0]],
                ]:
                    if debug:
                        print("merging", possible_links, missing_hap)
                    rescued_haps.add(possible_links[0] + missing_hap)
                else:
                    possible_links_also_missing = [
                        a for a in possible_links if a in left_missing
                    ]
                    if (
                        len(possible_links_also_missing) == 1
                        and left_missing == possible_links_also_missing
                        and len(right_missing) == 1
                    ):
                        if debug:
                            print(
                                "merging",
                                possible_links,
                                possible_links_also_missing,
                                missing_hap,
                            )
                        rescued_haps.add(possible_links_also_missing[0] + missing_hap)
        return rescued_haps

    @staticmethod
    def get_matching_hap(new_nodes, pos1, hap):
        """
        Parameters:
            new_nodes (dict): two positions and their associated partial haplotypes
            pos1 (int): positions of interest
            hap (str): full read bases
        Returns:
            hap1: corresponding segment in the read
            best_match (lst of str): list of best matches
            max_match (int): length of the best matches
        """
        candidates = list(new_nodes[pos1])
        pos2 = len(candidates[0])
        hap1 = hap[pos1 : (pos1 + pos2)]
        matches = {}
        for candidate_hap in candidates:
            match, mismatch, _ = VariantGraph.compare_two_haps(candidate_hap, hap1)
            if mismatch == 0 and match > 0:
                matches.setdefault(match, []).append(candidate_hap)
        best_match = None
        max_match = 0
        if matches != {}:
            all_matches = list(matches.keys())
            max_match = max(all_matches)
            best_match = matches[max_match]
        return hap1, best_match, max_match

    @staticmethod
    def compare_two_haps(hap1, hap2):
        """
        Calculate number of matching/mismatching/additional bases between two haplotypes
        """
        if len(hap1) != len(hap2):
            raise Exception("Two haplotypes are not equal length")
        match = 0
        mismatch = 0
        extra = 0
        for i, base1 in enumerate(hap1):
            base2 = hap2[i]
            if "x" not in [base1, base2]:
                if base1 == base2:
                    match += 1
                else:
                    mismatch += 1
            if base1 == "x" and base2 != "x":
                extra += 1
        return match, mismatch, extra

    @staticmethod
    def match_reads_and_haplotypes(haplotype_per_read, hap_list, min_match=1):
        """
        Get reads that support each haplotype and possible haplotypes of each read
        """
        support_reads = {}
        support_reads_nonunique = {}
        supporting_haps_per_read = {}
        for read_name in haplotype_per_read:
            read_hap = haplotype_per_read[read_name]
            matching_haplotypes = []
            for haplotype_to_extend in hap_list:
                match, mismatch, extend = VariantGraph.compare_two_haps(
                    haplotype_to_extend, read_hap
                )
                if mismatch == 0 and match >= min_match and extend >= 0:
                    matching_haplotypes.append(haplotype_to_extend)
                    support_reads_nonunique.setdefault(haplotype_to_extend, []).append(
                        read_hap
                    )
            supporting_haps_per_read.setdefault(read_name, matching_haplotypes)
            # unique match
            if len(matching_haplotypes) == 1:
                support_reads.setdefault(matching_haplotypes[0], []).append(read_hap)
            else:
                # allow for scenarios like a read "x11x" supporting both "21xx" and "xx12"
                for hap1 in matching_haplotypes:
                    has_overlap = False
                    for hap2 in matching_haplotypes:
                        if hap1 != hap2:
                            match, mismatch, extend = VariantGraph.compare_two_haps(
                                hap1, hap2
                            )
                            if match != 0 or mismatch != 0:
                                has_overlap = True
                                break
                    # has multiple supporting haplotypes but they have no overlap at all
                    # so all are unique
                    if has_overlap is False:
                        support_reads.setdefault(hap1, []).append(read_hap)
        return VariantGraph.HaplotypesReads(
            support_reads,
            support_reads_nonunique,
            supporting_haps_per_read,
        )

    def rm_add_edges(self, pos1, pos2, new_haps):
        """
        Add new nodes/edges between two positions and remove old edges/nodes
        Updates self.nodes, self.edge_per_position, self.edges
        Also updates self.node_positions for plotting
        """
        new_pos = pos1
        old_edge_per_position = copy.deepcopy(self.edge_per_position)

        prev_pos, _ = self.get_previous_pos(pos1)
        # rm edges, from prev_pos to pos2+1
        nstart = 0
        if prev_pos is not None:
            nstart = prev_pos
        for pos in range(nstart, pos2 + 1):
            if pos in self.edge_per_position:
                for node1, node2, _ in self.edge_per_position[pos]:
                    self.rm_edge(node1, node2)

        # pop positions before and next
        if prev_pos is not None:
            self.edge_per_position.pop(prev_pos, None)
        for n in range(pos1, pos2 + 1):
            self.edge_per_position.pop(n, None)
            old_nodes = [a for a in self.nodes if int(a.split("-")[1]) == n]
            for old_node in old_nodes:
                self.nodes.remove(old_node)

        for i, new_hap in enumerate(new_haps):
            new_node = f"{new_hap}-{pos1}"
            self.nodes.append(new_node)
            if len(new_haps) > 2:
                self.node_positions.setdefault(new_node, (new_pos * 1000, i * 5 + 10))
            else:
                self.node_positions.setdefault(new_node, (new_pos * 1000, (i + 1) * 10))
            # add edges after
            if pos2 in old_edge_per_position:
                old_edges = old_edge_per_position[pos2]
                for node1, node2, _ in old_edges:
                    if new_hap[-1] == node1[0]:
                        self.add_edge(new_node, node2)
                        self.edge_per_position.setdefault(pos1, [])
                        if (new_node, node2, None) not in self.edge_per_position[pos1]:
                            self.edge_per_position[pos1].append((new_node, node2, None))
            # add edges before
            if prev_pos is not None and prev_pos in old_edge_per_position:
                old_edges = old_edge_per_position[prev_pos]
                for node1, node2, _ in old_edges:
                    if new_hap[0] == node2[0]:
                        self.add_edge(node1, new_node)
                        self.edge_per_position.setdefault(prev_pos, [])
                        if (node1, new_node, None) not in self.edge_per_position[
                            prev_pos
                        ]:
                            self.edge_per_position[prev_pos].append(
                                (node1, new_node, None)
                            )

    def merge_edges_simple(self, pos1, pos2):
        """Merge edges when there are two haplotypes only"""
        new_haps = self.path(pos1, pos2)
        self.rm_add_edges(pos1, pos2, new_haps)

    def merge_edges(self, pos1, pos2, debug=False):
        """Merge edges in more complex scenario - more than two haplotypes"""
        success, new_haps, _, _ = self.subregion_assembly(pos1, pos2, debug=debug)
        if success is True:
            self.rm_add_edges(pos1, pos2, new_haps)
        return new_haps

    def initialize_graph(self, debug=False):
        """Initialize graph, add and prune edges"""
        self.get_edges(debug=debug)
        self.update_reads()

        for pos in self.edge_per_position:
            for node1, node2, _ in self.edge_per_position[pos]:
                self.add_edge(node1, node2)
                if node1 not in self.nodes:
                    self.nodes.append(node1)
                if node2 not in self.nodes:
                    self.nodes.append(node2)

        for i in range(5):
            for j in range(self.nvar):
                node = f"{i}-{j}"
                if node not in self.nodes:
                    hap = "x" * j + str(i) + "x" * (self.nvar - j - 1)
                    read_support = self.match_reads_and_haplotypes(
                        self.reads_original, [hap]
                    )
                    if hap in read_support.unique and len(read_support.unique[hap]) > 1:
                        self.nodes.append(node)

        for node1, node2 in self.edges:
            for n in (node1, node2):
                if n not in self.nodes:
                    self.nodes.append(n)
        self.get_positions()
        if debug:
            pprint(self.edge_per_position)

    def filter_candidates(self, hap_block, candidates, lenx, candidate_first=False):
        """
        Filter candidates against reads
        """
        read_support = {}
        reads = self.reads
        hap1 = hap_block
        for candidate in candidates:
            candidate_len = len(candidate)
            if candidate_first:
                hap2 = (
                    "x" * (lenx - candidate_len) + candidate + "x" * (self.nvar - lenx)
                )
            else:
                hap2 = (
                    "x" * (self.nvar - lenx) + candidate + "x" * (lenx - candidate_len)
                )
            for read in reads:
                read_seq = reads[read]
                match1, mismatch1, _ = VariantGraph.compare_two_haps(hap1, read_seq)
                match2, mismatch2, _ = VariantGraph.compare_two_haps(hap2, read_seq)
                if mismatch1 == 0 and mismatch2 == 0:
                    if match1 > 0 and match2 > 0:
                        read_support.setdefault(candidate, []).append(
                            min(match1, match2)
                        )
        if read_support == {}:
            return set()
        if len(read_support) == 1:
            return [list(read_support.keys())[0]]
        medians = []
        maxs = []
        for candidate in read_support:
            medians.append(np.median(read_support[candidate]))
            maxs.append(max(read_support[candidate]))
        for i in range(len(medians)):
            others = []
            for j in range(len(medians)):
                if j != i:
                    others.append(maxs[j])
            if medians[i] > max(others) and medians[i] >= 4:
                new_candidate = list(read_support.keys())[i]
                if len(read_support[new_candidate]) >= 4:
                    return [new_candidate]
        return candidates

    def extend_pivot_blocks(self, main_haps, debug=False):
        """
        Take candidate matching haplotypes and extend with shared bases
        """
        new_main_haps = []
        for hap in main_haps:
            identical_bases_right = ""
            identical_bases_left = ""
            prefix_x = 0
            suffix_x = self.nvar - 1
            if "x" in hap:
                for prefix_x in range(self.nvar):
                    if hap[prefix_x] != "x":
                        break
                for suffix_x in reversed(range(self.nvar)):
                    if hap[suffix_x] != "x":
                        break
                if debug:
                    print(f"extending {hap} {prefix_x} {suffix_x}")
                # extend right
                if suffix_x < self.nvar - 2:
                    next_pos = suffix_x + 1
                    # check that this is not a site very distant from next site
                    if (
                        suffix_x in self.dnhap
                        and self.dnhap[suffix_x] != -1
                        and suffix_x in self.edge_info
                        and self.edge_info[suffix_x] != []
                        and min(self.edge_info[suffix_x]) >= 2
                    ):
                        previous_pos, _ = self.get_previous_pos(next_pos)
                        _, _, dnext, dbefore = self.subregion_assembly(
                            previous_pos, next_pos, allow_x=False, debug=debug
                        )
                        if debug:
                            pprint(dnext)
                        hap_block = hap[previous_pos : suffix_x + 1]
                        if hap_block in dnext:  # and len(hap_block) >= 5:
                            candidates = dnext[hap_block]
                            # check candidates against reads
                            if len(candidates) > 1:
                                candidates = self.filter_candidates(
                                    hap, candidates, self.nvar - 1 - suffix_x
                                )
                            if candidates != set():
                                next_hap_len = len(list(candidates)[0])
                                i = 0
                                while True:
                                    if i >= next_hap_len:
                                        break
                                    bases_per_candidate = list(
                                        set([test_hap[i] for test_hap in candidates])
                                    )
                                    if len(bases_per_candidate) != 1:
                                        break
                                    identical_bases_right += bases_per_candidate[0]
                                    i += 1
                # extend left
                if prefix_x >= 2:
                    next_pos = prefix_x
                    _, nend = self.get_nstart_nend()
                    if next_pos < nend:
                        next_next_pos, _ = self.get_next_pos(next_pos)
                        hap_block = hap[prefix_x:next_next_pos]
                    else:
                        hap_block = hap[prefix_x:]
                    # check that this is not a site very distant from previous site
                    if (
                        (prefix_x - 1) in self.dnhap
                        and self.dnhap[prefix_x - 1] != -1
                        and (prefix_x - 1) in self.edge_info
                        and self.edge_info[prefix_x - 1] != []
                        and min(self.edge_info[prefix_x - 1]) >= 2
                    ):
                        previous_pos, _ = self.get_previous_pos(next_pos)
                        _, _, dnext, dbefore = self.subregion_assembly(
                            previous_pos, next_pos, allow_x=False, debug=debug
                        )
                        if debug:
                            pprint(dbefore)
                        if hap_block in dbefore:  # and len(hap_block) >= 5:
                            candidates = dbefore[hap_block]
                            # check candidates against reads
                            if len(candidates) > 1:
                                candidates = self.filter_candidates(
                                    hap,
                                    candidates,
                                    prefix_x,
                                    candidate_first=True,
                                )
                            if candidates != set():
                                previous_hap_len = len(list(candidates)[0])
                                i = previous_hap_len - 1
                                while True:
                                    if i < 0:
                                        break
                                    bases_per_candidate = list(
                                        set([test_hap[i] for test_hap in candidates])
                                    )
                                    if len(bases_per_candidate) != 1:
                                        break
                                    identical_bases_left += bases_per_candidate[0]
                                    i -= 1
                                identical_bases_left = identical_bases_left[::-1]
            extended_hap = (
                "x" * (prefix_x - len(identical_bases_left))
                + identical_bases_left
                + hap[prefix_x : suffix_x + 1]
                + identical_bases_right
                + "x" * (self.nvar - 1 - suffix_x - len(identical_bases_right))
            )
            if debug and extended_hap != hap:
                print(f"after extension is {extended_hap}")
            new_main_haps.append(extended_hap)
        return new_main_haps

    def assemble_haplotypes(self, debug=False, make_plot=True):
        """Stepwise assembly of haplotype blocks based on graph"""
        if make_plot:
            _, axs = plt.subplots(5)
            self.visualize(axs[0])

        self.summerize_edges()
        segments = self.get_segments()
        if debug:
            print(segments)

        # two haplotypes only, simple walk through the graph
        if len(segments) == 1 and list(segments.values())[0] == 2:
            segment_start, segment_end = list(segments.keys())[0]
            assembled_haps = self.path(segment_start, segment_end + 1)
            if make_plot:
                self.plotfig()
            return assembled_haps

        # identify simple parts of the graph with just two haplotypes
        for i, (segment_start, segment_end) in enumerate(segments):
            if segments[(segment_start, segment_end)] == 2:
                if i not in [0, len(segments) - 1]:
                    segment_start += 1
                if segment_end == self.nvar - 2:
                    segment_end += 1
                if segment_end > segment_start:
                    self.merge_edges_simple(
                        segment_start,
                        segment_end,
                    )
        if make_plot:
            self.visualize(axs[1])
        # identify more complex parts of the graph
        for i, (segment_start, segment_end) in enumerate(segments):
            if (
                segments[(segment_start, segment_end)] > 2
                and segment_end >= segment_start
            ):
                segment_end += 1
                _ = self.merge_edges(
                    segment_start,
                    segment_end,
                    debug=False,
                )
                if debug:
                    print(segment_start, segment_end, 10)
                    pprint(self.edge_per_position)
        if make_plot:
            self.visualize(axs[2])

        if debug:
            print("done with segments")

        nstart, nend = self.get_nstart_nend()
        if nstart == nend:
            return [a.split("-")[0] for a in self.nodes]

        # scan whole graph and merge when possible

        # generalize this special logic for future
        # if first site is far away from second, start from the second site
        # if 1 in self.edge_per_position and (
        #    self.edge_info[0] == [] or min(self.edge_info[0]) < 5
        # ):
        #    nstart = 1
        while nstart < nend - 1:
            if debug:
                print(f"scanning {nstart}-{nend}")
            success, final_haps, _, _ = self.subregion_assembly(
                nstart, nend, debug=debug
            )
            current_haps = self.get_hap_by_pos(nstart)
            pos = nstart + len(final_haps[0])
            block_ending_pos, _ = self.get_previous_pos(pos)
            if len(current_haps[0]) != len(final_haps[0]):
                self.rm_add_edges(nstart, block_ending_pos, final_haps)
            nstart = pos
            if success:
                break

        if make_plot:
            self.visualize(axs[3])

        if debug:
            print("done with scanning...")
            pprint(self.edge_per_position)

        # merge blocks by size
        self.merge_blocks_by_size(debug=debug)
        if make_plot:
            self.visualize(axs[4])
            self.plotfig()

        # final merge
        nstart, nend = self.get_nstart_nend()
        if nstart == nend:
            return [a.split("-")[0] for a in self.nodes]
        if debug:
            print(f"running final {nstart}-{nend}")
            pprint(self.nodes)
        success, final_haps, _, _ = self.subregion_assembly(
            nstart, nend, allow_x=True, debug=debug
        )

        return final_haps

    def order_blocks_by_size(self):
        """When merging blocks, longer ones will be done first"""
        hap_lens = set()
        for node in self.nodes:
            hap, pos = node.split("-")
            pos = int(pos)
            hap_len = len(hap)
            hap_lens.add((pos, hap_len))
        hap_lens = sorted(hap_lens, key=lambda x: x[0])
        block_order = []
        for i, node in enumerate(hap_lens):
            if i < len(hap_lens) - 1:
                len1 = node[1]
                next_node = hap_lens[i + 1]
                len2 = next_node[1]
                block_order.append(
                    (node[0], next_node[0], min(len1, len2), max(len1, len2))
                )
        block_order = sorted(block_order, key=lambda x: (0 - x[2], 0 - x[3]))
        return block_order

    def merge_blocks_by_size(self, debug=False):
        """Merge blocks, ordered by size"""
        while self.edge_per_position != {}:
            block_order = self.order_blocks_by_size()
            if len(block_order) == 1:
                break
            if debug:
                print("priority", block_order)
            i = 0
            for pos1, pos2, _, _ in block_order:
                i += 1
                success, new_haps, _, _ = self.subregion_assembly(
                    pos1, pos2, debug=debug
                )
                if success is True:
                    if debug:
                        print(f"merging blocks {pos1}-{pos2}")
                    self.rm_add_edges(pos1, pos2, new_haps)
                    break
            if i == len(block_order) and success is False:
                break

    def run(self, debug=False, make_plot=False):
        """Run the whole haplotype phasing process"""
        self.initialize_graph(debug=debug)
        final_haps = self.assemble_haplotypes(debug=debug, make_plot=make_plot)
        if debug:
            print("assembled_haps are", final_haps)

        final_haps = self.filter_low_support_haps(final_haps, min_count=4, debug=debug)
        highest_cn, main_haps_candidates = self.get_highest_cn(final_haps)

        # extend all possible haps
        extended_main_hap_candidates = []
        for main_haps_candidate in main_haps_candidates:
            extended = self.extend_pivot_blocks(main_haps_candidate, debug=debug)
            extended_main_hap_candidates.append(extended)
        if self.pivot_pos != -1:
            scores = []
            for extended_main_hap_candidate in extended_main_hap_candidates:
                pivot_ones = [
                    a for a in extended_main_hap_candidate if a[self.pivot_pos] != "x"
                ]
                nhap = len(pivot_ones)
                total_length = sum([len(a.strip("x")) for a in pivot_ones])

                represented_haps = []
                for hap1 in final_haps:
                    mismatches = []
                    for hap2 in pivot_ones:
                        match, mismatch, extend = VariantGraph.compare_two_haps(
                            hap1, hap2
                        )
                        mismatches.append(mismatch)
                    if 0 in mismatches:
                        represented_haps.append(hap1)
                num_represented_haps = len(represented_haps)

                scores.append(((nhap, num_represented_haps, total_length), pivot_ones))
            sort_candidates = sorted(
                scores,
                key=lambda x: x[0],
                reverse=True,
            )
            main_haps = sort_candidates[0][1]
        else:
            main_haps = extended_main_hap_candidates[0]
        if main_haps == []:
            main_haps = final_haps
        return main_haps, final_haps, highest_cn
