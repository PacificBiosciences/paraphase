# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


import pysam
import os
import copy
import numpy as np
from collections import Counter
from itertools import product
import re
import logging
import json
import traceback
from scipy.stats import poisson
from collections import namedtuple
from .haplotype_assembler import VariantGraph


class Phaser:
    clip_5p = r"^\d+S|^\d+H"
    clip_3p = r"\d+S$|\d+H$"
    deletion = r"\d+D"
    fields = [
        "total_cn",
        "gene_cn",
        "final_haplotypes",
        "two_copy_haplotypes",
        "alleles_final",
        "haplotype_links",
        "highest_total_cn",
        "assembled_haplotypes",
        "sites_for_phasing",
        "unique_supporting_reads",
        "het_sites_not_used_in_phasing",
        "homozygous_sites",
        "haplotype_details",
        "nonunique_supporting_reads",
        "read_details",
        "genome_depth",
        "region_depth",
        "sample_sex",
        "heterozygous_sites",
        "phase_region",
        "raw_alleles",
        "fusions_called",
    ]
    GeneCall = namedtuple(
        "GeneCall",
        fields,
        defaults=(None,) * len(fields),
    )
    CoverageSummary = namedtuple("CoverageSummary", "median percentile80")
    MEAN_BASE_QUAL = 25

    def __init__(
        self,
        sample_id,
        outdir,
        args=None,
        wgs_depth=None,
        genome_bam=None,
        sample_sex=None,
    ):
        self.outdir = outdir
        self.sample_id = sample_id
        self.homopolymer_sites = {}
        self.het_sites = []  # for phasing
        self.init_het_sites = []
        self.het_no_phasing = []
        self.homo_sites = []
        self.candidate_pos = set()
        self.mdepth = wgs_depth
        self.genome_bam = genome_bam
        self.sample_sex = sample_sex
        self.trusted_read_support = 20
        self.targeted = False
        self.genome_build = args.genome if args is not None else None
        if args is not None:
            self.targeted = args.targeted
        self.min_vaf = None
        if args is not None:
            self.min_vaf = args.min_variant_frequency
        self.min_haplotype_frequency = None
        if args is not None:
            self.min_haplotype_frequency = args.min_haplotype_frequency

    def set_parameter(self, config):
        self.gene = config["gene"]
        self.bam = os.path.join(
            self.outdir, self.sample_id + f"_{self.gene}_realigned.bam"
        )
        if os.path.exists(self.bam) is False:
            raise Exception(f"File {self.bam} not found.")
        self._bamh = pysam.AlignmentFile(self.bam, "rb")
        self.nchr = config["nchr"]
        self.ref = config["data"]["reference"]
        self._refh = pysam.FastaFile(self.ref)
        self.pivot_site = None
        if "pivot_site" in config:
            self.pivot_site = config["pivot_site"]
        self.nchr_old = config["nchr_old"]
        self.offset = int(self.nchr_old.split("_")[1]) - 1

        # define region boundary
        self.left_boundary = config.get("left_boundary")
        self.right_boundary = config.get("right_boundary")
        if self.left_boundary is None:
            self.left_boundary = int(self.nchr_old.split("_")[1])
        if self.right_boundary is None:
            self.right_boundary = int(self.nchr_old.split("_")[2])
        self.gene_start = config.get("gene_start")
        if self.gene_start is None:
            self.gene_start = self.left_boundary
        self.gene_end = config.get("gene_end")
        if self.gene_end is None:
            self.gene_end = self.right_boundary

        self.call_fusion = None
        if "call_fusion" in config:
            self.call_fusion = config["call_fusion"]
            fusion_json = config["data"].get("fusion_json")
            self.fusion_gene_def_variants = []
            with open(fusion_json) as f:
                self.fusion_gene_def_variants = json.load(f).get(self.gene)
        self.use_supplementary = False
        if "use_supplementary" in config or "is_tandem" in config:
            self.use_supplementary = True
        self.to_phase = False
        if "to_phase" in config or "is_tandem" in config:
            self.to_phase = True
        self.is_reverse = False
        if "is_reverse" in config:
            self.is_reverse = config["is_reverse"]
        self.is_palindrome = False
        if "is_palindrome" in config:
            self.is_palindrome = config["is_palindrome"]
        self.expect_cn2 = False
        if "expect_cn2" in config:
            self.expect_cn2 = True
        self.clip_3p_positions = []
        self.clip_5p_positions = []
        if "clip_3p_positions" in config:
            self.clip_3p_positions = config["clip_3p_positions"]
        self.clip_3p_positions = sorted(self.clip_3p_positions)
        if "clip_5p_positions" in config:
            self.clip_5p_positions = config["clip_5p_positions"]
        self.clip_5p_positions = sorted(self.clip_5p_positions)
        self.noisy_region = []
        if "noisy_region" in config:
            self.noisy_region = config["noisy_region"]
        self.add_sites = []
        if "add_sites" in config:
            self.add_sites = config["add_sites"]
        self.gene1_cn2 = False
        if "gene1_cn2" in config:
            self.gene1_cn2 = True
        self.match = {}
        self.gene2_region = config.get("gene2_region")
        if self.gene2_region is not None:
            self.position_match = config["data"].get("gene_position_match")
            if self.position_match is not None:
                with open(self.position_match) as f:
                    for line in f:
                        at = line.split()
                        self.match.setdefault(int(at[0]), int(at[1]))

        self.del1_reads = set()
        self.del1_reads_partial = set()
        self.del2_reads = set()
        self.del2_reads_partial = set()
        self.deletion1_size = config.get("deletion1_size")
        self.deletion2_size = config.get("deletion2_size")
        self.deletion1_name = config.get("deletion1_name")
        self.deletion2_name = config.get("deletion2_name")
        self.del2_3p_pos1 = config.get("del2_3p_pos1")
        self.del2_3p_pos2 = config.get("del2_3p_pos2")
        self.del2_5p_pos1 = config.get("del2_5p_pos1")
        self.del2_5p_pos2 = config.get("del2_5p_pos2")
        self.del1_3p_pos1 = config.get("del1_3p_pos1")
        self.del1_3p_pos2 = config.get("del1_3p_pos2")
        self.del1_5p_pos1 = config.get("del1_5p_pos1")
        self.del1_5p_pos2 = config.get("del1_5p_pos2")

    def get_regional_depth(self, bam_handle, query_region, ninterval=100):
        """Get depth of the query regions"""
        region_depth = []
        for region in query_region:
            depth = []
            nstep = max(1, int((region[1] - region[0]) / ninterval))
            for pos in range(region[0], region[1], nstep):
                site_depth = bam_handle.count(
                    self.nchr, pos - 1, pos, read_callback="all"
                )
                depth.append(site_depth)
            region_depth.append(
                self.CoverageSummary(np.median(depth), np.nanpercentile(depth, 80))
            )
        return region_depth

    def check_coverage_before_analysis(self):
        """check low coverage regions for enrichment data"""
        check_region_start = self.left_boundary
        check_region_end = self.right_boundary
        if self.clip_5p_positions != [] and self.clip_3p_positions != []:
            clip_5p_max = max(self.clip_5p_positions)
            clip_3p_min = min(self.clip_3p_positions)
            if clip_5p_max < clip_3p_min:
                check_region_start = clip_5p_max
                check_region_end = clip_3p_min
        elif self.clip_5p_positions != []:
            check_region_start = max(self.left_boundary, max(self.clip_5p_positions))
        elif self.clip_3p_positions != []:
            check_region_end = min(self.right_boundary, min(self.clip_3p_positions))
        depth1 = self.get_regional_depth(
            self._bamh, [[self.left_boundary, self.right_boundary]]
        )[0]
        self.region_avg_depth = depth1
        if (
            check_region_start != self.left_boundary
            or check_region_end != self.right_boundary
        ):
            depth2 = self.get_regional_depth(
                self._bamh, [[check_region_start, check_region_end]]
            )[0]
            if np.isnan(depth1.median) or depth2.median > depth1.median:
                self.region_avg_depth = depth2
        if np.isnan(self.region_avg_depth.median) or (
            self.region_avg_depth.median <= 8
            and self.region_avg_depth.percentile80 < 50
        ):
            logging.warning(
                "This region does not appear to have coverage. Will not attempt to phase haplotypes."
            )
            return False
        return True

    def get_range_in_other_gene(self, pos, search_range=200):
        """
        Find the correponding coordinates in the other gene
        """
        if pos in self.match:
            return self.match[pos]
        for i in range(search_range):
            new_pos = pos + i
            if new_pos in self.match:
                return self.match[new_pos]
        return None

    def get_homopolymer(self):
        """
        Get the homopolymer and dinucleotide sites
        Variants for phasing:
        SNVs, not matching excluded bases and sites labeled as 1
        indels, not in excluded sites
        """
        seq = self._refh.fetch(self.nchr_old).upper()
        nstart = self.offset
        drepeat = {}
        # start with 5mer
        for i in range(len(seq) - 5):
            for nu in ["A", "C", "G", "T"]:
                if seq[i : i + 5].count(nu) >= 5:
                    for pos in range(i, i + 5):
                        drepeat.setdefault(pos, nu)
        # add nearby 4mers (same base), separated by stretches of a different base
        while True:
            extend = False
            for i in range(len(seq) - 4):
                for nu in ["A", "C", "G", "T"]:
                    if seq[i : i + 4].count(nu) == 4:
                        for j in range(1, 5):
                            if (
                                i + 4 + j in drepeat
                                and drepeat[i + 4 + j] == nu
                                and len(set(seq[i + 4 : i + 4 + j])) == 1
                            ):
                                for pos in range(i, i + 4 + j):
                                    if pos not in drepeat:
                                        drepeat.setdefault(pos, nu)
                                        extend = True
                            if (
                                i - 1 - j in drepeat
                                and drepeat[i - 1 - j] == nu
                                and len(set(seq[i - j : i])) == 1
                            ):
                                for pos in range(i - j, i + 4):
                                    if pos not in drepeat:
                                        drepeat.setdefault(pos, nu)
                                        extend = True
            if extend is False:
                break
        drepeat = dict(sorted(drepeat.items()))
        drepeat_keys = list(drepeat.keys())
        for i in range(len(drepeat)):
            pos = drepeat_keys[i]
            pos_adjusted = pos + nstart + 1
            base = drepeat[pos]
            if i + 1 == len(drepeat) or (
                i + 1 < len(drepeat) and drepeat_keys[i + 1] != pos + 1
            ):
                self.homopolymer_sites.setdefault(pos_adjusted, base)
                self.homopolymer_sites.setdefault(pos + 1 + nstart + 1, f"{base},1")
            elif i == 0 or (i > 0 and drepeat_keys[i - 1] != pos - 1):
                self.homopolymer_sites.setdefault(pos - 1 + nstart + 1, base)
                self.homopolymer_sites.setdefault(pos_adjusted, "A,C,G,T")
            else:
                self.homopolymer_sites.setdefault(pos_adjusted, base)

        # dinucleotide
        drepeat = {}
        # start with 4 copies of a dimer
        for i in range(len(seq) - 8):
            for dimer in [a + b for a, b in product("ACGT", repeat=2) if a != b]:
                if seq[i : i + 8].count(dimer) == 4:
                    for pos in range(i, i + 8):
                        drepeat.setdefault(pos, []).append(dimer)
        # add nearby any 2 copies of the same dimer, separated by either stretches of
        # the same 2 bases or a different dimer
        counter = 0
        while counter < 3:
            counter += 1
            extend = False
            for i in range(len(seq) - 4):
                for dimer in [a + b for a, b in product("ACGT", repeat=2) if a != b]:
                    if seq[i : i + 4].count(dimer) == 2:
                        for j in range(8):
                            pos_to_check = i + 4 + j
                            to_add = False
                            if pos_to_check in drepeat and (
                                dimer in drepeat[pos_to_check]
                                or dimer[::-1] in drepeat[pos_to_check]
                            ):
                                if j <= 1:
                                    to_add = True
                                else:
                                    mid_seq = seq[i + 4 : i + 4 + j]
                                    reduced_nu = list(set(mid_seq))
                                    if len(reduced_nu) == 2:
                                        base1 = reduced_nu[0]
                                        base2 = reduced_nu[1]
                                        if (
                                            (
                                                base1 + base1 not in mid_seq
                                                and base2 + base2 not in mid_seq
                                            )
                                            or base1 + base2 == dimer
                                            or base2 + base1 == dimer
                                        ):
                                            to_add = True
                                if to_add is True:
                                    for pos in range(i, i + 4 + j):
                                        drepeat.setdefault(pos, [])
                                        if dimer not in drepeat[pos]:
                                            drepeat[pos].append(dimer)
                                        extend = True
                            pos_to_check = i - 1 - j
                            to_add = False
                            if pos_to_check in drepeat and (
                                dimer in drepeat[pos_to_check]
                                or dimer[::-1] in drepeat[pos_to_check]
                            ):
                                if j <= 1:
                                    to_add = True
                                else:
                                    mid_seq = seq[i - j : i]
                                    reduced_nu = list(set(mid_seq))
                                    if len(reduced_nu) == 2:
                                        base1 = reduced_nu[0]
                                        base2 = reduced_nu[1]
                                        if (
                                            (
                                                base1 + base1 not in mid_seq
                                                and base2 + base2 not in mid_seq
                                            )
                                            or base1 + base2 == dimer
                                            or base2 + base1 == dimer
                                        ):
                                            to_add = True
                                if to_add is True:
                                    for pos in range(i - j, i + 4):
                                        drepeat.setdefault(pos, [])
                                        if dimer not in drepeat[pos]:
                                            drepeat[pos].append(dimer)
                                        extend = True
            if extend is False:
                break
        drepeat = dict(sorted(drepeat.items()))
        for pos in drepeat:
            pos_adjusted = pos + nstart + 1
            if pos_adjusted not in self.homopolymer_sites:
                bases = []
                for nus in drepeat[pos]:
                    bases.append(nus[0])
                    bases.append(nus[1])
                bases = list(set(bases))  # + ["1"]
                self.homopolymer_sites.setdefault(pos_adjusted, ",".join(bases))

    @staticmethod
    def depth_prob(nread, haploid_depth):
        """Find probability of cn state based on depth"""
        prob = []
        for i in range(4):
            depthexpected = (i + 1) * haploid_depth
            pmf = poisson.pmf(nread, depthexpected)
            prob.append(pmf)
        sum_prob = sum(prob)
        if sum_prob == 0:
            return None
        post_prob = [float(a) / float(sum_prob) for a in prob]
        return post_prob

    def find_big_deletion(self, min_size=5000, min_count=3, padding=50):
        """Call big deletions from data"""
        del_reads = []
        for read in self._bamh.fetch(
            self.nchr, self.left_boundary, self.right_boundary
        ):
            read_cigar = read.cigarstring
            del_len = [int(a[:-1]) for a in re.findall(Phaser.deletion, read_cigar)]
            if del_len != []:
                longest_del_size = max(del_len)
                if longest_del_size >= min_size:
                    cigar_before_del = read_cigar.split(f"{longest_del_size}D")[0]
                    del_pos = read.reference_start + sum(
                        [
                            int(a[:-1])
                            for a in re.findall(r"\d+D|\d+M", cigar_before_del)
                        ]
                    )
                    del_reads.append((del_pos, del_pos + longest_del_size))
        counter = Counter(del_reads).most_common(2)
        if len(counter) > 0 and counter[0][1] >= min_count:
            nstart, nend = counter[0][0]
            self.deletion1_size = nend - nstart
            self.deletion1_name = f"{nstart}_del_{nend-nstart}"
            self.del1_3p_pos1 = nstart - padding
            self.del1_3p_pos2 = nstart + padding
            self.del1_5p_pos1 = nend - padding
            self.del1_5p_pos2 = nend + padding
        if len(counter) > 1 and counter[1][1] >= min_count:
            nstart, nend = counter[1][0]
            self.deletion2_size = nend - nstart
            self.deletion2_name = f"{nstart}_del_{nend-nstart}"
            self.del2_3p_pos1 = nstart - padding
            self.del2_3p_pos2 = nstart + padding
            self.del2_5p_pos1 = nend - padding
            self.del2_5p_pos2 = nend + padding

    @staticmethod
    def check_del(read, del_size, pos1, pos2):
        """Find reads having the 6.3kb deletion in its cigar string"""
        starting_pos = read.reference_start
        for (op_type, op_len) in read.cigartuples:
            if op_type == 2 and abs(op_len - del_size) < min(del_size * 0.1, 50):
                if pos1 - 3 <= starting_pos <= pos2 + 3:
                    return True
            if op_type in [0, 2, 7, 8]:
                starting_pos += op_len
        return False

    def get_long_del_reads(
        self,
        p3_pos1,
        p3_pos2,
        p5_pos1,
        p5_pos2,
        del_size,
        min_clip_len=300,
        min_extend=1000,
    ):
        """
        Find reads having big deletions. (Improve to general SVs for future)
        Could be softclipped at either side or have the deletion in cigar.
        Parameters:
            min_clip_len (int): minimum length for the soft-clip
        Returns: fully-spanning reads (set), partially spanning reads (set)
        """
        bamh = self._bamh
        p5_reads = set()
        p3_reads = set()
        del_reads = set()
        # 3 prime clip
        pos1 = p3_pos1
        pos2 = p3_pos2
        reference_start_cutoff = pos1 - min_extend
        for read in bamh.fetch(self.nchr, pos1, pos2):
            read_name = self.get_read_name(read)
            find_clip_3p = re.findall(self.clip_3p, read.cigarstring)
            if find_clip_3p != [] and pos1 < read.reference_end < pos2:
                if (
                    int(find_clip_3p[0][:-1]) >= min_clip_len
                    and read.reference_start < reference_start_cutoff
                ):
                    p3_reads.add(read_name)
            deletion_in_cigar = self.check_del(read, del_size, pos1, pos2)
            if deletion_in_cigar:
                del_reads.add(read_name)
        # 5 prime clip
        pos1 = p5_pos1
        pos2 = p5_pos2
        reference_end_cutoff = pos2 + min_extend
        for read in bamh.fetch(self.nchr, pos1, pos2):
            read_name = self.get_read_name(read)
            find_clip_5p = re.findall(self.clip_5p, read.cigarstring)
            if find_clip_5p != [] and pos1 < read.reference_start < pos2:
                if (
                    int(find_clip_5p[0][:-1]) >= min_clip_len
                    and read.reference_end > reference_end_cutoff
                ):
                    p5_reads.add(read_name)
        if del_reads != set() or (p3_reads != set() and p5_reads != set()):
            return (
                del_reads.union(p3_reads.intersection(p5_reads)),
                del_reads.union(p3_reads).union(p5_reads),
            )
        return set(), set()

    def get_pivot_site_index(self):
        """Return the index of the pivot site in list of het sites"""
        positions = [int(site.split("_")[0]) for site in self.het_sites]
        if self.pivot_site in positions:
            return positions.index(self.pivot_site), True
        return -1, False

    def get_read_name(self, read):
        """Rename reads when supplementary"""
        read_name = read.query_name
        if read.is_supplementary and self.use_supplementary:
            read_name = (
                read_name + f"_sup_{read.reference_start}_{read.reference_length}"
            )
        return read_name

    def get_read_names(self, read, partial_deletion_reads):
        """Add read names for supplementary alignments"""
        read_names = [read.query_name]
        if read.is_supplementary and self.use_supplementary:
            sup_name = (
                read.query_name + f"_sup_{read.reference_start}_{read.reference_length}"
            )
            read_names = [sup_name]
            if (
                sup_name in partial_deletion_reads
                and read.query_name in partial_deletion_reads
            ):
                read_names.append(read.query_name)
        return read_names

    def add_variant_around_clip(self):
        """add sites before 5p clip and after 3p clip"""
        if self.clip_5p_positions != []:  # and self.het_sites != []:
            for i in range(len(self.clip_5p_positions)):
                clip_pos = self.clip_5p_positions[i]
                if i > 0:
                    prev_clip_pos = self.clip_5p_positions[i - 1]
                    var_before_clip = [
                        a
                        for a in self.het_sites
                        if prev_clip_pos < int(a.split("_")[0]) < clip_pos
                    ]
                else:
                    var_before_clip = [
                        a for a in self.het_sites if int(a.split("_")[0]) < clip_pos
                    ]
                if var_before_clip == []:
                    var_pos = clip_pos - 30
                    ref_base = self._refh.fetch(
                        self.nchr_old, var_pos - self.offset - 1, var_pos - self.offset
                    ).upper()
                    var_base = [a for a in ["A", "C", "G", "T"] if a != ref_base][0]
                    new_var = f"{var_pos}_{ref_base}_{var_base}"
                    self.het_sites.append(new_var)
        if self.clip_3p_positions != []:  # and self.het_sites != []:
            for i in reversed(range(len(self.clip_3p_positions))):
                clip_pos = self.clip_3p_positions[i]
                if i == len(self.clip_3p_positions) - 1:
                    var_after_clip = [
                        a for a in self.het_sites if int(a.split("_")[0]) > clip_pos
                    ]
                else:
                    next_clip_pos = self.clip_3p_positions[i + 1]
                    var_after_clip = [
                        a
                        for a in self.het_sites
                        if clip_pos < int(a.split("_")[0]) < next_clip_pos
                    ]
                if var_after_clip == []:
                    var_pos = clip_pos + 30
                    ref_base = self._refh.fetch(
                        self.nchr_old, var_pos - self.offset - 1, var_pos - self.offset
                    ).upper()
                    var_base = [a for a in ["A", "C", "G", "T"] if a != ref_base][0]
                    new_var = f"{var_pos}_{ref_base}_{var_base}"
                    self.het_sites.append(new_var)

    def get_haplotypes_from_reads(
        self,
        exclude_reads=[],
        min_mapq=5,
        min_clip_len=50,
        clip_buffer=20,
        check_clip=False,
        partial_deletion_reads=[],
        kept_sites=[],
        add_sites=[],
        homo_sites=[],
        multi_allelic_sites={},
    ):
        """
        Go through reads and get bases at sites of interest.
        Two rounds, with variant site filtering in between.
        """
        raw_read_haps = self.get_haplotypes_from_reads_step(
            exclude_reads,
            min_mapq,
            min_clip_len,
            clip_buffer,
            check_clip,
            partial_deletion_reads,
            multi_allelic_sites,
        )
        absent_base_per_site = self.remove_var(raw_read_haps, kept_sites, homo_sites)
        if self.het_sites != []:
            for var in add_sites:
                if var not in self.het_sites:
                    self.het_sites.append(var)
        # add sites before 5p clip and after 3p clip
        self.add_variant_around_clip()
        self.het_sites = sorted(self.het_sites)
        raw_read_haps = self.get_haplotypes_from_reads_step(
            exclude_reads,
            min_mapq,
            min_clip_len,
            clip_buffer,
            check_clip,
            partial_deletion_reads,
            multi_allelic_sites,
            absent_base_per_site,
        )
        return raw_read_haps

    def get_haplotypes_from_reads_step(
        self,
        exclude_reads=[],
        min_mapq=5,
        min_clip_len=50,
        clip_buffer=20,
        check_clip=False,
        partial_deletion_reads=[],
        multi_allelic_sites={},
        absent_base_per_site={},
    ):
        """
        Go through reads and get bases at sites of interest
        Returns:
            read_haps (dict of str:list): collapse each read into just the positions
            of interest. 1 corresponds to ref, 2 corresponds to alt
        """
        het_sites = self.het_sites
        read_haps = {}
        nvar = len(het_sites)
        for dsnp_index, allele_site in enumerate(het_sites):
            snp_position_gene1, allele1, allele2, *at = allele_site.split("_")
            snp_position = int(snp_position_gene1)
            reads_with_flanking_indels = []
            for pileupcolumn in self._bamh.pileup(
                self.nchr,
                snp_position - 2,
                snp_position,
                truncate=True,
                min_base_quality=self.MEAN_BASE_QUAL,
            ):
                # require that the base on the read is not flanked by any indels
                if pileupcolumn.reference_pos == snp_position - 2:
                    for read in pileupcolumn.pileups:
                        if read.indel != 0 or read.is_del:
                            read_names = self.get_read_names(
                                read.alignment, partial_deletion_reads
                            )
                            for read_name in read_names:
                                reads_with_flanking_indels.append(read_name)
                if pileupcolumn.reference_pos == snp_position - 1:
                    for read in pileupcolumn.pileups:
                        read_names = self.get_read_names(
                            read.alignment, partial_deletion_reads
                        )
                        for read_name in read_names:
                            if (
                                not read.is_del
                                and not read.is_refskip
                                and not read.alignment.is_secondary
                                and read.alignment.mapping_quality >= min_mapq
                                and read_name not in exclude_reads
                                and read.indel == 0
                                and read_name not in reads_with_flanking_indels
                            ):
                                read_seq = read.alignment.query_sequence
                                start_pos = read.query_position
                                end_pos = start_pos + 1
                                if end_pos < len(read_seq):
                                    hap = read_seq[start_pos:end_pos]
                                    if read_name not in read_haps:
                                        read_haps.setdefault(read_name, ["x"] * nvar)
                                    if hap.upper() == allele1.upper() or (
                                        snp_position in multi_allelic_sites
                                        and hap.upper()
                                        == multi_allelic_sites[snp_position]
                                    ):
                                        if allele_site not in absent_base_per_site or (
                                            allele_site in absent_base_per_site
                                            and absent_base_per_site[allele_site] != "1"
                                        ):
                                            read_haps[read_name][dsnp_index] = "1"
                                    elif hap.upper() == allele2.upper():
                                        if allele_site not in absent_base_per_site or (
                                            allele_site in absent_base_per_site
                                            and absent_base_per_site[allele_site] != "2"
                                        ):
                                            read_haps[read_name][dsnp_index] = "2"

        # for softclips starting at a predefined position, mark sites as 0 instead of x
        if check_clip and (
            self.clip_3p_positions != [] or self.clip_5p_positions != []
        ):
            for dsnp_index, allele_site in enumerate(het_sites):
                snp_position_gene1, allele1, allele2, *at = allele_site.split("_")
                snp_position = int(snp_position_gene1)
                for clip_position in sorted(self.clip_3p_positions):
                    if snp_position > clip_position:
                        for read in self._bamh.fetch(
                            self.nchr,
                            clip_position - clip_buffer,
                            clip_position + clip_buffer,
                        ):
                            read_name = self.get_read_name(read)
                            if read_name not in read_haps:
                                read_haps.setdefault(read_name, ["x"] * nvar)
                            if abs(read.reference_end - clip_position) < clip_buffer:
                                find_clip_3p = re.findall(
                                    self.clip_3p, read.cigarstring
                                )
                                if (
                                    find_clip_3p != []
                                    and int(find_clip_3p[0][:-1]) >= min_clip_len
                                ):
                                    read_haps[read_name][dsnp_index] = "0"
                for clip_position in sorted(self.clip_5p_positions, reverse=True):
                    if snp_position < clip_position:
                        for read in self._bamh.fetch(
                            self.nchr,
                            clip_position - clip_buffer,
                            clip_position + clip_buffer,
                        ):
                            read_name = self.get_read_name(read)
                            if read_name not in read_haps:
                                read_haps.setdefault(read_name, ["x"] * nvar)
                            if abs(read.reference_start - clip_position) < clip_buffer:
                                find_clip_5p = re.findall(
                                    self.clip_5p, read.cigarstring
                                )
                                if (
                                    find_clip_5p != []
                                    and int(find_clip_5p[0][:-1]) >= min_clip_len
                                ):
                                    read_haps[read_name][dsnp_index] = "0"
        return read_haps

    def remove_var(self, raw_read_haps, kept_sites, homo_sites):
        """remove variants that are not present after checking each read-haplotype"""
        bases_per_site = {}
        sites_to_remove = []
        absent_base_per_site = {}
        for i in range(len(self.het_sites)):
            for read, hap in raw_read_haps.items():
                base = hap[i]
                bases_per_site.setdefault(i, []).append(base)

        for pos in bases_per_site:
            bases = bases_per_site[pos]
            bases_x = bases.count("x")
            bases_0 = bases.count("0")
            bases_ref = bases.count("1")
            bases_alt = bases.count("2")
            this_var = self.het_sites[pos]
            if bases_x == len(bases) - bases_0:
                sites_to_remove.append(this_var)
            elif bases_ref + bases_alt == len(bases) - bases_x - bases_0 and (
                bases_alt <= 3 or bases_ref <= 3
            ):
                if bases_alt <= 3 and bases_ref <= 3:
                    sites_to_remove.append(this_var)
                else:
                    if this_var not in kept_sites:
                        sites_to_remove.append(this_var)
                    if this_var in homo_sites:
                        if bases_alt <= 3:
                            absent_base_per_site.setdefault(this_var, "2")
                        elif bases_ref <= 3:
                            absent_base_per_site.setdefault(this_var, "1")
        for var in sites_to_remove:
            if var in self.het_sites:
                self.het_sites.remove(var)
            if var in self.init_het_sites:
                self.init_het_sites.remove(var)
        return absent_base_per_site

    def allow_del_bases(self, pos):
        """
        During variant calling, allow "bases in deletions" for positions in
        two long deletions
        """
        if (
            self.del2_reads_partial != set()
            and self.del2_3p_pos1 <= pos <= self.del2_5p_pos2
        ):
            return True
        if (
            self.del1_reads_partial != set()
            and self.del1_3p_pos1 <= pos <= self.del1_5p_pos2
        ):
            return True
        return False

    def process_indel(self, pos, ref_seq, var_seq):
        """Translate pysam indel seq into real sequence"""
        if "+" in var_seq:
            ins_base = var_seq.split(re.findall(r"\+\d+", var_seq)[0])[1]
            indel_size = len(ins_base)
            var_seq = ref_seq + ins_base
        else:
            del_len = int(re.findall(r"\-\d+", var_seq)[0][1:])
            indel_size = del_len
            var_seq = ref_seq
            offset_pos = pos - self.offset
            ref_seq = self._refh.fetch(
                self.nchr_old,
                offset_pos - 1,
                offset_pos + del_len,
            ).upper()
        return ref_seq, var_seq, indel_size

    def get_candidate_pos(
        self,
        regions_to_check=[],
        min_read_support=5,
        min_vaf=0.11,
        white_list={},
        allow_del_region=[],
    ):
        """
        Get all polymorphic sites in the region, update self.candidate_pos
        """
        min_variant_frequency = min_vaf
        if self.min_vaf is not None:
            min_variant_frequency = self.min_vaf
        bamh = self._bamh
        pileups_raw = {}
        for pileupcolumn in bamh.pileup(
            self.nchr,
            self.left_boundary,
            self.right_boundary,
            truncate=True,
        ):
            pos = pileupcolumn.pos + 1
            pileups_raw.setdefault(
                pos,
                [a.upper() for a in pileupcolumn.get_query_sequences(add_indels=True)],
            )
        variants = {}
        variants_no_phasing = {}
        for pos in pileups_raw:
            all_bases = pileups_raw[pos]
            total_depth = len(all_bases)
            del_bases_count = all_bases.count("*")
            # get reference base
            offset_pos = pos - self.offset
            ref_seq_genome = self._refh.fetch(
                self.nchr_old, offset_pos - 1, offset_pos
            ).upper()

            if total_depth >= min_read_support and (
                del_bases_count < min_read_support
                or self.allow_del_bases(pos)
                or (
                    allow_del_region != []
                    and allow_del_region[0] < pos < allow_del_region[1]
                )
            ):
                all_bases = [a for a in all_bases if a != "*"]
                counter = Counter(all_bases)
                # include multi-allelic sites
                bases = counter.most_common(3)
                # homozygous
                if len(counter) == 1 or (
                    len(counter) >= 2
                    and bases[0][1]
                    > max(len(all_bases) - min_read_support, len(all_bases) * 0.85)
                ):
                    var_seq = bases[0][0]
                    ref_seq = ref_seq_genome
                    if var_seq != ref_seq:
                        # SNV and indels
                        if "-" not in var_seq and "+" not in var_seq:
                            # "homo" sites in large deletions should be put back into het sites
                            if (
                                self.allow_del_bases(pos)
                                and del_bases_count >= min_read_support
                                and pos not in self.homopolymer_sites
                            ):
                                variants.setdefault(pos, []).append((ref_seq, var_seq))
                            elif pos not in self.homopolymer_sites or (
                                pos in self.homopolymer_sites
                                and var_seq
                                not in self.homopolymer_sites[pos].split(",")
                            ):
                                self.homo_sites.append(f"{pos}_{ref_seq}_{var_seq}")
                        elif pos not in self.homopolymer_sites:
                            ref_seq, var_seq, indel_size = self.process_indel(
                                pos, ref_seq, var_seq
                            )
                            if indel_size < 25:
                                self.homo_sites.append(f"{pos}_{ref_seq}_{var_seq}")
                elif len(counter) >= 2:
                    found_ref = ref_seq_genome in [a[0] for a in bases]
                    if found_ref or pos in white_list:
                        for var_seq, var_count in bases:
                            ref_seq = ref_seq_genome
                            if var_seq != ref_seq and (
                                (
                                    var_count >= min_read_support
                                    and var_count >= total_depth * min_variant_frequency
                                )
                                or (
                                    self.targeted is False
                                    and var_count >= self.trusted_read_support
                                )
                            ):
                                # SNV
                                if "-" not in var_seq and "+" not in var_seq:
                                    if pos not in self.homopolymer_sites:
                                        variants.setdefault(pos, []).append(
                                            (ref_seq, var_seq)
                                        )
                                    else:
                                        prohibited_bases = self.homopolymer_sites[
                                            pos
                                        ].split(",")
                                        if var_seq not in prohibited_bases:
                                            if "1" in prohibited_bases:
                                                variants.setdefault(pos, []).append(
                                                    (ref_seq, var_seq)
                                                )
                                            else:
                                                variants_no_phasing.setdefault(
                                                    pos, (ref_seq, var_seq)
                                                )
                                # indels
                                elif pos not in self.homopolymer_sites:
                                    ref_seq, var_seq, indel_size = self.process_indel(
                                        pos, ref_seq, var_seq
                                    )
                                    if indel_size < 25:
                                        variants_no_phasing.setdefault(
                                            pos, (ref_seq, var_seq)
                                        )

        # exclude variants caused by shifted softclips of the big deletions
        excluded_variants = []
        for region in regions_to_check:
            var_to_check = [a for a in variants if region[0] < a < region[1]]
            excluded_variants += var_to_check
        for pos in variants:
            # for now, filter out multi-allelic sites, except in white list
            if pos not in excluded_variants:
                if len(variants[pos]) == 1:
                    ref_seq, var_seq = variants[pos][0]
                    self.candidate_pos.add(f"{pos}_{ref_seq}_{var_seq}")
                elif pos in white_list:
                    for ref_seq, var_seq in variants[pos]:
                        if var_seq != white_list[pos]:
                            self.candidate_pos.add(f"{pos}_{ref_seq}_{var_seq}")
                            break

        excluded_variants = []
        for region in regions_to_check:
            var_to_check = [a for a in variants_no_phasing if region[0] < a < region[1]]
            excluded_variants += var_to_check
        for pos in variants_no_phasing:
            if pos not in excluded_variants:
                ref_seq, var_seq = variants_no_phasing[pos]
                self.het_no_phasing.append(f"{pos}_{ref_seq}_{var_seq}")

    def remove_noisy_sites(self):
        """remove variants in predefined noisy sites"""
        problematic_sites = set()
        for site in self.het_sites:
            for region in self.noisy_region:
                if region[0] <= int(site.split("_")[0]) <= region[1]:
                    problematic_sites.add(site)
        for site in problematic_sites:
            if site in self.het_sites:
                self.het_sites.remove(site)

    @staticmethod
    def simplify_read_haps(read_haps):
        """Simplify read haplotypes for output"""
        haplotypes_to_reads = {}
        reads_to_haplotypes = {}
        for read in read_haps:
            hap = read_haps[read]
            haplotypes_to_reads.setdefault("".join(hap), []).append(read)
            reads_to_haplotypes.setdefault(read, "".join(hap))
        return haplotypes_to_reads, reads_to_haplotypes

    def check_variants_in_haplotypes(self, variant):
        """
        For variants not used in phasing, check which haplotypes they are in.
        """
        dreads = {}
        var_pos, ref, alt = variant.split("_")
        var_pos = int(var_pos)
        var_size = len(alt) - len(ref)
        if var_size < 0:
            indel_N = "N" * abs(var_size)
            indel_base_in_read = f"{ref[0]}{var_size}{indel_N}"
        elif var_size > 0:
            indel_base_in_read = f"{ref[0]}+{var_size}{alt[1:]}"
        ref_base = ref[0]
        alt_base = alt[0] if var_size == 0 else indel_base_in_read
        for pileupcolumn in self._bamh.pileup(
            self.nchr,
            var_pos - 1,
            var_pos,
            truncate=True,
        ):
            read_names = pileupcolumn.get_query_names()
            bases = [
                a.upper() for a in pileupcolumn.get_query_sequences(add_indels=True)
            ]
            for i, read in enumerate(read_names):
                if bases[i] == ref_base:
                    dreads.setdefault(read, "ref")
                elif bases[i] == alt_base:
                    dreads.setdefault(read, "alt")
                else:
                    dreads.setdefault(read, ".")
        return dreads

    @staticmethod
    def get_start_end(hap):
        """get range of positions that are not x"""
        haplen = len(hap)
        for nstart, base in enumerate(hap):
            if base != "x":
                break
        for nend in reversed(range(haplen)):
            if hap[nend] != "x":
                break
        return nstart, nend

    def get_hap_variant_ranges(self, hap):
        """get boundaries of (partial) haplotypes"""
        nstart, nend = self.get_start_end(hap)
        if nstart == 0:
            nstart_previous_pos = self.left_boundary
        else:
            nstart_previous = nstart - 1
            nstart_previous_pos = int(self.het_sites[nstart_previous].split("_")[0]) + 1
        if nend == len(hap) - 1:
            nend_next_pos = self.right_boundary
        else:
            nend_next = nend + 1
            nend_next_pos = int(self.het_sites[nend_next].split("_")[0]) - 1
        return nstart_previous_pos, nend_next_pos

    def output_variants_in_haplotypes(self, haps, reads, nonunique, known_del={}):
        """
        Summarize all variants in each haplotype.
        Output all variants and their genotypes.
        Haplotypes are different length, so a range (boundary) is reported
        """
        het_sites = self.het_sites
        haplotype_variants = {}
        haplotype_info = {}
        var_no_phasing = copy.deepcopy(self.het_no_phasing)
        for hap, hap_name in haps.items():
            haplotype_variants.setdefault(hap_name, [])
        # het sites not used in phasing
        if reads != {}:
            for var in var_no_phasing:
                genotypes = []
                var_reads = self.check_variants_in_haplotypes(var)
                haps_with_variant = []
                for hap, hap_name in haps.items():
                    hap_reads = reads[hap]
                    hap_reads_nonunique = [a for a in nonunique if hap in nonunique[a]]
                    genotype = self.get_genotype_in_hap(
                        var_reads, hap_reads, hap_reads_nonunique
                    )
                    genotypes.append(genotype)
                    if genotype == "1":
                        haps_with_variant.append(hap_name)
                if haps_with_variant == []:
                    self.het_no_phasing.remove(var)
                else:
                    for hap_name in haps_with_variant:
                        haplotype_variants[hap_name].append(var)
        # het sites and homo sites
        for hap, hap_name in haps.items():
            # find boundary for confident variant calling
            hap_bound_start, hap_bound_end = self.get_hap_variant_ranges(hap)
            is_truncated = []
            # het sites
            for i in range(len(hap)):
                if hap[i] == "2":
                    haplotype_variants[hap_name].append(het_sites[i])
                elif hap[i] in known_del:
                    del_name = known_del[hap[i]]
                    if del_name not in haplotype_variants[hap_name]:
                        haplotype_variants[hap_name].append(del_name)
            # homo sites
            filtered_homo_sites = self.homo_sites
            # check known deletions
            del1_name = None
            del2_name = None
            if known_del != {}:
                known_del_names = list(known_del.values())
                del1_name = known_del_names[0]
                if len(known_del_names) > 1:
                    del2_name = known_del_names[1]
            if del1_name in haplotype_variants[hap_name]:
                pos1 = self.del1_3p_pos1
                pos2 = self.del1_5p_pos2
                filtered_homo_sites = [
                    a
                    for a in filtered_homo_sites
                    if int(a.split("_")[0]) < pos1 or int(a.split("_")[0]) > pos2
                ]
            if del2_name in haplotype_variants[hap_name]:
                pos1 = self.del2_3p_pos1
                pos2 = self.del2_5p_pos2
                filtered_homo_sites = [
                    a
                    for a in filtered_homo_sites
                    if int(a.split("_")[0]) < pos1 or int(a.split("_")[0]) > pos2
                ]
            # haps with clips
            clip_position_5p = self.get_5pclip_from_hap(hap)
            if clip_position_5p is not None and clip_position_5p != 0:
                filtered_homo_sites = [
                    a
                    for a in filtered_homo_sites
                    if int(a.split("_")[0]) > clip_position_5p
                ]
                hap_bound_start = max(hap_bound_start, clip_position_5p)
                haplotype_variants[hap_name].append(f"{clip_position_5p}_clip_5p")
                if clip_position_5p > self.gene_start:
                    if "5p" not in is_truncated:
                        is_truncated.append("5p")
            clip_position_3p = self.get_3pclip_from_hap(hap)
            if clip_position_3p is not None and clip_position_3p != 0:
                filtered_homo_sites = [
                    a
                    for a in filtered_homo_sites
                    if int(a.split("_")[0]) < clip_position_3p
                ]
                hap_bound_end = min(hap_bound_end, clip_position_3p)
                haplotype_variants[hap_name].append(f"{clip_position_3p}_clip_3p")
                if clip_position_3p < self.gene_end:
                    if "3p" not in is_truncated:
                        is_truncated.append("3p")

            haplotype_variants[hap_name] += filtered_homo_sites

            hap_bound_start = max(hap_bound_start, self.left_boundary)
            hap_bound_end = min(hap_bound_end, self.right_boundary)
            boundary_gene2 = None
            if self.gene2_region is not None:
                bound1_in_other_gene = self.get_range_in_other_gene(hap_bound_start)
                bound2_in_other_gene = self.get_range_in_other_gene(hap_bound_end)
                if None not in [bound1_in_other_gene, bound2_in_other_gene]:
                    boundary_gene2 = [
                        min(bound1_in_other_gene, bound2_in_other_gene),
                        max(bound1_in_other_gene, bound2_in_other_gene),
                    ]
                else:
                    boundary_gene2 = [bound1_in_other_gene, bound2_in_other_gene]
            var_tmp = haplotype_variants[hap_name]
            var_tmp1 = [
                a
                for a in var_tmp
                if hap_bound_start <= int(a.split("_")[0]) <= hap_bound_end
            ]
            var_tmp1 = list(set(var_tmp1))
            var_tmp2 = sorted(var_tmp1, key=lambda x: int(x.split("_")[0]))
            haplotype_info.setdefault(
                hap_name,
                {
                    "variants": var_tmp2,
                    "boundary": [hap_bound_start, hap_bound_end],
                    "boundary_gene2": boundary_gene2,
                    "is_truncated": is_truncated,
                },
            )

        return haplotype_info

    def get_5pclip_from_hap(self, hap):
        """Given a haplotype, get its 5p clip position"""
        assert len(self.het_sites) == len(hap)
        clips_not_present = []
        for i, base in enumerate(hap):
            if i < len(hap) - 1:
                next_base = hap[i + 1]
                site_before = int(self.het_sites[i].split("_")[0])
                site_after = int(self.het_sites[i + 1].split("_")[0])
                for clip_position in self.clip_5p_positions:
                    if site_before < clip_position < site_after:
                        if base == "0" and next_base not in ["0", "x"]:
                            return clip_position
                        if next_base not in ["0", "x"]:
                            if base not in ["0", "x"] or (
                                base == "x" and site_after - clip_position < 5000
                            ):
                                clips_not_present.append(clip_position)
        if clips_not_present == self.clip_5p_positions:
            return 0
        return None

    def get_3pclip_from_hap(self, hap):
        """Given a haplotype, get its 3p clip position"""
        assert len(self.het_sites) == len(hap)
        clips_not_present = []
        for i, base in enumerate(hap):
            if i < len(hap) - 1:
                next_base = hap[i + 1]
                site_before = int(self.het_sites[i].split("_")[0])
                site_after = int(self.het_sites[i + 1].split("_")[0])
                for clip_position in self.clip_3p_positions:
                    if site_before < clip_position < site_after:
                        if next_base == "0" and base not in ["0", "x"]:
                            return clip_position
                        if base not in ["0", "x"]:
                            if next_base not in ["0", "x"] or (
                                next_base == "x" and clip_position - site_before < 5000
                            ):
                                clips_not_present.append(clip_position)
        if clips_not_present == self.clip_3p_positions:
            return 0
        return None

    def get_genotype_in_hap(self, var_reads, hap_reads, hap_reads_nonunique):
        """For a given variant, return its status in a haplotype"""
        hap_reads_contain_var = [var_reads[a] for a in hap_reads if a in var_reads]
        if len(hap_reads_contain_var) < 3:
            hap_reads_contain_var += [
                var_reads[a] for a in hap_reads_nonunique if a in var_reads
            ]
        if len(hap_reads_contain_var) >= 3:
            hap_reads_contain_var_counter = Counter(hap_reads_contain_var).most_common(
                2
            )
            if len(hap_reads_contain_var_counter) == 1 or hap_reads_contain_var_counter[
                1
            ][1] <= min(hap_reads_contain_var_counter[0][1] * 0.15, 2):
                if hap_reads_contain_var_counter[0][0] == "alt":
                    return "1"
                elif hap_reads_contain_var_counter[0][0] == "ref":
                    return "0"
        return "."

    @staticmethod
    def update_reads_for_deletions(
        raw_read_haps, het_sites, n1, n2, del_reads_partial, base, del_name
    ):
        """
        For reads carrying known big deletions, update read haplotype to
        reflect the deletion. This is needed for downstream phasing
        """
        pos1 = -1
        pos2 = -1
        for i, var in enumerate(het_sites):
            if int(var.split("_")[0]) > n1:
                pos1 = i
                break
        for i, var in enumerate(het_sites):
            if int(var.split("_")[0]) > n2:
                pos2 = i
                break
        if pos1 != -1 and pos2 != -1:
            if pos1 < pos2:
                for read in del_reads_partial:
                    if read in raw_read_haps:
                        hap = list(raw_read_haps[read])
                        for i in range(pos1, pos2):
                            hap[i] = base
                        raw_read_haps[read] = "".join(hap)
            elif pos1 == pos2:
                het_sites.insert(pos1, del_name)
                for read in raw_read_haps:
                    hap = list(raw_read_haps[read])
                    hap.insert(pos1, "x")
                    if read in del_reads_partial:
                        hap[pos1] = base
                    elif (
                        hap[pos1 - 1] == "0"
                        and pos1 - 1 >= 0
                        and hap[pos1 + 1] == "0"
                        and pos1 + 1 < len(hap)
                    ):
                        hap[pos1] = "0"
                    else:
                        flanking_left = hap[max(0, pos1 - 2) : pos1]
                        flanking_right = hap[
                            min(pos1 + 1, len(hap)) : min(pos1 + 3, len(hap))
                        ]
                        if "x" not in flanking_left and "x" not in flanking_right:
                            hap[pos1] = "1"
                    raw_read_haps[read] = "".join(hap)
        return raw_read_haps, het_sites

    def get_read_counts(self, uniquely_supporting_haps):
        """
        Get unique supporting read counts for each haplotype
        over a region where all haplotypes have defined bases
        """
        if uniquely_supporting_haps == {}:
            return {}
        nhap = len(uniquely_supporting_haps)
        nvar = len(self.het_sites)
        # depth per site for haplotype
        hap_bases = []
        for i in range(nhap):
            hap_bases.append([])
        j = 0
        for hap in uniquely_supporting_haps:
            reads = uniquely_supporting_haps[hap]
            for i in range(nvar):
                bases = [a[i] for a in reads]
                hap_bases[j].append(len(bases) - bases.count("x") - bases.count("0"))
            j += 1
        ranges = []
        for i in range(nvar):
            if min([hap_bases[a][i] for a in range(nhap)]) >= 5:
                for j in range(nvar):
                    if j > i:
                        if min([hap_bases[a][j] for a in range(nhap)]) < 5:
                            ranges.append([i, j])
                            break
                    if j == nvar - 1 and j > i:
                        ranges.append([i, nvar - 1])
        if ranges == []:
            return None
        longest_range = sorted(ranges, key=lambda x: x[1] - x[0], reverse=True)[0]
        mid = int((longest_range[1] + longest_range[0]) / 2)
        nstart = max(mid - 1, longest_range[0])
        nend = min(mid + 1, longest_range[1])
        if nend == nstart:
            return None

        read_count = {}
        for hap in uniquely_supporting_haps:
            reads = uniquely_supporting_haps[hap]
            lreads = [a for a in reads if a[nstart:nend] != "x" * (nend - nstart)]
            read_count.setdefault(hap, len(lreads))
        return read_count

    def phase_haps(self, raw_read_haps, debug=False):
        """
        Assemble and evaluate haplotypes
        """
        total_depth = self.region_avg_depth.median
        min_support = 4
        if (
            self.targeted
            and self.min_haplotype_frequency is not None
            and total_depth is not None
        ):
            min_support = max(min_support, total_depth * self.min_haplotype_frequency)
        het_sites = self.het_sites
        haplotypes_to_reads, raw_read_haps = self.simplify_read_haps(raw_read_haps)

        ass_haps = []
        original_haps = []
        nvar = len(het_sites)
        hcn = 0
        if nvar == 0:
            return ([], [], 0, {}, {}, {}, None)
        elif nvar == 1:
            ass_haps = ["1", "2"]
            original_haps = ["1", "2"]
        else:
            pivot_index, _ = self.get_pivot_site_index()
            hap_graph = VariantGraph(
                raw_read_haps, pivot_index, figure_id=self.sample_id
            )
            ass_haps, original_haps, hcn = hap_graph.run(
                debug=debug, make_plot=debug, min_count=min_support
            )

        if ass_haps == []:
            return (ass_haps, original_haps, hcn, {}, {}, raw_read_haps, None)

        (
            uniquely_supporting_reads,
            nonuniquely_supporting_reads,
            read_counts,
        ) = self.get_read_support(raw_read_haps, haplotypes_to_reads, ass_haps)

        # remove spurious ones
        ass_haps = self.adjust_spurious_haplotypes(uniquely_supporting_reads)
        (
            uniquely_supporting_reads,
            nonuniquely_supporting_reads,
            read_counts,
        ) = self.get_read_support(raw_read_haps, haplotypes_to_reads, ass_haps)

        # remove low-support ones
        read_counts = [
            len(uniquely_supporting_reads[a]) for a in uniquely_supporting_reads
        ]
        read_counts = sorted(read_counts)
        # more stringent if one copy is low but others are high
        if (
            len(read_counts) > 2
            and read_counts[0] <= 4
            and read_counts[1] >= 12
            and "x" not in "".join(ass_haps)
            and min_support == 4
        ):
            min_support = 5

        ass_haps = [
            a
            for a in uniquely_supporting_reads
            if len(uniquely_supporting_reads[a]) >= min_support
        ]
        (
            uniquely_supporting_reads,
            nonuniquely_supporting_reads,
            read_counts,
        ) = self.get_read_support(raw_read_haps, haplotypes_to_reads, ass_haps)

        return (
            ass_haps,
            original_haps,
            hcn,
            uniquely_supporting_reads,
            nonuniquely_supporting_reads,
            raw_read_haps,
            read_counts,
        )

    def phase_haps_catch_error(self, raw_read_haps, debug=False):
        """Run phase_haps, catch error"""
        try:
            (
                ass_haps,
                original_haps,
                hcn,
                uniquely_supporting_reads,
                nonuniquely_supporting_reads,
                raw_read_haps,
                read_counts,
            ) = self.phase_haps(raw_read_haps, debug=debug)
        except Exception:
            logging.warning(
                "Did not phase haplotypes successfully, possibly due to low coverage. See error message below"
            )
            traceback.print_exc()
            return (
                self.GeneCall(
                    sites_for_phasing=self.het_sites,
                    het_sites_not_used_in_phasing=self.het_no_phasing,
                    homozygous_sites=self.homo_sites,
                    read_details={a: "".join(b) for a, b in raw_read_haps.items()},
                    genome_depth=self.mdepth,
                    region_depth=self.region_avg_depth._asdict(),
                    sample_sex=self.sample_sex,
                    heterozygous_sites=self.init_het_sites,
                ),
                None,
            )
        return (
            None,
            (
                ass_haps,
                original_haps,
                hcn,
                uniquely_supporting_reads,
                nonuniquely_supporting_reads,
                raw_read_haps,
                read_counts,
            ),
        )

    def get_read_support(self, raw_read_haps, haplotypes_to_reads, ass_haps):
        """Find uniquely and nonuniquely supporting reads for given haplotypes"""
        read_support = VariantGraph.match_reads_and_haplotypes(raw_read_haps, ass_haps)
        uniquely_supporting_haps = read_support.unique
        read_counts = self.get_read_counts(uniquely_supporting_haps)

        uniquely_supporting_reads = {}
        for hap in ass_haps:
            uniquely_supporting_reads.setdefault(hap, [])
        for hap in uniquely_supporting_haps:
            for read_hap in uniquely_supporting_haps[hap]:
                uniquely_supporting_reads[hap] += haplotypes_to_reads[read_hap]
        for hap in uniquely_supporting_haps:
            uniquely_supporting_reads[hap] = sorted(
                list(set(uniquely_supporting_reads[hap]))
            )

        nonuniquely_supporting_reads = {}
        for read in read_support.by_read:
            num_matches = len(read_support.by_read[read])
            if num_matches > 1:
                nonuniquely_supporting_reads.setdefault(
                    read, read_support.by_read[read]
                )
        return (
            uniquely_supporting_reads,
            nonuniquely_supporting_reads,
            read_counts,
        )

    def compare_depth(self, haplotypes, ass_haps, loose=False, stringent=False):
        """
        For each haplotype, identify the variants where it's different
        from other haplotypes. Check depth at those variant sites and
        see if the depth suggests twice coverage.
        """
        if haplotypes is None or len(haplotypes) == 1:
            return []
        two_cp_haps = []
        bamh = self._bamh
        hap_name_to_seq = {v: k for k, v in ass_haps.items()}
        bound_start = []
        bound_end = []
        for hap, hap_info in haplotypes.items():
            bound = hap_info["boundary"]
            hap_clip_5p = self.get_5pclip_from_hap(hap_name_to_seq[hap])
            hap_clip_3p = self.get_3pclip_from_hap(hap_name_to_seq[hap])
            if hap_clip_5p is not None and hap_clip_5p != 0:
                bound_start.append(hap_clip_5p)
            else:
                bound_start.append(bound[0])
            if hap_clip_3p is not None and hap_clip_3p != 0:
                bound_end.append(hap_clip_3p)
            else:
                bound_end.append(bound[1])
        nstart = max(bound_start)
        nend = min(bound_end)
        variants = set()
        for hap in haplotypes:
            vars = haplotypes[hap]["variants"]
            for var in vars:
                if len(var.split("_")) == 3:
                    pos, ref, alt = var.split("_")
                    pos = int(pos)
                    if nstart < pos < nend and var in self.het_sites:
                        variants.add(var)

        for hap in haplotypes:
            sites = {}
            other_haps = [a for a in haplotypes.keys() if a != hap]
            other_cn = len(other_haps)
            this_hap_var = haplotypes[hap]["variants"]
            other_haps_var = []
            for a in other_haps:
                other_haps_var += haplotypes[a]["variants"]
            for var in variants:
                pos, ref, alt = var.split("_")
                pos = int(pos)
                if var in this_hap_var and var not in other_haps_var:
                    sites.setdefault(pos, alt)
                elif var not in this_hap_var and other_haps_var.count(var) == other_cn:
                    sites.setdefault(pos, ref)

            counts = []
            for pos in sites:
                hap_base = sites[pos]
                for pileupcolumn in bamh.pileup(
                    self.nchr,
                    pos - 1,
                    pos,
                    truncate=True,
                    min_base_quality=self.MEAN_BASE_QUAL,
                ):
                    bases = [a.upper() for a in pileupcolumn.get_query_sequences()]
                    base_num = bases.count(hap_base)
                    counts.append([base_num, len(bases) - base_num])

            probs = []
            nsites = len(sites)
            for n1, n2 in counts:
                this_site_prob = self.depth_prob(n1, n2 / other_cn)
                probs.append(this_site_prob)
            probs_fil = [a for a in probs if a is not None and a[0] < 0.25]
            if stringent is True:
                if len(probs_fil) >= nsites * 0.8 and nsites >= 5:
                    two_cp_haps.append(hap)
            elif len(probs_fil) >= nsites * 0.6 and nsites >= 5:
                two_cp_haps.append(hap)
            elif loose is True:
                if len(probs_fil) >= nsites * 0.5 and nsites >= 5:
                    two_cp_haps.append(hap)

        # there can only be one such haplotype
        if len(two_cp_haps) > 1:
            return []
        return two_cp_haps

    def adjust_spurious_haplotypes(self, uniquely_supporting_reads, flanking_bp=10):
        """Identify spurious haplotypes caused by locally misaligned reads"""
        passing_haplotypes = list(uniquely_supporting_reads.keys())
        suspicious_hap_pair = []
        lhap = uniquely_supporting_reads.keys()
        for hap1 in lhap:
            for hap2 in lhap:
                if hap1 != hap2:
                    nmatch = 0
                    nmismatch = 0
                    mismatch_sites = []
                    for i, base1 in enumerate(hap1):
                        base2 = hap2[i]
                        if "x" not in [base1, base2]:
                            if base1 == base2:
                                nmatch += 1
                            elif base1 in ["1", "2"] and base2 in ["1", "2"]:
                                nmismatch += 1
                                mismatch_sites.append(self.het_sites[i])
                    if nmatch >= 5 and nmismatch == 1 and len(mismatch_sites) == 1:
                        mismatch_pos = int(mismatch_sites[0].split("_")[0])
                        hap1_reads = uniquely_supporting_reads[hap1]
                        hap2_reads = uniquely_supporting_reads[hap2]
                        hap_pair = None
                        if len(hap1_reads) <= 5 and len(hap2_reads) >= 6:
                            hap_pair = [hap2, hap1]
                        elif len(hap2_reads) <= 5 and len(hap1_reads) >= 6:
                            hap_pair = [hap1, hap2]
                        if (
                            hap_pair is not None
                            and [
                                hap_pair,
                                mismatch_pos,
                            ]
                            not in suspicious_hap_pair
                        ):
                            suspicious_hap_pair.append([hap_pair, mismatch_pos])

        for hap_pair, mismatch_pos in suspicious_hap_pair:
            hap1, hap2 = hap_pair
            hap1_reads = uniquely_supporting_reads[hap1]
            hap2_reads = uniquely_supporting_reads[hap2]
            hap1_reads_at_pos = []
            hap2_reads_at_pos = []
            for pileupcolumn in self._bamh.pileup(
                self.nchr,
                mismatch_pos - 1,
                mismatch_pos,
                truncate=True,
                min_base_quality=self.MEAN_BASE_QUAL,
            ):
                for read in pileupcolumn.pileups:
                    if not read.is_del and not read.is_refskip:
                        read_name = self.get_read_name(read.alignment)
                        read_seq = read.alignment.query_sequence
                        read_pos = read.query_position
                        if (
                            read_name in hap1_reads + hap2_reads
                            and read_pos >= flanking_bp
                            and read_pos + flanking_bp < len(read_seq)
                        ):
                            start_pos = read_pos - flanking_bp
                            end_pos = read_pos + flanking_bp
                            if read_name in hap1_reads:
                                hap1_reads_at_pos.append(read_seq[start_pos:end_pos])
                            if read_name in hap2_reads:
                                hap2_reads_at_pos.append(read_seq[start_pos:end_pos])

            if set(hap1_reads_at_pos).intersection(set(hap2_reads_at_pos)) != set():
                if hap2 in passing_haplotypes:
                    passing_haplotypes.remove(hap2)
        return passing_haplotypes

    @staticmethod
    def check_linking_read(aln1, aln2, reverse=False):
        """Determine the direction of links between two alignments"""
        if reverse is False:
            if aln1[-1] != "x" and aln2[0] != "x":
                return "1-2"
            if aln2[-1] != "x" and aln1[0] != "x":
                return "2-1"
        else:
            if aln1[-1] != "x" and aln2[-1] != "x":
                return "0-0"
            if aln1[0] != "x" and aln2[0] != "x":
                return "0-0"
        return None

    def phase_alleles(
        self,
        uniq_reads,
        nonuniquely_supporting_reads,
        raw_read_haps,
        ass_haps,
        reverse=False,
        min_read=2,
        haps_to_exclude=[],
    ):
        """
        Phase haplotypes into alleles using read evidence
        """
        new_reads = {}
        # unique
        for hap in uniq_reads:
            if hap not in haps_to_exclude:
                for read in uniq_reads[hap]:
                    short_name = read.split("_sup")[0]
                    new_reads.setdefault(short_name, []).append({read: [hap]})
        # nonunique
        for read, supported_haps in nonuniquely_supporting_reads.items():
            new_supported_haps = [a for a in supported_haps if a not in haps_to_exclude]
            short_name = read.split("_sup")[0]
            new_reads.setdefault(short_name, []).append({read: new_supported_haps})

        (
            nondirected_links,
            directed_links,
            directed_links_loose,
        ) = self.get_directed_links(new_reads, raw_read_haps, ass_haps, reverse)

        read_links = {}
        for hap_link in nondirected_links:
            if len(nondirected_links[hap_link]) >= min_read:
                hap1, hap2 = hap_link.split("-")
                read_links.setdefault(hap1, []).append(hap2)
                read_links.setdefault(hap2, []).append(hap1)
        read_links = dict(
            sorted(read_links.items(), key=lambda item: len(item[1]), reverse=True)
        )
        alleles = Phaser.get_alleles_from_links(read_links, ass_haps.values())
        linked_haps = [a for a in alleles]
        alleles = self.filter_alleles(alleles, ass_haps, haps_to_exclude)

        return (
            alleles,
            read_links,
            {a: len(b) for a, b in directed_links.items()},
            {a: Counter(b) for a, b in directed_links_loose.items()},
            linked_haps,
        )

    @staticmethod
    def filter_alleles(alleles, ass_haps, haps_to_exclude):
        """Filter out possibly incorrectly phased alleles"""
        # scenarios to filter out
        # 1. more than two alleles
        if len(alleles) > 2:
            alleles = []
        # 2. two alleles sharing haplotypes
        if len(alleles) == 2:
            if set(alleles[0]).intersection(set(alleles[1])) != set():
                alleles = []
        # 3. not all haplotypes are included in alleles
        haps_in_alleles = []
        for allele in alleles:
            for hap in allele:
                haps_in_alleles.append(hap)
        haps_considered = [
            hap_name
            for (hap_seq, hap_name) in ass_haps.items()
            if hap_seq not in haps_to_exclude
        ]
        if len(alleles) == 2 and len(haps_in_alleles) < len(haps_considered):
            alleles = []
        if len(alleles) == 1 and len(haps_in_alleles) < len(haps_considered):
            if len(haps_in_alleles) == len(haps_considered) - 1 and len(
                haps_in_alleles
            ) in [2, 3]:
                # 1/2 or 1/3 are okay, add the single haplotype as the second allele
                alleles.append([a for a in haps_considered if a not in haps_in_alleles])
            else:
                alleles = []
        return alleles

    def get_directed_links(self, new_reads, raw_read_haps, ass_haps, reverse):
        """Get links between haplotypes from reads"""
        nondirected_links = {}
        directed_links = {}
        directed_links_loose = {}
        for read, hap_info in new_reads.items():
            nsegment = len(hap_info)
            if nsegment >= 2:
                links_found = set()
                for i in range(nsegment):
                    for j in range(i + 1, nsegment):
                        hap_info1 = hap_info[i]
                        hap_info2 = hap_info[j]
                        read1, haps1 = list(hap_info1.items())[0]
                        read2, haps2 = list(hap_info2.items())[0]
                        aln1 = raw_read_haps[read1]
                        aln2 = raw_read_haps[read2]
                        check_link = self.check_linking_read(aln1, aln2, reverse)
                        if len(haps1) == 1 and len(haps2) == 1:
                            hap1 = haps1[0]
                            hap2 = haps2[0]
                            if hap1 != hap2:
                                hap1_renamed = ass_haps[hap1]
                                hap2_renamed = ass_haps[hap2]
                                link_to_add = f"{hap1_renamed}-{hap2_renamed}"
                                if link_to_add not in links_found:
                                    links_found.add(link_to_add)
                                    nondirected_links.setdefault(
                                        link_to_add, []
                                    ).append(1)
                                if check_link is not None and check_link in [
                                    "1-2",
                                    "0-0",
                                ]:
                                    directed_links.setdefault(
                                        f"{hap1_renamed}-{hap2_renamed}", []
                                    ).append(1)
                                    directed_links_loose.setdefault(
                                        hap1_renamed, []
                                    ).append(hap2_renamed)
                                elif check_link is not None and check_link in [
                                    "2-1",
                                    "0-0",
                                ]:
                                    directed_links.setdefault(
                                        f"{hap2_renamed}-{hap1_renamed}", []
                                    ).append(1)
                                    directed_links_loose.setdefault(
                                        hap2_renamed, []
                                    ).append(hap1_renamed)
                        elif len(haps1) == 1 or len(haps2) == 1:
                            for hap1 in haps1:
                                for hap2 in haps2:
                                    if hap1 != hap2:
                                        hap1_renamed = ass_haps[hap1]
                                        hap2_renamed = ass_haps[hap2]
                                        if check_link is not None and check_link in [
                                            "1-2",
                                            "0-0",
                                        ]:
                                            directed_links_loose.setdefault(
                                                hap1_renamed, []
                                            ).append(hap2_renamed)
                                        elif check_link is not None and check_link in [
                                            "2-1",
                                            "0-0",
                                        ]:
                                            directed_links_loose.setdefault(
                                                hap2_renamed, []
                                            ).append(hap1_renamed)
        return nondirected_links, directed_links, directed_links_loose

    @staticmethod
    def get_alleles_from_links(read_links, total_haps):
        """phase alleles from read_links between reads"""
        alleles = []
        if read_links != {}:
            alleles = [[list(read_links.keys())[0]] + list(read_links.values())[0]]
            for hap1 in read_links:
                for hap2 in read_links[hap1]:
                    hap1_in = sum([hap1 in a for a in alleles])
                    hap2_in = sum([hap2 in a for a in alleles])
                    if hap1_in == 0 and hap2_in == 0:
                        alleles.append([hap1, hap2])
                    elif hap1_in == 0:
                        for a in alleles:
                            if hap2 in a:
                                a_index = alleles.index(a)
                                if hap1 not in alleles[a_index]:
                                    alleles[a_index].append(hap1)
                                if hap2 not in alleles[a_index]:
                                    alleles[a_index].append(hap2)
                    else:
                        for a in alleles:
                            if hap1 in a:
                                a_index = alleles.index(a)
                                if hap2 not in alleles[a_index]:
                                    alleles[a_index].append(hap2)
                                if hap1 not in alleles[a_index]:
                                    alleles[a_index].append(hap1)
        # merge alleles
        while True:
            to_merge = []
            for hap in total_haps:
                hap_found_in_alleles = [hap in a for a in alleles]
                if hap_found_in_alleles.count(True) > 1:
                    to_merge.append(hap)
                    break
            if to_merge == []:
                break
            new_alleles = []
            hap = to_merge[0]
            merged = []
            for each_allele in alleles:
                if hap not in each_allele:
                    new_alleles.append(each_allele)
                else:
                    merged += each_allele
            merged = sorted(list(set(merged)))
            new_alleles.append(merged)
            alleles = new_alleles
        return alleles

    def add_homo_sites(self, min_no_var_region_size=10000, max_homo_var_to_add=10):
        """add some homozygous sites to het sites in long homozygous regions"""
        # if self.het_sites == []:
        #    return []
        self.het_sites = sorted(self.het_sites, key=lambda x: int(x.split("_")[0]))
        homo_sites_to_add = []
        het_pos = []
        for var in self.het_sites:
            pos = int(var.split("_")[0])
            if self.left_boundary < pos < self.right_boundary:
                het_pos.append(pos)
        if het_pos == []:
            for var in self.homo_sites:
                pos, ref, alt = var.split("_")
                if len(ref) == 1 and len(alt) == 1:
                    homo_sites_to_add.append(var)
            if self.homo_sites == []:
                full_range = self.right_boundary - self.left_boundary
                interval_size = int(full_range / 4)
                positions = [
                    self.left_boundary + interval_size,
                    self.left_boundary + interval_size * 2,
                    self.left_boundary + interval_size * 3,
                ]
                for var_pos in positions:
                    ref_base = self._refh.fetch(
                        self.nchr_old, var_pos - self.offset - 1, var_pos - self.offset
                    ).upper()
                    var_base = [a for a in ["A", "C", "G", "T"] if a != ref_base][0]
                    new_var = f"{var_pos}_{ref_base}_{var_base}"
                    homo_sites_to_add.append(new_var)
        else:
            min_pos = min(het_pos)
            max_pos = max(het_pos)
            if min_pos - self.left_boundary > min_no_var_region_size:
                for var in self.homo_sites:
                    pos, ref, alt = var.split("_")
                    if int(pos) < min_pos and len(ref) == 1 and len(alt) == 1:
                        homo_sites_to_add.append(var)
            if self.right_boundary - max_pos > min_no_var_region_size:
                for var in self.homo_sites:
                    pos, ref, alt = var.split("_")
                    if int(pos) > max_pos and len(ref) == 1 and len(alt) == 1:
                        homo_sites_to_add.append(var)
            for i in range(len(self.het_sites) - 1):
                interval_start = int(self.het_sites[i].split("_")[0])
                interval_end = int(self.het_sites[i + 1].split("_")[0])
                if interval_end - interval_start > min_no_var_region_size:
                    for var in self.homo_sites:
                        pos, ref, alt = var.split("_")
                        if (
                            interval_start < int(pos) < interval_end
                            and len(ref) == 1
                            and len(alt) == 1
                        ):
                            homo_sites_to_add.append(var)
        if homo_sites_to_add != []:
            homo_sites_to_add = sorted(
                homo_sites_to_add, key=lambda x: int(x.split("_")[0])
            )
            homo_sites_to_add_size = len(homo_sites_to_add)
            homo_sites_to_add = homo_sites_to_add[
                :: int(np.ceil(homo_sites_to_add_size / max_homo_var_to_add))
            ]
            self.het_sites += homo_sites_to_add
            self.het_sites = sorted(self.het_sites)
            self.remove_noisy_sites()
        return homo_sites_to_add

    def find_fusion(self, ass_haps):
        """Call fusion based on haplotypes"""
        # update two-copy haplotypes
        two_cp_haps = self.update_twp_cp_in_fusion_cases(ass_haps)
        fusions_called = {}
        for hap, hap_name in ass_haps.items():
            if hap.endswith("x") is False and hap.startswith("x") is False:
                if (hap.endswith("0") is False and hap.startswith("0") is True) or (
                    hap.endswith("0") is True and hap.startswith("0") is False
                ):
                    new_hap, all_sites = self.new_hap_for_breakpoint(hap)
                    fusion_breakpoint_index = self.get_fusion_breakpoint_index(
                        hap, new_hap
                    )
                    if fusion_breakpoint_index is not None:
                        bp1 = int(all_sites[fusion_breakpoint_index].split("_")[0])
                        bp2 = self.get_range_in_other_gene(bp1, search_range=1000)
                        bp3 = int(all_sites[fusion_breakpoint_index - 1].split("_")[0])
                        bp4 = self.get_range_in_other_gene(bp3, search_range=1000)
                        if bp2 is not None:
                            if bp1 < bp2:
                                fusion_breakpoint = (
                                    (bp3, bp1),
                                    (bp4, bp2),
                                )
                            else:
                                fusion_breakpoint = (
                                    (bp4, bp2),
                                    (bp3, bp1),
                                )
                            fusions_called.setdefault(hap_name, {})
                            fusion_type = self.get_fusion_type(hap)
                            fusions_called[hap_name].setdefault("type", fusion_type)
                            fusions_called[hap_name].setdefault("sequence", new_hap)
                            fusions_called[hap_name].setdefault(
                                "breakpoint", fusion_breakpoint
                            )
        return two_cp_haps, fusions_called

    def get_fusion_type(self, hap):
        """Fusion type: deletion or duplication"""
        fusion_type = None
        if self.call_fusion == "5p":
            if hap.endswith("0") is False and hap.startswith("0") is True:
                fusion_type = "duplication"
            elif hap.endswith("0") is True and hap.startswith("0") is False:
                fusion_type = "deletion"
        elif self.call_fusion == "3p":
            if hap.endswith("0") is False and hap.startswith("0") is True:
                fusion_type = "deletion"
            elif hap.endswith("0") is True and hap.startswith("0") is False:
                fusion_type = "duplication"
        return fusion_type

    @staticmethod
    def update_twp_cp_in_fusion_cases(ass_haps):
        """Update two-copy haplotypes based on the presence of gene/paralogs"""
        two_cp_haps = []
        if True not in [a.startswith("x") or a.endswith("x") for a in ass_haps]:
            gene1s = [
                a
                for a in ass_haps
                if a.endswith("0") is False and a.startswith("0") is False
            ]
            gene2s = [
                a
                for a in ass_haps
                if a.endswith("0") is True and a.startswith("0") is True
            ]
            fusions = [
                a
                for a in ass_haps
                if (a.endswith("0") is False and a.startswith("0") is True)
                or (a.endswith("0") is True and a.startswith("0") is False)
            ]
            if fusions == [] and len(ass_haps) < 4:
                if len(gene1s) == 1 and ass_haps[gene1s[0]] not in two_cp_haps:
                    two_cp_haps.append(ass_haps[gene1s[0]])
                if len(gene2s) == 1 and ass_haps[gene2s[0]] not in two_cp_haps:
                    two_cp_haps.append(ass_haps[gene2s[0]])
            # homozygous fusion
            elif len(fusions) == 1 and len(ass_haps) == 1:
                two_cp_haps.append(ass_haps[fusions[0]])
        return two_cp_haps

    def new_hap_for_breakpoint(self, hap):
        """
        Get the haplotype sequence for breakpoint identification
        This is ideally based on PSVs defined in self.fusion_gene_def_variants
        """
        new_hap = ""
        if self.fusion_gene_def_variants != []:
            all_sites = self.fusion_gene_def_variants
            for var_site in all_sites:
                base = "1"
                if var_site in self.homo_sites:
                    base = "2"
                elif var_site in self.het_sites:
                    base = hap[self.het_sites.index(var_site)]
                new_hap += base
        else:
            all_sites = sorted(
                self.homo_sites + self.het_sites, key=lambda x: int(x.split("_")[0])
            )
            if self.clip_5p_positions != []:
                all_sites = [
                    a
                    for a in all_sites
                    if int(a.split("_")[0]) > max(self.clip_5p_positions)
                ]
            if self.clip_3p_positions != []:
                all_sites = [
                    a
                    for a in all_sites
                    if int(a.split("_")[0]) < min(self.clip_3p_positions)
                ]
            for var_site in all_sites:
                if var_site in self.homo_sites:
                    new_hap += "2"
                elif var_site in self.het_sites:
                    new_hap += hap[self.het_sites.index(var_site)]
        return new_hap, all_sites

    @staticmethod
    def get_fusion_breakpoint_index(hap, new_hap):
        """Infer the switch from gene1 sequence to gene2 sequence"""
        # 2s to 1s
        if hap.startswith("0") is True and hap.endswith("0") is False:
            counts = []
            for i, _ in enumerate(new_hap):
                counts.append(
                    new_hap[:i].count("2") + new_hap[i:].count("1"),
                )
            bp_index = counts.index(max(counts))
            if bp_index == 0 or bp_index == len(counts) - 1:
                return None
            return bp_index
        # 1s to 2s
        if hap.startswith("0") is False and hap.endswith("0") is True:
            counts = []
            for i, _ in enumerate(new_hap):
                counts.append(
                    new_hap[:i].count("1") + new_hap[i:].count("2"),
                )
            bp_index = counts.index(max(counts))
            if bp_index == 0 or bp_index == len(counts) - 1:
                return None
            return bp_index
        return None

    def get_cn2_haplotype(
        self, read_counts, ass_haps, min_cn1_read_count=10, prob_cutoff=0.05
    ):
        """
        Check if the haplotype with the highest depth has twice the reads
        of the haplotype with the second highest depth
        """
        if read_counts is None:
            return []
        if len(read_counts) < 2:
            return []
        two_cp_haps = []
        haps = list(read_counts.keys())
        counts = list(read_counts.values())
        max_count = max(counts)
        cp2_hap = haps[counts.index(max_count)]
        others_max = sorted(counts, reverse=True)[1]
        probs = self.depth_prob(max_count, others_max)
        if (
            probs is not None
            and probs[0] < prob_cutoff
            and others_max >= min_cn1_read_count
        ):
            two_cp_haps.append(ass_haps[cp2_hap])
        return two_cp_haps

    def call(self):
        """Main function to phase haplotypes and call copy numbers"""
        if self.check_coverage_before_analysis() is False:
            return self.GeneCall(
                genome_depth=self.mdepth,
                region_depth=self.region_avg_depth._asdict(),
                sample_sex=self.sample_sex,
                phase_region=f"{self.genome_build}:{self.nchr}:{self.left_boundary}-{self.right_boundary}",
            )
        self.get_homopolymer()
        self.find_big_deletion()

        if self.deletion1_size is not None:
            self.del1_reads, self.del1_reads_partial = self.get_long_del_reads(
                self.del1_3p_pos1,
                self.del1_3p_pos2,
                self.del1_5p_pos1,
                self.del1_5p_pos2,
                self.deletion1_size,
            )
        if self.deletion2_size is not None:
            self.del2_reads, self.del2_reads_partial = self.get_long_del_reads(
                self.del2_3p_pos1,
                self.del2_3p_pos2,
                self.del2_5p_pos1,
                self.del2_5p_pos2,
                self.deletion2_size,
            )

        # scan for polymorphic sites
        regions_to_check = []
        if self.del2_reads_partial != set():
            regions_to_check += [
                [self.del2_3p_pos1, self.del2_3p_pos2],
                [self.del2_5p_pos1, self.del2_5p_pos2],
            ]
        if self.del1_reads_partial != set():
            regions_to_check += [
                [self.del1_3p_pos1, self.del1_3p_pos2],
                [self.del1_5p_pos1, self.del1_5p_pos2],
            ]
        self.get_candidate_pos(regions_to_check=regions_to_check)

        self.het_sites = sorted(list(self.candidate_pos))
        self.remove_noisy_sites()
        self.init_het_sites = [a for a in self.het_sites]
        homo_sites_to_add = self.add_homo_sites()

        # add pivot site
        if self.pivot_site is not None:
            if self.pivot_site not in [
                int(a.split("_")[0]) for a in self.het_sites + self.add_sites
            ]:
                ref_base = self._refh.fetch(
                    self.nchr_old,
                    self.pivot_site - self.offset - 1,
                    self.pivot_site - self.offset,
                ).upper()
                var_base = [a for a in ["A", "C", "G", "T"] if a != ref_base][0]
                pivot_var = f"{self.pivot_site}_{ref_base}_{var_base}"
                self.add_sites.append(pivot_var)

        raw_read_haps = self.get_haplotypes_from_reads(
            check_clip=True,
            kept_sites=homo_sites_to_add,
            add_sites=self.add_sites,
            homo_sites=homo_sites_to_add,
        )

        het_sites = self.het_sites
        known_del = {}
        if self.del1_reads_partial != set():
            raw_read_haps, het_sites = self.update_reads_for_deletions(
                raw_read_haps,
                het_sites,
                self.del1_3p_pos1,
                self.del1_5p_pos2,
                self.del1_reads_partial,
                "3",
                self.deletion1_name,
            )
            known_del.setdefault("3", self.deletion1_name)
        if self.del2_reads_partial != set():
            raw_read_haps, het_sites = self.update_reads_for_deletions(
                raw_read_haps,
                het_sites,
                self.del2_3p_pos1,
                self.del2_5p_pos2,
                self.del2_reads_partial,
                "4",
                self.deletion2_name,
            )
            known_del.setdefault("4", self.deletion2_name)
        self.het_sites = het_sites

        simple_call, phase_result = self.phase_haps_catch_error(raw_read_haps)
        if simple_call is not None:
            return simple_call
        (
            ass_haps,
            original_haps,
            hcn,
            uniquely_supporting_reads,
            nonuniquely_supporting_reads,
            raw_read_haps,
            read_counts,
        ) = phase_result
        tmp = {}
        for i, hap in enumerate(ass_haps):
            tmp.setdefault(hap, f"{self.gene}_hap{i+1}")
        ass_haps = tmp

        haplotypes = None
        if ass_haps != {}:
            haplotypes = self.output_variants_in_haplotypes(
                ass_haps,
                uniquely_supporting_reads,
                nonuniquely_supporting_reads,
                known_del=known_del,
            )

        two_cp_haps = []
        if len(ass_haps) == 1 and self.init_het_sites == []:
            two_cp_haps.append(list(ass_haps.values())[0])
        else:
            if self.targeted:
                if two_cp_haps == []:
                    # check if one haplotype has more reads than others
                    two_cp_haps = self.get_cn2_haplotype(
                        read_counts, ass_haps, min_cn1_read_count=15
                    )
            else:
                if (
                    len(ass_haps) == 3
                    and self.expect_cn2 is False
                    and self.gene != "BPY2"
                ) or (self.gene == "BPY2" and len(ass_haps) < 3):
                    # check if one haplotype has twice the depth of others
                    # at variant sites unique to it
                    two_cp_haps = self.compare_depth(
                        haplotypes, ass_haps, stringent=True
                    )
                    if (
                        two_cp_haps == []
                        and read_counts is not None
                        and len(read_counts) >= 2
                    ):
                        # check if one haplotype has more reads than others
                        two_cp_haps = self.get_cn2_haplotype(read_counts, ass_haps)

        # call fusion
        fusions_called = None
        if self.call_fusion is not None:
            two_cp_haps, fusions_called = self.find_fusion(ass_haps)

        # check gene1 haplotypes and update to cn2 if assume gene1 is never cn1
        # only for targeted mode
        if self.targeted and self.gene1_cn2 and two_cp_haps == []:
            region_length = self.right_boundary - self.left_boundary
            snp_count = len(self.het_sites)
            if snp_count > region_length * 0.008:
                gene1_haps = []
                gene2_haps = []
                for hap_seq, hap_name in ass_haps.items():
                    # this is assuming gene2 is very different from gene1
                    if hap_seq.count("2") > len(hap_seq) * 0.7:
                        gene2_haps.append(hap_name)
                    else:
                        gene1_haps.append(hap_name)
                if len(gene1_haps) == 1:
                    two_cp_haps.append(gene1_haps[0])

        total_cn = len(ass_haps) + len(two_cp_haps)

        # fully homozygous
        if ass_haps == {} and self.init_het_sites == []:
            total_cn = 2

        # two pairs of identical copies
        if (
            total_cn == 2
            and self.expect_cn2 is False
            and self.gene != "BPY2"
            and self.call_fusion is None
        ):
            if self.mdepth is not None:
                prob = self.depth_prob(int(self.region_avg_depth.median), self.mdepth)
                if prob[0] < 0.75:
                    total_cn = 4
                    if two_cp_haps == [] and ass_haps != {}:
                        two_cp_haps = list(ass_haps.values())

        # correct CN for palindrome genes
        # if self.sample_sex is not None:
        #    if self.is_palindrome:
        #        if self.sample_sex == "female" and total_cn < 4:
        #            total_cn = None
        #        elif self.sample_sex == "male" and total_cn < 2:
        #            total_cn = None
        if total_cn is not None and total_cn <= 1:
            total_cn = None

        # phase
        alleles = []
        linked_haps = []
        hap_links = {}
        if self.to_phase is True:
            (
                alleles,
                hap_links,
                _,
                _,
                linked_haps,
            ) = self.phase_alleles(
                uniquely_supporting_reads,
                nonuniquely_supporting_reads,
                raw_read_haps,
                ass_haps,
                reverse=self.is_reverse,
            )
        # all haplotypes are phased into one allele
        if (
            len(alleles) == 1
            and sorted(alleles[0]) == sorted(ass_haps.values())
            and (
                self.sample_sex is not None
                and self.sample_sex == "male"
                and ("X" in self.nchr or "Y" in self.nchr)
            )
            is False
        ):
            alleles = []
            linked_haps = []

        self.close_handle()

        return self.GeneCall(
            total_cn,
            None,
            ass_haps,
            two_cp_haps,
            alleles,
            hap_links,
            hcn,
            original_haps,
            self.het_sites,
            uniquely_supporting_reads,
            self.het_no_phasing,
            self.homo_sites,
            haplotypes,
            nonuniquely_supporting_reads,
            raw_read_haps,
            self.mdepth,
            self.region_avg_depth._asdict(),
            self.sample_sex,
            self.init_het_sites,
            f"{self.genome_build}:{self.nchr}:{self.left_boundary}-{self.right_boundary}",
            linked_haps,
            fusions_called,
        )

    def close_handle(self):
        self._bamh.close()
        self._refh.close()
