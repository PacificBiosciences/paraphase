# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


import pysam
from pprint import pprint
import os
import numpy as np
import copy
from collections import Counter
import re
from scipy.stats import poisson
from .haplotype_assembler import VariantGraph


class Phaser:

    clip_5p = r"^\d+S|^\d+H"
    clip_3p = r"\d+S$|\d+H$"
    deletion = r"\d+D"

    def __init__(self, sample_id, outdir, config, wgs_depth=None):
        self.outdir = outdir
        self.sample_id = sample_id
        self.bam = os.path.join(outdir, self.sample_id + "_realigned.bam")
        if os.path.exists(self.bam) is False:
            raise Exception(f"File {self.bam} not found.")
        self._bamh = pysam.AlignmentFile(self.bam, "rb")
        self.homopolymer_file = config["data"]["homopolymer"]
        self.homopolymer_sites = {}
        self.het_sites = []  # for phasing
        self.het_no_phasing = []
        self.homo_sites = []
        self.candidate_pos = set()
        self.mdepth = wgs_depth
        self.nchr = config["coordinates"]["hg38"]["nchr"]
        self.ref = config["data"]["reference"]
        self._refh = pysam.FastaFile(self.ref)
        self.left_boundary = config["coordinates"]["hg38"]["left_boundary"]
        self.right_boundary = config["coordinates"]["hg38"]["right_boundary"]
        self.pivot_site = config["coordinates"]["hg38"]["pivot_site"]
        self.nchr_old = config["coordinates"]["hg38"]["nchr_old"]
        self.offset = int(self.nchr_old.split("_")[1]) - 1

    def get_homopolymer(self):
        """Parse the homopolymer site file"""
        with open(self.homopolymer_file) as f:
            for line in f:
                at = line.split()
                self.homopolymer_sites.setdefault(int(at[1]), at[-1])

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

    @staticmethod
    def check_del(read, del_size):
        """Find reads having the 6.3kb deletion in its cigar string"""
        del_len = [int(a[:-1]) for a in re.findall(Phaser.deletion, read.cigarstring)]
        if del_len != [] and abs(max(del_len) - del_size) < 50:
            return True
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
            find_clip_3p = re.findall(self.clip_3p, read.cigarstring)
            if find_clip_3p != [] and pos1 < read.reference_end < pos2:
                if (
                    int(find_clip_3p[0][:-1]) >= min_clip_len
                    and read.reference_start < reference_start_cutoff
                ):
                    p3_reads.add(read.query_name)
            if self.check_del(read, del_size):
                del_reads.add(read.query_name)
        # 5 prime clip
        pos1 = p5_pos1
        pos2 = p5_pos2
        reference_end_cutoff = pos2 + min_extend
        for read in bamh.fetch(self.nchr, pos1, pos2):
            find_clip_5p = re.findall(self.clip_5p, read.cigarstring)
            if find_clip_5p != [] and pos1 < read.reference_start < pos2:
                if (
                    int(find_clip_5p[0][:-1]) >= min_clip_len
                    and read.reference_end > reference_end_cutoff
                ):
                    p5_reads.add(read.query_name)
            if self.check_del(read, del_size):
                del_reads.add(read.query_name)
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

    def get_haplotypes_from_reads(self, het_sites, exclude_reads=[], min_mapq=5):
        """
        Go through reads and get bases at sites of interest
        Returns:
            read_haps (dict of str:list): collapse each read into just the positions
            of interest. 1 corresponds to ref, 2 corresponds to alt
        """
        read_haps = {}
        nvar = len(het_sites)
        for dsnp_index, allele_site in enumerate(het_sites):
            snp_position_gene1, allele1, allele2, *at = allele_site.split("_")
            snp_position = int(snp_position_gene1)
            for pileupcolumn in self._bamh.pileup(
                self.nchr,
                snp_position - 1,
                snp_position,
                truncate=True,
                min_base_quality=30,
            ):
                for read in pileupcolumn.pileups:
                    read_name = read.alignment.query_name
                    if (
                        not read.is_del
                        and not read.is_refskip
                        and not read.alignment.is_secondary
                        and read.alignment.mapping_quality >= min_mapq
                        and read_name not in exclude_reads
                    ):
                        read_seq = read.alignment.query_sequence
                        start_pos = read.query_position
                        end_pos = start_pos + 1
                        if end_pos < len(read_seq):
                            hap = read_seq[start_pos:end_pos]
                            if read_name not in read_haps:
                                read_haps.setdefault(read_name, ["x"] * nvar)
                            if hap.upper() == allele1.upper():
                                read_haps[read_name][dsnp_index] = "1"
                            elif hap.upper() == allele2.upper():
                                read_haps[read_name][dsnp_index] = "2"
        return read_haps

    def allow_del_bases(self, pos):
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
            )
        return ref_seq, var_seq, indel_size

    def get_candidate_pos(self, regions_to_check=[], min_read_support=5, min_vaf=0.11):
        """
        Get all polymorphic sites in the region, update self.candidate_pos
        """
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
            ref_seq = self._refh.fetch(self.nchr_old, offset_pos - 1, offset_pos)

            if total_depth >= min_read_support and (
                del_bases_count < min_read_support or self.allow_del_bases(pos)
            ):
                all_bases = [a for a in all_bases if a != "*"]
                counter = Counter(all_bases)
                bases = counter.most_common(2)
                if (
                    len(counter) >= 2
                    and bases[1][1] >= min_read_support
                    and bases[1][1] / total_depth > min_vaf
                ):
                    var_seq = None
                    found_ref = ref_seq in [a[0] for a in bases]
                    for base in bases:
                        if base[0] != ref_seq:
                            var_seq = base[0]
                            break
                    if found_ref is True and var_seq is not None:
                        # SNV
                        if "-" not in var_seq and "+" not in var_seq:
                            if pos not in self.homopolymer_sites:
                                variants.setdefault(pos, (ref_seq, var_seq))
                            else:
                                prohibited_bases = self.homopolymer_sites[pos].split(
                                    ","
                                )
                                if var_seq not in prohibited_bases:
                                    if "1" in prohibited_bases:
                                        variants.setdefault(pos, (ref_seq, var_seq))
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
                                variants_no_phasing.setdefault(pos, (ref_seq, var_seq))
                elif len(counter) == 1 or (
                    len(counter) >= 2
                    and bases[0][1] > len(all_bases) - min_read_support
                ):
                    var_seq = bases[0][0]
                    if var_seq != ref_seq:
                        # SNV and indels
                        if "-" not in var_seq and "+" not in var_seq:
                            # "homo" sites in large deletions should be put back into het sites
                            if (
                                self.allow_del_bases(pos)
                                and del_bases_count >= min_read_support
                                and pos not in self.homopolymer_sites
                            ):
                                variants.setdefault(pos, (ref_seq, var_seq))
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
        # exclude variants caused by shifted softclips of the big deletions
        excluded_variants = []
        for region in regions_to_check:
            var_to_check = [a for a in variants if region[0] < a < region[1]]
            excluded_variants += var_to_check
        for pos in variants:
            if pos not in excluded_variants:
                ref_seq, var_seq = variants[pos]
                self.candidate_pos.add(f"{pos}_{ref_seq}_{var_seq}")

        excluded_variants = []
        for region in regions_to_check:
            var_to_check = [a for a in variants_no_phasing if region[0] < a < region[1]]
            excluded_variants += var_to_check
        for pos in variants_no_phasing:
            if pos not in excluded_variants:
                ref_seq, var_seq = variants_no_phasing[pos]
                self.het_no_phasing.append(f"{pos}_{ref_seq}_{var_seq}")

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

    def output_variants_in_haplotypes(self, haps, reads, nonunique):
        """
        Summarize all variants in each haplotype.
        Output all variants and their genotypes.
        Haplotypes are different length, so a range (boundary) is reported
        """
        het_sites = self.het_sites
        haplotype_variants = {}
        dvar = {}
        var_no_phasing = copy.deepcopy(self.het_no_phasing)
        for hap_index, hap in enumerate(haps):
            hap_name = f"hap{hap_index}"
            haplotype_variants.setdefault(hap_name, [])
        # het sites not used in phasing
        if reads != {}:
            for var in var_no_phasing:
                genotypes = []
                var_reads = self.check_variants_in_haplotypes(var)
                haps_with_variant = []
                for hap_index, hap in enumerate(haps):
                    hap_name = f"hap{hap_index}"
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
                    dvar.setdefault(var, genotypes)
        # het sites and homo sites
        for hap_index, hap in enumerate(haps):
            hap_name = f"hap{hap_index}"
            for i in range(len(hap)):
                if hap[i] == "2":
                    haplotype_variants[hap_name].append(het_sites[i])
            # need some coordinate check if there is long deletion
            haplotype_variants[hap_name] += self.homo_sites

            var_nstart, var_nend = self.get_hap_variant_ranges(hap)
            var_tmp = haplotype_variants[hap_name]
            var_tmp1 = [
                a for a in var_tmp if var_nstart <= int(a.split("_")[0]) <= var_nend
            ]
            var_tmp1 = list(set(var_tmp1))
            haplotype_variants[hap_name] = sorted(
                var_tmp1,
                key=lambda x: int(x.split("_")[0]),
            )
            haplotype_variants[hap_name].append((var_nstart, var_nend))

        # summary per variant
        all_haps = haps
        nhap = len(all_haps)
        for var in self.homo_sites:
            dvar.setdefault(var, ["1"] * nhap)
        for i, var in enumerate(het_sites):
            dvar.setdefault(var, ["."] * nhap)
            for hap_index in range(len(all_haps)):
                hap = all_haps[hap_index]
                if hap[i] == "2":
                    dvar[var][hap_index] = "1"
                elif hap[i] == "1":
                    dvar[var][hap_index] = "0"

        return haplotype_variants, {
            var: "|".join(dvar[var]) for var in dict(sorted(dvar.items()))
        }

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
                        "x" not in hap[(pos1 - 2) : pos1]
                        and "x" not in hap[(pos1 + 1) : (pos1 + 3)]
                    ):
                        hap[pos1] = "1"
                    raw_read_haps[read] = "".join(hap)
        return raw_read_haps, het_sites

    def get_read_counts(self, uniquely_supporting_haps):
        """
        Get unique supporting read counts for each haplotype
        """
        nhap = len(uniquely_supporting_haps)
        nvar = len(self.het_sites)
        hap_bases = []
        for i in range(nhap):
            hap_bases.append([])
        j = 0
        for hap in uniquely_supporting_haps:
            reads = uniquely_supporting_haps[hap]
            for i in range(nvar):
                bases = [a[i] for a in reads]
                hap_bases[j].append(len(bases) - bases.count("x"))
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
            ass_haps, original_haps, hcn = hap_graph.run(debug=debug, make_plot=debug)

        read_support = VariantGraph.match_reads_and_haplotypes(raw_read_haps, ass_haps)

        uniquely_supporting_haps = read_support.unique
        read_counts = self.get_read_counts(uniquely_supporting_haps)

        uniquely_supporting_reads = {}
        for hap in uniquely_supporting_haps:
            for read_hap in uniquely_supporting_haps[hap]:
                uniquely_supporting_reads.setdefault(hap, [])
                uniquely_supporting_reads[hap] += haplotypes_to_reads[read_hap]
        for hap in uniquely_supporting_haps:
            uniquely_supporting_reads[hap] = list(set(uniquely_supporting_reads[hap]))

        nonuniquely_supporting_reads = {}
        for read in read_support.by_read:
            num_matches = len(read_support.by_read[read])
            if num_matches > 1:
                nonuniquely_supporting_reads.setdefault(
                    read, read_support.by_read[read]
                )
        return (
            ass_haps,
            original_haps,
            hcn,
            uniquely_supporting_reads,
            nonuniquely_supporting_reads,
            raw_read_haps,
            read_counts,
        )

    def close_handle(self):
        self._bamh.close()
        self._refh.close()
