# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


import pysam
import os
import copy
import numpy as np
from collections import Counter
import re
import logging
from scipy.stats import poisson
from .haplotype_assembler import VariantGraph


class Phaser:

    clip_5p = r"^\d+S|^\d+H"
    clip_3p = r"\d+S$|\d+H$"
    deletion = r"\d+D"

    def __init__(self, sample_id, outdir, wgs_depth=None, genome_bam=None):
        self.outdir = outdir
        self.sample_id = sample_id
        self.homopolymer_sites = {}
        self.het_sites = []  # for phasing
        self.het_no_phasing = []
        self.homo_sites = []
        self.candidate_pos = set()
        self.mdepth = wgs_depth
        self.genome_bam = genome_bam

    def set_parameter(self, config):
        self.gene = config["gene"]
        self.bam = os.path.join(
            self.outdir, self.sample_id + f"_{self.gene}_realigned.bam"
        )
        if os.path.exists(self.bam) is False:
            raise Exception(f"File {self.bam} not found.")
        self._bamh = pysam.AlignmentFile(self.bam, "rb")
        self.homopolymer_file = config["data"]["homopolymer"]
        self.nchr = config["coordinates"]["hg38"]["nchr"]
        self.ref = config["data"]["reference"]
        self._refh = pysam.FastaFile(self.ref)
        self.left_boundary = config["coordinates"]["hg38"]["left_boundary"]
        self.right_boundary = config["coordinates"]["hg38"]["right_boundary"]
        self.pivot_site = None
        if "pivot_site" in config["coordinates"]["hg38"]:
            self.pivot_site = config["coordinates"]["hg38"]["pivot_site"]
        self.nchr_old = config["coordinates"]["hg38"]["nchr_old"]
        self.offset = int(self.nchr_old.split("_")[1]) - 1
        self.use_supplementary = False
        if "use_supplementary" in config:
            self.use_supplementary = config["use_supplementary"]
        self.clip_3p_positions = []
        self.clip_5p_positions = []
        if "clip_3p_positions" in config["coordinates"]["hg38"]:
            self.clip_3p_positions = config["coordinates"]["hg38"]["clip_3p_positions"]
        if "clip_5p_positions" in config["coordinates"]["hg38"]:
            self.clip_5p_positions = config["coordinates"]["hg38"]["clip_5p_positions"]
        self.noisy_region = []
        if "noisy_region" in config["coordinates"]["hg38"]:
            self.noisy_region = config["coordinates"]["hg38"]["noisy_region"]

    def get_regional_depth(self, bam_handle, query_region, ninterval=100):
        """Get depth of the query regions"""
        region_depth = []
        for region in query_region:
            depth = []
            nstep = max(1, int((region[1] - region[0]) / ninterval))
            for pos in range(region[0], region[1], nstep):
                for pileupcolumn in bam_handle.pileup(
                    self.nchr, pos - 1, pos, truncate=True
                ):
                    site_depth = pileupcolumn.get_num_aligned()
                    depth.append(site_depth)
            region_depth.append(np.median(depth))
        return region_depth

    def check_coverage_before_analysis(self):
        """check low coverage regions for enrichment data"""
        region_depth = self.get_regional_depth(
            self._bamh, [[self.left_boundary, self.right_boundary]]
        )[0]
        if np.isnan(region_depth) or region_depth < 10:
            logging.warning(
                "This region does not appear to have coverage. Will not attempt to phase haplotypes."
            )
            return False
        return True

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
            read_name = self.get_read_name(read)
            find_clip_3p = re.findall(self.clip_3p, read.cigarstring)
            if find_clip_3p != [] and pos1 < read.reference_end < pos2:
                if (
                    int(find_clip_3p[0][:-1]) >= min_clip_len
                    and read.reference_start < reference_start_cutoff
                ):
                    p3_reads.add(read_name)
            if self.check_del(read, del_size):
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
            if self.check_del(read, del_size):
                del_reads.add(read_name)
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
        read_name = read.query_name
        if read.is_supplementary and self.use_supplementary:
            read_name = (
                read_name + f"_sup_{read.reference_start}_{read.reference_length}"
            )
        return read_name

    def get_read_names(self, read, partial_deletion_reads):
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

    def get_haplotypes_from_reads(
        self,
        het_sites,
        exclude_reads=[],
        min_mapq=5,
        min_clip_len=50,
        check_clip=False,
        partial_deletion_reads=[],
    ):
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
                min_base_quality=29,
            ):
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

        # for softclips starting at a predefined position, mark sites as 0 instead of x
        if check_clip:
            for dsnp_index, allele_site in enumerate(het_sites):
                snp_position_gene1, allele1, allele2, *at = allele_site.split("_")
                snp_position = int(snp_position_gene1)
                for clip_position in sorted(self.clip_3p_positions):
                    if snp_position > clip_position:
                        for read in self._bamh.fetch(
                            self.nchr, clip_position - 10, clip_position + 10
                        ):
                            read_name = self.get_read_name(read)
                            if read_name in read_haps:
                                if abs(read.reference_end - clip_position) < 20:
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
                            self.nchr, clip_position - 10, clip_position + 10
                        ):
                            read_name = self.get_read_name(read)
                            if read_name in read_haps:
                                if abs(read.reference_start - clip_position) < 20:
                                    find_clip_5p = re.findall(
                                        self.clip_5p, read.cigarstring
                                    )
                                    if (
                                        find_clip_5p != []
                                        and int(find_clip_5p[0][:-1]) >= min_clip_len
                                    ):
                                        read_haps[read_name][dsnp_index] = "0"
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

    def remove_noisy_sites(self):
        """remove variants in predefined noisy sites"""
        problematic_sites = []
        for site in self.het_sites:
            for region in self.noisy_region:
                if region[0] <= int(site.split("_")[0]) <= region[1]:
                    problematic_sites.append(site)
        for site in problematic_sites:
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

    def output_variants_in_haplotypes(self, haps, reads, nonunique, two_cp_haps=[]):
        """
        Summarize all variants in each haplotype.
        Output all variants and their genotypes.
        Haplotypes are different length, so a range (boundary) is reported
        """
        het_sites = self.het_sites
        haplotype_variants = {}
        haplotype_info = {}
        dvar = {}
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
                    dvar.setdefault(var, genotypes)
        # het sites and homo sites
        for hap, hap_name in haps.items():
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
            var_tmp2 = sorted(var_tmp1, key=lambda x: int(x.split("_")[0]))
            haplotype_info.setdefault(
                hap_name, {"variants": var_tmp2, "boundary": [var_nstart, var_nend]}
            )

        # summary per variant
        all_haps = haps
        nhap = len(all_haps)
        for var in self.homo_sites:
            dvar.setdefault(var, ["1"] * nhap)
        for i, var in enumerate(het_sites):
            dvar.setdefault(var, [])
            for hap, hap_name in haps.items():
                base_call = "."
                if hap[i] == "2":
                    base_call = "1"
                elif hap[i] == "1":
                    base_call = "0"
                dvar[var].append(base_call)
                if hap_name in two_cp_haps:
                    dvar[var].append(base_call)

        return haplotype_info, {
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
                    elif hap[pos1 - 1] == "0" and hap[pos1 + 1] == "0":
                        hap[pos1] = "0"
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
        for hap in ass_haps:
            uniquely_supporting_reads.setdefault(hap, [])
        for hap in uniquely_supporting_haps:
            for read_hap in uniquely_supporting_haps[hap]:
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

    def compare_depth(self, haplotypes, loose=False):
        """
        For each haplotype, identify the variants where it's different
        from other haplotypes. Check depth at those variant sites and
        see if the depth suggests twice coverage.
        """
        two_cp_haps = []
        bamh = self._bamh
        boundaries = [haplotypes[a]["boundary"] for a in haplotypes]
        nstart = max([a[0] for a in boundaries])
        nend = min(a[1] for a in boundaries)
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
                    self.nchr, pos - 1, pos, truncate=True, min_base_quality=29
                ):
                    bases = [a.upper() for a in pileupcolumn.get_query_sequences()]
                    base_num = bases.count(hap_base)
                    counts.append([base_num, len(bases) - base_num])

            probs = []
            nsites = len(sites)
            for n1, n2 in counts:
                probs.append(self.depth_prob(n1, n2 / other_cn))
            probs_fil = [a for a in probs if a[0] < 0.25]
            if len(probs_fil) >= nsites * 0.6 and nsites >= 5:
                two_cp_haps.append(hap)
            elif loose is True:
                if len(probs_fil) >= nsites * 0.5 and nsites >= 5:
                    two_cp_haps.append(hap)

        return two_cp_haps

    def close_handle(self):
        self._bamh.close()
        self._refh.close()
