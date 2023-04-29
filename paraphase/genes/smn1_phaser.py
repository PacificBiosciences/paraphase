# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


from collections import namedtuple
import copy
import json
from ..phaser import Phaser


class Smn1Phaser(Phaser):
    fields = (
        "smn1_cn",
        "smn2_cn",
        "smn2_del78_cn",
        "smn1_read_number",
        "smn2_read_number",
        "smn2_del78_read_number",
        "highest_total_cn",
        "smn1_haplotypes",
        "smn2_haplotypes",
        "smn2_del78_haplotypes",
        "assembled_haplotypes",
        "two_copy_haplotypes",
        "sites_for_phasing",
        "unique_supporting_reads",
        "het_sites_not_used_in_phasing",
        "homozygous_sites",
        "haplotype_details",
        "variant_genotypes",
        "nonunique_supporting_reads",
        "read_details",
        "final_haplotypes",
        "genome_depth",
    )
    GeneCall = namedtuple(
        "GeneCall",
        fields,
        defaults=(None,) * len(fields),
    )
    HaplotypeInfo = namedtuple("HaplotypeInfo", "variants boundary haplogroup")

    def __init__(
        self, sample_id, outdir, genome_depth=None, genome_bam=None, sample_sex=None
    ):
        Phaser.__init__(self, sample_id, outdir, genome_depth, genome_bam, sample_sex)
        self.has_smn1 = False
        self.has_smn2 = False
        self.smn1_reads = set()
        self.smn2_reads = set()
        self.smn1_reads_splice = None
        self.smn2_reads_splice = None
        self.smn1_del_reads = set()
        self.smn1_del_reads_partial = set()
        self.smn2_del_reads = set()
        self.smn2_del_reads_partial = set()

    def set_parameter(self, config):
        super().set_parameter(config)
        self.deletion1_size = config["deletion1_size"]
        self.deletion2_size = config["deletion2_size"]
        self.del2_3p_pos1 = config["del2_3p_pos1"]
        self.del2_3p_pos2 = config["del2_3p_pos2"]
        self.del2_5p_pos1 = config["del2_5p_pos1"]
        self.del2_5p_pos2 = config["del2_5p_pos2"]
        self.del1_3p_pos1 = config["del1_3p_pos1"]
        self.del1_3p_pos2 = config["del1_3p_pos2"]
        self.del1_5p_pos1 = config["del1_5p_pos1"]
        self.del1_5p_pos2 = config["del1_5p_pos2"]
        self.strip_c_region_start = config["strip_c_region_start"]
        self.strip_c_region_end = config["strip_c_region_end"]
        with open(config["data"]["known_haplotypes"]) as f:
            self.known_haps = json.load(f)

    def check_smn1_smn2_presence(self):
        """
        Get number of reads at splice variant sites and
        check if there is any SMN1/SMN2
        """
        bamh = self._bamh
        pivot_site = self.pivot_site
        for pileupcolumn in bamh.pileup(
            self.nchr, pivot_site - 1, pivot_site, truncate=True
        ):
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    read_seq = pileupread.alignment.query_sequence
                    read_pos = pileupread.query_position
                    base1 = read_seq[read_pos].upper()
                    if base1 == "C":
                        self.smn1_reads.add(pileupread.alignment.query_name)
                    elif base1 == "T" and read_seq[read_pos : read_pos + 2] != "TC":
                        self.smn2_reads.add(pileupread.alignment.query_name)

        self.smn1_reads_splice = len(self.smn1_reads)
        self.smn2_reads_splice = len(self.smn2_reads)
        if self.smn1_reads_splice >= 2:
            self.has_smn1 = True
        if len(self.smn2_del_reads_partial) + self.smn2_reads_splice >= 2:
            self.has_smn2 = True

    def allow_del_bases(self, pos):
        """
        During variant calling, allow calling bases at positions within big
        deletions
        """
        if (
            self.smn2_del_reads_partial != set()
            and self.del2_3p_pos1 <= pos <= self.del2_5p_pos2
        ):
            return True
        if (
            self.smn1_del_reads_partial != set()
            and self.del1_3p_pos1 <= pos <= self.del1_5p_pos2
        ):
            return True
        return False

    def output_variants_in_haplotypes(
        self, smn1_haps, smn2_haps, smn2_del_haps, reads, nonunique, two_cp_haps
    ):
        """
        Summarize all variants in each haplotype.
        Output all variants and their genotypes.
        Need to consider that haplotypes are different length
        """
        het_sites = self.het_sites
        haplotype_info = {}
        haplotype_variants = {}
        dvar = {}
        var_no_phasing = copy.deepcopy(self.het_no_phasing)
        for smn_haps in [smn1_haps, smn2_haps, smn2_del_haps]:
            for hap, hap_name in smn_haps.items():
                haplotype_variants.setdefault(hap_name, [])
        # het sites not used in phasing
        if reads != {}:
            for var in var_no_phasing:
                genotypes = []
                var_reads = self.check_variants_in_haplotypes(var)
                haps_with_variant = []
                for smn_haps in [smn1_haps, smn2_haps, smn2_del_haps]:
                    for hap, hap_name in smn_haps.items():
                        hap_reads = reads[hap]
                        hap_reads_nonunique = [
                            a for a in nonunique if hap in nonunique[a]
                        ]
                        genotype = self.get_genotype_in_hap(
                            var_reads, hap_reads, hap_reads_nonunique
                        )
                        genotypes.append(genotype)
                        if hap_name in two_cp_haps:
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
        for hap, hap_name in smn1_haps.items():
            for i in range(len(hap)):
                if hap[i] == "2":
                    haplotype_variants[hap_name].append(het_sites[i])
                elif (
                    hap[i] == "4"
                    and "70917700_smn1_del" not in haplotype_variants[hap_name]
                ):
                    haplotype_variants[hap_name].append("70917700_smn1_del")
            if "70917700_smn1_del" in haplotype_variants[hap_name]:
                pos1 = self.del1_3p_pos1
                pos2 = self.del1_5p_pos2
                for var in self.homo_sites:
                    pos = int(var.split("_")[0])
                    if pos < pos1 or pos > pos2:
                        haplotype_variants[hap_name].append(var)
            else:
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
        for hap, hap_name in smn2_haps.items():
            for i in range(len(hap)):
                if hap[i] == "2":
                    haplotype_variants[hap_name].append(het_sites[i])
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
        for hap, hap_name in smn2_del_haps.items():
            for i in range(len(hap)):
                if hap[i] == "2":
                    haplotype_variants[hap_name].append(het_sites[i])
                elif (
                    hap[i] == "3"
                    and "70948286_smn2_del" not in haplotype_variants[hap_name]
                ):
                    haplotype_variants[hap_name].append("70948286_smn2_del")
            pos1 = self.del2_3p_pos1
            pos2 = self.del2_5p_pos2
            for var in self.homo_sites:
                pos = int(var.split("_")[0])
                if pos < pos1 or pos > pos2:
                    haplotype_variants[hap_name].append(var)

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

        # assign haplogroups
        for hap_name in haplotype_variants:
            haplogroup = None
            if "del" not in hap_name:
                query_hap = haplotype_variants[hap_name]
                gene = hap_name[:4]
                haplogroup, candidates1, candidates2 = self.assign_hap_to_group(
                    query_hap, gene
                )
            else:
                haplogroup = "smn2_del_exon78"
            haplotype_info.setdefault(
                hap_name,
                self.HaplotypeInfo(
                    haplotype_variants[hap_name][:-1],
                    haplotype_variants[hap_name][-1],
                    haplogroup,
                )._asdict(),
            )

        # summarize variants
        all_haps = (
            list(smn1_haps.keys()) + list(smn2_haps.keys()) + list(smn2_del_haps.keys())
        )
        nhap = len(all_haps) + len(two_cp_haps)
        for var in self.homo_sites:
            dvar.setdefault(var, ["1"] * nhap)
        for i, var in enumerate(het_sites):
            dvar.setdefault(var, [])
            for smn_haps in [smn1_haps, smn2_haps, smn2_del_haps]:
                for hap, hap_name in smn_haps.items():
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

    def assign_haps_to_gene(self, ass_haps):
        """Assign assembled haplotypes to smn1/smn2/smn2del"""
        splice_index, found_splice = self.get_pivot_site_index()
        smn1_haps = []
        smn2_haps = []
        smn2_del_haps = []

        for hap in ass_haps:
            if "3" in hap:
                smn2_del_haps.append(hap)

        if found_splice is True:
            for a in ass_haps:
                if a[splice_index] == "1":
                    smn1_haps.append(a)
                elif "3" not in a:
                    smn2_haps.append(a)
        elif self.has_smn1 is False or self.smn2_reads_splice < 2:
            for a in ass_haps:
                if "3" not in a:
                    if self.has_smn1 is True:
                        smn1_haps.append(a)
                    else:
                        smn2_haps.append(a)
        return smn1_haps, smn2_haps, smn2_del_haps

    def adjust_smn1_cn(self, smn1_cn, smn2_cn, hcn, ass_haps, read_counts, smn1_haps):
        """
        Adjust SMN1 CN if there is only one haplotype found.
        """
        new_smn1_cn = smn1_cn
        two_cp_haps = []
        if smn1_cn is None:
            return None, two_cp_haps
        if self.mdepth is not None:
            if self.mdepth < 20 and smn1_cn == 1:
                return None, two_cp_haps
            genome_depth = self.mdepth
            haploid_depth = genome_depth / 2
            depth1 = self.smn1_reads_splice
            copy_number_probs = self.depth_prob(depth1, haploid_depth)
            if smn1_cn == 1:
                # when there is only one haplotype, check if depth is
                # consistent with haploid depth
                copy_one_prob = copy_number_probs[0]
                if copy_one_prob < 0.05:
                    two_cp_haps = list(smn1_haps.values())
                    return 2, two_cp_haps
                if copy_one_prob < 0.25:
                    new_smn1_cn = None
            elif smn1_cn == 2:
                # scenario where two two-copy alleles are identical
                copy_four_prob = copy_number_probs[3]
                if copy_four_prob > 0.75:
                    two_cp_haps = list(smn1_haps.values())
                    return 4, two_cp_haps
            if smn1_cn in [2, 3] and read_counts is not None:
                # check if one smn1 haplotype has more reads than others
                haps = list(read_counts.keys())
                counts = list(read_counts.values())
                max_count = max(counts)
                cp2_hap = haps[counts.index(max_count)]
                if cp2_hap in smn1_haps:
                    others_max = sorted(counts, reverse=True)[1]
                    probs = self.depth_prob(max_count, others_max)
                    if probs[0] < 0.0005 and others_max >= 10:
                        two_cp_haps.append(smn1_haps[cp2_hap])
                        return smn1_cn + 1, two_cp_haps

        if smn1_cn == 1 and smn2_cn >= 1:
            haploid_depth = self.smn2_reads_splice / smn2_cn
            if haploid_depth >= 10:
                depth1 = self.smn1_reads_splice
                copy_number_probs = self.depth_prob(depth1, haploid_depth)
                # when there is only one haplotype, check if depth is
                # consistent with haploid depth
                copy_one_prob = copy_number_probs[0]
                if copy_one_prob < 0.25:
                    new_smn1_cn = None
                    if smn2_cn > 1 and copy_one_prob < 0.05:
                        two_cp_haps = list(smn1_haps.values())
                        return 2, two_cp_haps
        # if we see more haplotypes in other regions of the gene
        if smn1_cn == 1 and hcn > len(ass_haps):
            if self.has_smn2 is False:
                two_cp_haps = list(smn1_haps.values())
                return 2, two_cp_haps
            else:
                new_smn1_cn = None
        return new_smn1_cn, two_cp_haps

    def adjust_smn2_cn(self, smn1_cn, smn2_cn, smn2_haps):
        """
        Adjust SMN2 CN if there is only one haplotype found.
        """
        new_smn2_cn = smn2_cn
        two_cp_haps = []
        if smn2_cn is None:
            return None, two_cp_haps
        if self.mdepth is not None:
            if self.mdepth < 20 and smn2_cn == 1:
                return None, two_cp_haps
            genome_depth = self.mdepth
            haploid_depth = genome_depth / 2
            depth1 = self.smn2_reads_splice
            copy_number_probs = self.depth_prob(depth1, haploid_depth)
            # if smn1 cn is 3, then no need to adjust smn2 cn
            if (
                smn2_cn == 1
                and len(self.smn2_del_reads_partial) <= 1
                and (smn1_cn is None or smn1_cn <= 2)
            ):
                # when there is only one haplotype, check if depth is
                # consistent with haploid depth
                copy_one_prob = copy_number_probs[0]
                if copy_one_prob < 0.25:
                    new_smn2_cn = None
                    if copy_one_prob < 0.05:
                        two_cp_haps = list(smn2_haps.values())
                        return 2, two_cp_haps
        if smn1_cn is not None:
            if smn2_cn == 1 and smn1_cn == 2 and len(self.smn2_del_reads_partial) <= 1:
                haploid_depth = self.smn1_reads_splice / smn1_cn
                if haploid_depth > 10:
                    depth1 = self.smn2_reads_splice
                    copy_number_probs = self.depth_prob(depth1, haploid_depth)
                    # when there is only one haplotype, check if depth is
                    # consistent with haploid depth
                    copy_one_prob = copy_number_probs[0]
                    if copy_one_prob < 0.25:
                        new_smn2_cn = None
                        if copy_one_prob < 0.05:
                            two_cp_haps = list(smn2_haps.values())
                            return 2, two_cp_haps
        return new_smn2_cn, two_cp_haps

    def get_best_match(
        self,
        known_haps,
        query_hap,
        min_overlap_len=12000,
        strip_c=False,
        max_mismatch=15,
    ):
        """
        Get best match to a haplogroup
        """
        vars_all = query_hap[:-1]
        nstart, nend = query_hap[-1]
        haplogroup_candidates = {}
        for known_hap_read in known_haps:
            known_hap_vars_all = known_haps[known_hap_read][:-2]
            known_hap_nstart, known_hap_nend = known_haps[known_hap_read][-2]
            new_nstart = max(nstart, known_hap_nstart)
            new_nend = min(nend, known_hap_nend)
            if strip_c:
                new_nstart = max(new_nstart, self.strip_c_region_start)
                new_nend = min(new_nend, self.strip_c_region_end)
            if new_nend - new_nstart >= min_overlap_len:
                haplogroup = known_haps[known_hap_read][-1]
                if strip_c:
                    haplogroup = haplogroup.strip("c")
                var1 = set(
                    [
                        a
                        for a in vars_all
                        if new_nstart <= int(a.split("_")[0]) <= new_nend
                    ]
                )
                var2 = set(
                    [
                        a
                        for a in known_hap_vars_all
                        if new_nstart <= int(a.split("_")[0]) <= new_nend
                    ]
                )
                nmismatch = len(var1 - var2) + len(var2 - var1)
                haplogroup_candidates.setdefault(haplogroup, []).append(nmismatch)
        if haplogroup_candidates == {}:
            return None, {}
        tmp = {k: sorted(v) for k, v in haplogroup_candidates.items()}
        haplogroup_candidates = dict(sorted(tmp.items(), key=lambda item: min(item[1])))
        best_two_mismatches = list(haplogroup_candidates.values())[:2]
        best_match = None
        if (
            len(best_two_mismatches[0]) >= 2
            and (
                best_two_mismatches[0][1] <= max_mismatch
                or best_two_mismatches[0][0] <= 2
            )
            and best_two_mismatches[0][0] < best_two_mismatches[1][0] - 2
            and best_two_mismatches[0][1] < best_two_mismatches[1][0]
        ):
            best_match = list(haplogroup_candidates.keys())[0]
        return best_match, haplogroup_candidates

    def assign_hap_to_group(self, query_hap, gene):
        """
        Assign haplotype to a haplogroup
        """
        known_haps = self.known_haps[gene]
        best_match1, candidates1 = self.get_best_match(known_haps, query_hap)
        best_match2, candidates2 = self.get_best_match(
            known_haps, query_hap, strip_c=True, max_mismatch=10
        )
        if best_match1 is not None:
            best_match = best_match1
        elif best_match2 is not None:
            best_match = best_match2
        else:
            best_match = None
        if best_match is not None and best_match == "S1-9":
            if "70917700_smn1_del" in query_hap:
                best_match = "S1-9d"
        return best_match, candidates1, candidates2

    def call(self):
        """
        Main function that calls SMN1/SMN2 copy number and variants
        """
        if self.check_coverage_before_analysis() is False:
            return self.GeneCall()
        self.get_homopolymer()
        # find known deletions
        self.smn2_del_reads, self.smn2_del_reads_partial = self.get_long_del_reads(
            self.del2_3p_pos1,
            self.del2_3p_pos2,
            self.del2_5p_pos1,
            self.del2_5p_pos2,
            self.deletion2_size,
        )
        self.smn1_del_reads, self.smn1_del_reads_partial = self.get_long_del_reads(
            self.del1_3p_pos1,
            self.del1_3p_pos2,
            self.del1_5p_pos1,
            self.del1_5p_pos2,
            self.deletion1_size,
        )
        self.check_smn1_smn2_presence()
        if self.has_smn1 is False and self.has_smn2 is False:
            raise Exception("Not enough reads to analyze this region.")

        # scan for polymorphic sites
        regions_to_check = []
        if self.smn2_del_reads_partial != set():
            regions_to_check += [
                [self.del2_3p_pos1, self.del2_3p_pos2],
                [self.del2_5p_pos1, self.del2_5p_pos2],
            ]
        if self.smn1_del_reads_partial != set():
            regions_to_check += [
                [self.del1_3p_pos1, self.del1_3p_pos2],
                [self.del1_5p_pos1, self.del1_5p_pos2],
            ]
        self.get_candidate_pos(regions_to_check=regions_to_check)

        self.het_sites = sorted(list(self.candidate_pos))
        self.remove_noisy_sites()
        raw_read_haps = self.get_haplotypes_from_reads(add_sites=["70951946_C_T"])

        # update reads for those overlapping known deletions
        het_sites = self.het_sites
        if len(self.smn2_del_reads_partial) > 1:
            raw_read_haps, het_sites = self.update_reads_for_deletions(
                raw_read_haps,
                het_sites,
                self.del2_3p_pos1,
                self.del2_5p_pos2,
                self.smn2_del_reads_partial,
                "3",
                "70948286_smn2_del",
            )
        if len(self.smn1_del_reads_partial) > 1:
            raw_read_haps, het_sites = self.update_reads_for_deletions(
                raw_read_haps,
                het_sites,
                self.del1_3p_pos1,
                self.del1_5p_pos2,
                self.smn1_del_reads_partial,
                "4",
                "70917700_smn1_del",
            )
        self.het_sites = het_sites

        # assemble haplotypes
        (
            ass_haps,
            original_haps,
            hcn,
            uniquely_supporting_reads,
            nonuniquely_supporting_reads,
            raw_read_haps,
            read_counts,
        ) = self.phase_haps(raw_read_haps)

        # process and assign haplotypes
        smn1_cn = None
        smn2_cn = None
        smn2_del_cn = 0
        smn1_haps, smn2_haps, smn2_del_haps = self.assign_haps_to_gene(ass_haps)
        if ass_haps == [] and self.het_sites == []:
            if self.has_smn1 is True and self.has_smn2 is False:
                smn1_cn = 1
                smn2_cn = 0
            elif self.has_smn1 is False and self.has_smn2 is True:
                smn1_cn = 0
                smn2_cn = 1
        else:
            smn1_cn = len(smn1_haps)
            smn2_cn = len(smn2_haps)
            smn2_del_cn = len(smn2_del_haps)

        tmp = {}
        for i, hap in enumerate(smn1_haps):
            tmp.setdefault(hap, f"smn1hap{i+1}")
        smn1_haps = tmp
        tmp = {}
        for i, hap in enumerate(smn2_haps):
            tmp.setdefault(hap, f"smn2hap{i+1}")
        smn2_haps = tmp
        tmp = {}
        for i, hap in enumerate(smn2_del_haps):
            tmp.setdefault(hap, f"smn2del78hap{i+1}")
        smn2_del_haps = tmp

        # update cn based on depth
        smn1_cn_old = smn1_cn
        smn1_cn, two_cp_haps = self.adjust_smn1_cn(
            smn1_cn, smn2_cn, hcn, ass_haps, read_counts, smn1_haps
        )
        smn2_cn, two_cp_haps_smn2 = self.adjust_smn2_cn(smn1_cn_old, smn2_cn, smn2_haps)
        two_cp_haps += two_cp_haps_smn2

        # summarize variants
        haplotypes = None
        dvar = None
        if self.het_sites != []:
            haplotypes, dvar = self.output_variants_in_haplotypes(
                smn1_haps,
                smn2_haps,
                smn2_del_haps,
                uniquely_supporting_reads,
                nonuniquely_supporting_reads,
                two_cp_haps,
            )

        self.close_handle()

        return self.GeneCall(
            smn1_cn,
            smn2_cn,
            smn2_del_cn,
            self.smn1_reads_splice,
            self.smn2_reads_splice,
            len(self.smn2_del_reads_partial),
            hcn,
            smn1_haps,
            smn2_haps,
            smn2_del_haps,
            original_haps,
            two_cp_haps,
            self.het_sites,
            uniquely_supporting_reads,
            self.het_no_phasing,
            self.homo_sites,
            haplotypes,
            dvar,
            nonuniquely_supporting_reads,
            raw_read_haps,
            {**smn1_haps, **smn2_haps, **smn2_del_haps},
            self.mdepth,
        )
