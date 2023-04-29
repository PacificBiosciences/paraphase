# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


from collections import namedtuple
import copy
from ..phaser import Phaser


class RccxPhaser(Phaser):
    fields = (
        "total_cn",
        "final_haplotypes",
        "two_copy_haplotypes",
        "starting_hap",
        "ending_hap",
        "deletion_hap",
        "phasing_success",
        "alleles_final",
        "annotated_alleles",
        "hap_variants",
        "hap_links",
        "highest_total_cn",
        "assembled_haplotypes",
        "sites_for_phasing",
        "unique_supporting_reads",
        "het_sites_not_used_in_phasing",
        "homozygous_sites",
        "haplotype_details",
        "variant_genotypes",
        "nonunique_supporting_reads",
        "read_details",
        "genome_depth",
    )
    GeneCall = namedtuple(
        "GeneCall",
        fields,
        defaults=(None,) * len(fields),
    )

    def __init__(
        self, sample_id, outdir, genome_depth=None, genome_bam=None, sample_sex=None
    ):
        Phaser.__init__(self, sample_id, outdir, genome_depth, genome_bam, sample_sex)
        self.has_gene1 = False
        self.has_gene2 = False
        self.gene1_reads = set()
        self.gene2_reads = set()
        self.del1_reads = set()
        self.del1_reads_partial = set()
        self.del2_reads = set()
        self.del2_reads_partial = set()

    def set_parameter(self, config):
        super().set_parameter(config)
        self.variant_def = config["data"]["snp_file"]
        self.known_variants = {}
        with open(self.variant_def) as f:
            for line in f:
                split_line = line.split()
                self.known_variants.setdefault(
                    "_".join([split_line[1], split_line[2], split_line[3]]),
                    split_line[-1],
                )
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
                elif (
                    hap[i] == "3"
                    and "32043718_del120" not in haplotype_variants[hap_name]
                ):
                    haplotype_variants[hap_name].append("32043718_del120")
                elif (
                    hap[i] == "4"
                    and "32017431_del6367" not in haplotype_variants[hap_name]
                ):
                    haplotype_variants[hap_name].append("32017431_del6367")
            if "32017431_del6367" in haplotype_variants[hap_name]:
                pos1 = self.del1_3p_pos1
                pos2 = self.del1_5p_pos2
                for var in self.homo_sites:
                    pos = int(var.split("_")[0])
                    if pos < pos1 or pos > pos2:
                        haplotype_variants[hap_name].append(var)
            elif "32043718_del120" in haplotype_variants[hap_name]:
                pos1 = self.del2_3p_pos1
                pos2 = self.del2_5p_pos2
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

    @staticmethod
    def annotate_var(allele_var):
        """annotate an allele with variants"""
        annotated_allele = None
        if len(allele_var) == 2:
            if [] in allele_var:
                annotated_allele = "WT"
            else:
                tmp = sorted(allele_var, key=lambda x: len(x))
                annotated_allele = ",".join(tmp[0])
        elif len(allele_var) == 1:
            if allele_var == [[]]:
                annotated_allele = "pseudogene_deletion"
            else:
                annotated_allele = "deletion_" + ",".join(allele_var[0])
        elif len(allele_var) == 3:
            tmp = sorted(allele_var, key=lambda x: len(x))
            if tmp[0] == []:
                if tmp[1] == []:
                    annotated_allele = "gene_duplication"
                elif abs(len(tmp[1]) - len(tmp[2])) <= 1 and len(tmp[2]) >= 6:
                    annotated_allele = "pseudogene_duplication"
                else:
                    annotated_allele = "duplication_WT_plus_" + ",".join(tmp[1])
            else:
                if abs(len(tmp[1]) - len(tmp[2])) <= 1 and len(tmp[2]) >= 6:
                    annotated_allele = ",".join(tmp[0]) + "_pseudogene_duplication"
                else:
                    annotated_allele = (
                        "duplication_" + ",".join(tmp[0]) + "_plus_" + ",".join(tmp[1])
                    )

        return annotated_allele

    def annotate_alleles(
        self,
        successful_phasing,
        new_alleles,
        hap_variants,
        ending_copies,
        ass_haps,
        two_cp_haplotypes,
    ):
        """annotate the allele type"""
        annotated_alleles = []
        if successful_phasing:
            allele1 = new_alleles[0]
            allele2 = new_alleles[1]
            allele1_var = [hap_variants[a] for a in allele1]
            allele2_var = [hap_variants[a] for a in allele2]
            for allele_var in [allele1_var, allele2_var]:
                annotated_allele = self.annotate_var(allele_var)
                annotated_alleles.append(annotated_allele)
        else:
            if (
                len(ending_copies) == 2
                and len(ass_haps) == 4
                and two_cp_haplotypes == []
            ):
                for allele_var in [hap_variants[a] for a in ending_copies]:
                    annotated_allele = None
                    if allele_var == []:
                        annotated_allele = "WT"
                    else:
                        annotated_allele = ",".join(allele_var)
                    annotated_alleles.append(annotated_allele)
        return annotated_alleles

    def update_alleles(
        self,
        new_alleles,
        haplotypes,
        final_haps,
        single_copies,
        starting_copies,
        ending_copies,
    ):
        """Update phased alleles"""
        two_cp_haplotypes = []
        successful_phasing = False
        # the deletion haplotype will be reported as an allele
        if len(single_copies) == 1 and len(final_haps) < 5:
            if (
                len(new_alleles) == 1
                and len(new_alleles[0]) == len(final_haps) - 1
                and single_copies not in new_alleles
            ):
                new_alleles.append(single_copies)
            elif (
                len(new_alleles) == 1
                and len(new_alleles[0]) < len(final_haps) - 1
                and single_copies not in new_alleles
                and len(starting_copies) == 1
                and len(ending_copies) == 1
            ):
                new_alleles = []
            if new_alleles == []:
                new_alleles.append(single_copies)
                remaining_hap = [
                    a for a in final_haps.values() if a not in single_copies
                ]
                if len(remaining_hap) == len(final_haps) - 1:
                    new_alleles.append(remaining_hap)
        # deletions on each allele
        elif (
            len(single_copies) == 2
            and len(new_alleles) == 0
            and len(final_haps) == 2
            and len(starting_copies) == 0
            and len(ending_copies) == 0
        ):
            new_alleles = [[single_copies[0]], [single_copies[1]]]
            successful_phasing = True
        elif single_copies == []:
            # homozygous, each haplotype has cn 2
            if (
                len(final_haps) == 2
                and len(ending_copies) == 1
                and len(starting_copies) == 1
                and len(new_alleles) == 1
            ):
                two_cp_haplotypes = list(final_haps.values())
                new_alleles.append(new_alleles[0])
                successful_phasing = True
            # depth-based adjustment when found 3 haplotypes or <2 ending haplotypes
            if haplotypes is not None:
                two_cp_hap_candidate = self.compare_depth(haplotypes, loose=True)
                if len(ending_copies) == 1 and len(starting_copies) == 2:
                    if two_cp_hap_candidate == ending_copies:
                        two_cp_haplotypes = two_cp_hap_candidate
                        if len(final_haps) == 3:
                            new_alleles = [
                                [starting_copies[0], ending_copies[0]],
                                [starting_copies[1], ending_copies[0]],
                            ]
                            successful_phasing = True
                elif len(starting_copies) == 1 and len(ending_copies) == 2:
                    if two_cp_hap_candidate == starting_copies:
                        two_cp_haplotypes = two_cp_hap_candidate
                        if len(final_haps) == 3:
                            new_alleles = [
                                [starting_copies[0], ending_copies[0]],
                                [starting_copies[0], ending_copies[1]],
                            ]
                            successful_phasing = True

            # add missing links when there is no two-cp haplotypes
            if two_cp_haplotypes == [] and len(ending_copies) <= 2:
                # add the missing link in cn=4
                if (
                    len(final_haps) in [3, 4]
                    and len(new_alleles) == 1
                    and len(new_alleles[0]) == 2
                ):
                    remaining_hap = [
                        a for a in final_haps.values() if a not in new_alleles[0]
                    ]
                    if len(remaining_hap) == len(final_haps) - 2:
                        new_alleles.append(remaining_hap)
                # add the missing link in cn=5
                if len(final_haps) == 5:
                    if (
                        len(new_alleles) == 1
                        and len(new_alleles[0]) == 2
                        and (
                            (
                                new_alleles[0][0] in starting_copies
                                and new_alleles[0][1] in ending_copies
                            )
                            or (
                                new_alleles[0][1] in starting_copies
                                and new_alleles[0][0] in ending_copies
                            )
                        )
                    ):
                        remaining_hap = [
                            a for a in final_haps.values() if a not in new_alleles[0]
                        ]
                        if len(remaining_hap) == 3:
                            new_alleles.append(remaining_hap)
                    if len(new_alleles) == 2:
                        allele1 = (
                            new_alleles[0][0] in starting_copies
                            and new_alleles[0][1] in ending_copies
                        ) or (
                            new_alleles[0][1] in starting_copies
                            and new_alleles[0][0] in ending_copies
                        )
                        allele2 = (
                            new_alleles[1][0] in starting_copies
                            and new_alleles[1][1] in ending_copies
                        ) or (
                            new_alleles[1][1] in starting_copies
                            and new_alleles[1][0] in ending_copies
                        )
                        if allele1 is True and allele2 is False:
                            remaining_hap = [
                                a
                                for a in final_haps.values()
                                if a not in new_alleles[0]
                            ]
                            if len(remaining_hap) == 3:
                                new_alleles = [new_alleles[0], remaining_hap]
                        elif allele1 is False and allele2 is True:
                            remaining_hap = [
                                a
                                for a in final_haps.values()
                                if a not in new_alleles[1]
                            ]
                            if len(remaining_hap) == 3:
                                new_alleles = [new_alleles[1], remaining_hap]

        # check wrong phasing
        wrong_allele = False
        for allele in new_alleles:
            hp_set = set(allele)
            for hp in hp_set:
                if (
                    hp in ending_copies
                    and allele.count(hp) > 1
                    and hp not in two_cp_haplotypes
                ):
                    wrong_allele = True
                    break
        if wrong_allele:
            new_alleles = []

        if len(new_alleles) == 2:
            if sorted(new_alleles[0] + new_alleles[1]) == sorted(
                list(final_haps.values())
            ):
                successful_phasing = True

        return successful_phasing, new_alleles, two_cp_haplotypes

    def call(self):
        if self.check_coverage_before_analysis() is False:
            return self.GeneCall()
        self.get_homopolymer()
        self.del2_reads, self.del2_reads_partial = self.get_long_del_reads(
            self.del2_3p_pos1,
            self.del2_3p_pos2,
            self.del2_5p_pos1,
            self.del2_5p_pos2,
            self.deletion2_size,
        )
        self.del1_reads, self.del1_reads_partial = self.get_long_del_reads(
            self.del1_3p_pos1,
            self.del1_3p_pos2,
            self.del1_5p_pos1,
            self.del1_5p_pos2,
            self.deletion1_size,
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

        # add last snp outside of repeat
        var_found = False
        for var in self.candidate_pos:
            pos = int(var.split("_")[0])
            if pos > self.clip_3p_positions[0]:
                var_found = True
                break
        if var_found is False and self.candidate_pos != set():
            self.candidate_pos.add("32046300_G_A")
        # add last snp outside of repeat, 5prime
        var_found = False
        for var in self.candidate_pos:
            pos = int(var.split("_")[0])
            if pos < self.clip_5p_positions[0]:
                var_found = True
                break
        if var_found is False and self.candidate_pos != set():
            self.candidate_pos.add("32013265_A_T")

        self.het_sites = sorted(list(self.candidate_pos))
        self.remove_noisy_sites()

        raw_read_haps = self.get_haplotypes_from_reads(
            check_clip=True,
            partial_deletion_reads=self.del1_reads_partial,
            kept_sites=["32046300_G_A", "32013265_A_T"],
        )

        het_sites = self.het_sites
        if self.del2_reads_partial != set():
            raw_read_haps, het_sites = self.update_reads_for_deletions(
                raw_read_haps,
                het_sites,
                self.del2_3p_pos1,
                self.del2_5p_pos2,
                self.del2_reads_partial,
                "3",
                "32043718_del120",
            )
        if self.del1_reads_partial != set():
            raw_read_haps, het_sites = self.update_reads_for_deletions(
                raw_read_haps,
                het_sites,
                self.del1_3p_pos1,
                self.del1_5p_pos2,
                self.del1_reads_partial,
                "4",
                "32017431_del6367",
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

        tmp = {}
        for i, hap in enumerate(ass_haps):
            hap_name = f"hap{i+1}"
            tmp.setdefault(hap, hap_name)
        final_haps = tmp
        # get haps that extend into tnxb
        ending_copies = [
            final_haps[a]
            for a in ass_haps
            if a[0] not in ["0", "x"] and a[-1] not in ["0", "x"]
        ]
        starting_copies = [
            final_haps[a] for a in ass_haps if a[0] == "0" and a[-1] == "0"
        ]
        single_copies = [
            final_haps[a] for a in ass_haps if a[0] == "0" and a[-1] not in ["0", "x"]
        ]

        haplotypes = None
        dvar = None
        if final_haps != {}:
            haplotypes, dvar = self.output_variants_in_haplotypes(
                final_haps,
                uniquely_supporting_reads,
                nonuniquely_supporting_reads,
            )

        # phase haplotypes into alleles
        (alleles, hap_links, _, _,) = self.phase_alleles(
            uniquely_supporting_reads,
            nonuniquely_supporting_reads,
            raw_read_haps,
            final_haps,
            reverse=self.is_reverse,
        )

        successful_phasing, alleles, two_cp_haplotypes = self.update_alleles(
            alleles,
            haplotypes,
            final_haps,
            single_copies,
            starting_copies,
            ending_copies,
        )

        # annotate haplotypes by checking the diff sites
        # output variants carried by each haplotype
        hap_variants = {}
        if haplotypes is not None:
            for hap, hap_info in haplotypes.items():
                hap_variants.setdefault(hap, [])
                for var in hap_info["variants"]:
                    if var in self.known_variants:
                        hap_variants[hap].append(self.known_variants[var])

        total_cn = len(ass_haps) + len(two_cp_haplotypes)
        if ass_haps == [] and self.het_sites == []:
            # homozygous, feed all reads to call variants
            total_cn = 2
        if total_cn < 2 or len(ending_copies) > 2:
            total_cn = None

        annotated_alleles = self.annotate_alleles(
            successful_phasing,
            alleles,
            hap_variants,
            ending_copies,
            ass_haps,
            two_cp_haplotypes,
        )

        self.close_handle()

        return self.GeneCall(
            total_cn,
            final_haps,
            two_cp_haplotypes,
            starting_copies,
            ending_copies,
            single_copies,
            successful_phasing,
            alleles,
            annotated_alleles,
            hap_variants,
            hap_links,
            hcn,
            original_haps,
            self.het_sites,
            uniquely_supporting_reads,
            self.het_no_phasing,
            self.homo_sites,
            haplotypes,
            dvar,
            nonuniquely_supporting_reads,
            raw_read_haps,
            self.mdepth,
        )
