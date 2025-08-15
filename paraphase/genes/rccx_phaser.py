# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>

from collections import namedtuple
from ..phaser import Phaser


class RccxPhaser(Phaser):
    new_fields = (
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
    )
    new_fields += tuple(Phaser.fields[5:])
    GeneCall = namedtuple(
        "GeneCall",
        new_fields,
        defaults=(None,) * len(new_fields),
    )

    def __init__(
        self,
        sample_id,
        outdir,
        args=None,
        genome_depth=None,
        genome_bam=None,
        sample_sex=None,
    ):
        Phaser.__init__(
            self, sample_id, outdir, args, genome_depth, genome_bam, sample_sex
        )
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
        self.white_list = config["white_list"]
        self.deletion1_size = config["deletion1_size"]
        self.deletion2_size = config["deletion2_size"]
        self.deletion1_name = config["deletion1_name"]
        self.deletion2_name = config["deletion2_name"]
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
        hcn,
    ):
        """Update phased alleles"""
        two_cp_haplotypes = []
        successful_phasing = False
        # update the case with homozygous deletion
        if (
            len(final_haps) == 1
            and len(single_copies) == 1
            and self.init_het_sites == []
        ):
            two_cp_haplotypes.append(list(final_haps.values())[0])
            new_alleles = [[single_copies[0]], [single_copies[0]]]
            successful_phasing = True
        # the deletion haplotype will be reported as an allele
        elif len(single_copies) == 1 and len(final_haps) < 5:
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
                two_cp_hap_candidate = self.compare_depth(
                    haplotypes, final_haps, loose=True
                )
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
                    and hcn == len(final_haps)
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
            # if both copies are starting or ending
            if len(allele) == 2:
                hp1, hp2 = allele
                if (hp1 in starting_copies and hp2 in starting_copies) or (
                    hp1 in ending_copies and hp2 in ending_copies
                ):
                    wrong_allele = True
                    break
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
            return self.GeneCall(
                genome_depth=self.mdepth,
                region_depth=self.region_avg_depth._asdict(),
                sample_sex=self.sample_sex,
                phase_region=f"{self.genome_build}:{self.nchr}:{self.left_boundary}-{self.right_boundary}",
            )
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
        self.get_candidate_pos(
            regions_to_check=regions_to_check,
            white_list=self.white_list,
        )

        self.het_sites = sorted(list(self.candidate_pos))
        self.remove_noisy_sites()
        self.init_het_sites = [a for a in self.het_sites]
        homo_sites_to_add = self.add_homo_sites()

        raw_read_haps = self.get_haplotypes_from_reads(
            check_clip=True,
            partial_deletion_reads=self.del1_reads_partial,
            kept_sites=homo_sites_to_add,
            add_sites=self.add_sites,
            homo_sites=homo_sites_to_add,
            multi_allelic_sites=self.white_list,
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
                self.deletion2_name,
            )
        if self.del1_reads_partial != set():
            raw_read_haps, het_sites = self.update_reads_for_deletions(
                raw_read_haps,
                het_sites,
                self.del1_3p_pos1,
                self.del1_5p_pos2,
                self.del1_reads_partial,
                "4",
                self.deletion1_name,
            )
        self.het_sites = het_sites

        # assemble haplotypes
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
            hap_name = f"{self.gene}_hap{i+1}"
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
        if final_haps != {}:
            haplotypes = self.output_variants_in_haplotypes(
                final_haps,
                uniquely_supporting_reads,
                nonuniquely_supporting_reads,
                known_del={
                    "4": self.deletion1_name,
                    "3": self.deletion2_name,
                },
            )

        # phase haplotypes into alleles
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
            final_haps,
            reverse=self.is_reverse,
        )

        successful_phasing, alleles, two_cp_haplotypes = self.update_alleles(
            linked_haps,
            haplotypes,
            final_haps,
            single_copies,
            starting_copies,
            ending_copies,
            hcn,
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
        if ass_haps == [] and self.init_het_sites == []:
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
            nonuniquely_supporting_reads,
            raw_read_haps,
            self.mdepth,
            self.region_avg_depth._asdict(),
            self.sample_sex,
            self.init_het_sites,
            f"{self.genome_build}:{self.nchr}:{self.left_boundary}-{self.right_boundary}",
            alleles,
        )
