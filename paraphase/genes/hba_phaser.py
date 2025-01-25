# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>

import copy
from collections import namedtuple
import pysam
from ..phaser import Phaser


class HbaPhaser(Phaser):
    new_fields = copy.deepcopy(Phaser.fields)
    new_fields.remove("gene_cn")
    new_fields.insert(3, "genotype")
    new_fields.insert(4, "sv_called")
    new_fields.insert(5, "surrounding_region_depth")
    GeneCall = namedtuple(
        "GeneCall",
        new_fields,
        defaults=(None,) * len(new_fields),
    )

    def __init__(
        self,
        sample_id,
        outdir,
        args,
        genome_depth=None,
        genome_bam=None,
        sample_sex=None,
    ):
        Phaser.__init__(
            self, sample_id, outdir, args, genome_depth, genome_bam, sample_sex
        )

    def set_parameter(self, config):
        super().set_parameter(config)
        self.depth_region = config["depth_region"]

    def call(self):
        if self.check_coverage_before_analysis() is False:
            return self.GeneCall()
        genome_bamh = pysam.AlignmentFile(self.genome_bam, "rb")
        surrounding_region_depth = self.get_regional_depth(
            genome_bamh, self.depth_region
        )[0].median
        genome_bamh.close()

        self.get_homopolymer()
        self.get_candidate_pos()
        self.het_sites = sorted(list(self.candidate_pos))
        self.remove_noisy_sites()
        self.init_het_sites = [a for a in self.het_sites]
        raw_read_haps = self.get_haplotypes_from_reads(check_clip=True)

        (
            ass_haps,
            original_haps,
            hcn,
            uniquely_supporting_reads,
            nonuniquely_supporting_reads,
            raw_read_haps,
            read_counts,
        ) = self.phase_haps(raw_read_haps)

        count_hba1 = 0
        count_hba2 = 0
        count_unknown = 0
        count_3p7del = 0
        count_3p7dup = 0
        count_4p2del = 0
        count_4p2dup = 0
        count_homology = 0
        tmp = {}
        for i, hap in enumerate(ass_haps):
            clip_5p = self.get_5pclip_from_hap(hap)
            clip_3p = self.get_3pclip_from_hap(hap)
            if clip_3p is None or clip_5p is None:
                count_unknown += 1
                new_hap_name = f"{self.gene}_unknownhap{count_unknown}"
            elif clip_5p == 0 and clip_3p == 0:
                count_hba2 += 1
                new_hap_name = f"{self.gene}_hba2hap{count_hba2}"
            elif clip_5p == 0:
                if clip_3p == self.clip_3p_positions[0]:
                    count_3p7del += 1
                    new_hap_name = f"{self.gene}_3p7delhap{count_3p7del}"
                elif clip_3p == self.clip_3p_positions[1]:
                    count_4p2dup += 1
                    new_hap_name = f"{self.gene}_4p2duphap{count_4p2dup}"
            elif clip_3p == 0:
                if clip_5p == self.clip_5p_positions[0]:
                    count_3p7dup += 1
                    new_hap_name = f"{self.gene}_3p7duphap{count_3p7dup}"
                elif clip_5p == self.clip_5p_positions[1]:
                    count_4p2del += 1
                    new_hap_name = f"{self.gene}_4p2delhap{count_4p2del}"
            else:
                if (
                    clip_5p == self.clip_5p_positions[0]
                    and clip_3p == self.clip_3p_positions[0]
                ):
                    count_hba1 += 1
                    new_hap_name = f"{self.gene}_hba1hap{count_hba1}"
                if (
                    clip_5p == self.clip_5p_positions[1]
                    and clip_3p == self.clip_3p_positions[1]
                ):
                    count_homology += 1
                    new_hap_name = f"{self.gene}_homologyhap{count_homology}"
            tmp.setdefault(hap, new_hap_name)
        ass_haps = tmp

        haplotypes = None
        if self.het_sites != []:
            haplotypes = self.output_variants_in_haplotypes(
                ass_haps,
                uniquely_supporting_reads,
                nonuniquely_supporting_reads,
            )

        sv_called = []
        for hap in ass_haps.values():
            if "del" in hap or "dup" in hap:
                sv_called.append(hap)

        two_cp_haps = []
        if count_3p7del == 0 and count_4p2del == 1:
            if count_hba1 == 1 and count_hba2 == 1:
                two_cp_haps += [a for a in ass_haps.values() if "hba1" in a]
                count_hba1 += 1
        elif count_3p7del == 0 and count_4p2del == 0:
            if count_hba1 == 1 and count_hba2 == 2:
                two_cp_haps += [a for a in ass_haps.values() if "hba1" in a]
                count_hba1 += 1
            elif count_hba1 == 2 and count_hba2 == 1:
                two_cp_haps += [a for a in ass_haps.values() if "hba2" in a]
                count_hba2 += 1
            elif count_hba1 == 1 and count_hba2 == 1:
                # two possible scenarios:
                # two pairs of identical copies, or big deletion on one allele
                # use depth
                if self.mdepth is None:
                    # assume two pairs of identical copies
                    two_cp_haps += [
                        a for a in ass_haps.values() if "hba1" in a or "hba2" in a
                    ]
                    count_hba1 += 1
                    count_hba2 += 1
                else:
                    prob = self.depth_prob(
                        int(surrounding_region_depth), self.mdepth / 2
                    )
                    # print(surrounding_region_depth, self.mdepth / 2, prob)
                    if prob[0] < 0.999:
                        two_cp_haps += [
                            a for a in ass_haps.values() if "hba1" in a or "hba2" in a
                        ]
                        count_hba1 += 1
                        count_hba2 += 1
        elif (
            count_3p7del == 1
            and count_hba1 == 0
            and count_hba2 == 0
            and count_3p7dup == 0
            and count_4p2del == 0
            and count_4p2dup == 0
        ):
            count_3p7del += 1
            two_cp_haps += [a for a in ass_haps.values() if "del" in a]

        total_cn = (
            len(ass_haps)
            + len(two_cp_haps)
            - count_homology
            - count_4p2del
            - count_unknown
        )
        if self.init_het_sites == []:
            total_cn = 2
        if total_cn <= 1:
            total_cn = None

        genotype = None
        if count_4p2del == 0 and count_4p2dup == 0:
            if count_3p7del == 1 and total_cn == 3:
                genotype = "-a/aa"
            if count_3p7del == 1 and total_cn == 4 and count_3p7dup == 1:
                genotype = "-a/aaa"
            elif count_3p7del == 2 and total_cn == 2:
                genotype = "-a/-a"
            if count_3p7del == 0 and total_cn == 4 and count_3p7dup == 0:
                genotype = "aa/aa"
            if count_3p7del == 0 and total_cn == 5 and count_3p7dup == 1:
                genotype = "aa/aaa"
            if count_3p7del == 0 and total_cn == 6 and count_3p7dup == 2:
                genotype = "aaa/aaa"
            if count_3p7del == 0 and total_cn == 2 and count_3p7dup == 0:
                genotype = "--/aa"
        # 4.2
        elif count_4p2del > 0 and count_4p2dup == 0:
            if count_3p7del == 1 and count_hba1 == 1 and count_hba2 == 0:
                # 3.7/4.2
                genotype = "-a/-a"
            elif count_hba2 == 1 and count_hba1 == 2:
                if count_3p7dup == 0:
                    genotype = "-a/aa"
                else:
                    # anti3.7/4.2
                    genotype = "-a/aaa"
            elif count_hba2 == 0 and count_hba1 == 2:
                genotype = "-a/-a"
        elif count_4p2del == 0 and count_4p2dup > 0:
            # we dont consider anti4.2/anti4.2 for now. Should be very rare
            if count_3p7del == 1 and count_hba1 == 1 and count_hba2 == 1:
                # 3.7/anti4.2
                genotype = "aaa/-a"
            elif count_hba2 == 2 and count_hba1 == 2:
                if count_3p7dup == 0:
                    genotype = "aaa/aa"
                else:
                    # anti3.7/anti4.2
                    genotype = "aaa/aaa"
        elif count_4p2del > 0 and count_4p2dup > 0:
            # 4.2/anti4.2
            if count_hba1 == 2 and count_hba2 == 1:
                genotype = "aaa/-a"

        # phase
        alleles = []
        hap_links = {}
        if self.to_phase is True:
            (
                alleles,
                hap_links,
                _,
                _,
            ) = self.phase_alleles(
                uniquely_supporting_reads,
                nonuniquely_supporting_reads,
                raw_read_haps,
                ass_haps,
                reverse=self.is_reverse,
            )
            # case where all haplotypes are phased into one allele
            if len(alleles) == 1:
                if len(alleles[0]) == len(ass_haps):
                    alleles = []
        self.close_handle()

        return self.GeneCall(
            total_cn,
            ass_haps,
            two_cp_haps,
            genotype,
            sv_called,
            surrounding_region_depth,
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
        )
