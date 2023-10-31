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
        self, sample_id, outdir, genome_depth=None, genome_bam=None, sample_sex=None
    ):
        Phaser.__init__(self, sample_id, outdir, genome_depth, genome_bam, sample_sex)

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
        count_del = 0
        count_dup = 0
        tmp = {}
        for i, hap in enumerate(ass_haps):
            if "x" in [hap[0], hap[-1]]:
                count_unknown += 1
                new_hap_name = f"hba_unknown_hap{count_unknown}"
            elif hap[0] == "0" and hap[-1] == "0":
                count_hba1 += 1
                new_hap_name = f"hba1_hap{count_hba1}"
            elif hap[0] == "0" and hap[-1] != "0":
                count_dup += 1
                new_hap_name = f"hba_dup_hap{count_dup}"
            elif hap[0] != "0" and hap[-1] == "0":
                count_del += 1
                new_hap_name = f"hba_del_hap{count_del}"
            else:
                count_hba2 += 1
                new_hap_name = f"hba2_hap{count_hba2}"
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
        genotype = None
        for hap in ass_haps.values():
            if "del" in hap or "dup" in hap:
                sv_called.append(hap)
        two_cp_haps = []
        if count_unknown == 0 and count_del == 0:
            if count_hba1 == 1 and count_hba2 == 2:
                two_cp_haps += [a for a in ass_haps.values() if "hba1" in a]
            elif count_hba1 == 2 and count_hba2 == 1:
                two_cp_haps += [a for a in ass_haps.values() if "hba2" in a]
            elif count_hba1 == 1 and count_hba2 == 1:
                # two possible scenarios:
                # two pairs of identical copies, or big deletion on one allele
                # use depth
                if self.mdepth is not None:
                    prob = self.depth_prob(
                        int(surrounding_region_depth), self.mdepth / 2
                    )
                    # print(surrounding_region_depth, self.mdepth / 2, prob)
                    if prob[0] < 0.999:
                        two_cp_haps += [
                            a for a in ass_haps.values() if "hba1" in a or "hba2" in a
                        ]

        total_cn = len(ass_haps) + len(two_cp_haps)
        if self.het_sites == [] or total_cn <= 1:
            total_cn = None
        if count_del == 1 and total_cn == 3:
            genotype = "-a/aa"
        if count_del == 1 and total_cn == 4 and count_dup == 1:
            genotype = "-a/aaa"
        elif count_del == 2 and total_cn == 2:
            genotype = "-a/-a"
        if count_del == 0 and total_cn == 4 and count_dup == 0:
            genotype = "aa/aa"
        if count_del == 0 and total_cn == 5 and count_dup == 1:
            genotype = "aa/aaa"
        if count_del == 0 and total_cn == 6 and count_dup == 2:
            genotype = "aaa/aaa"
        if count_del == 0 and total_cn == 2 and count_dup == 0:
            genotype = "--/aa"

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
