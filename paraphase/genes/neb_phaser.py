# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


import copy
from collections import namedtuple
from ..phaser import Phaser


class NebPhaser(Phaser):
    new_fields = copy.deepcopy(Phaser.fields)
    new_fields.remove("gene_cn")
    new_fields.insert(5, "repeat_name")
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

    def call(self):
        if self.check_coverage_before_analysis() is False:
            return self.GeneCall(
                genome_depth=self.mdepth,
                region_depth=self.region_avg_depth._asdict(),
                sample_sex=self.sample_sex,
                phase_region=f"{self.genome_build}:{self.nchr}:{self.left_boundary}-{self.right_boundary}",
            )
        self.get_homopolymer()
        self.get_candidate_pos()
        self.het_sites = sorted(list(self.candidate_pos))
        self.remove_noisy_sites()
        self.init_het_sites = [a for a in self.het_sites]
        homo_sites_to_add = self.add_homo_sites()
        raw_read_haps = self.get_haplotypes_from_reads(
            kept_sites=homo_sites_to_add,
            homo_sites=homo_sites_to_add,
        )

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
        if self.het_sites != []:
            haplotypes = self.output_variants_in_haplotypes(
                ass_haps,
                uniquely_supporting_reads,
                nonuniquely_supporting_reads,
            )

        # assign tri 1 2 3
        tri1 = []
        tri2 = []
        tri3 = []
        for hap, hap_name in ass_haps.items():
            nsite = min(int(len(hap) / 2), 10)
            start_seq = hap[:nsite]
            end_seq = hap[0 - nsite :]
            if start_seq.count("1") >= start_seq.count("2") and end_seq.count(
                "1"
            ) >= end_seq.count("2"):
                tri1.append(hap_name)
            elif start_seq.count("1") < start_seq.count("2") and end_seq.count(
                "1"
            ) < end_seq.count("2"):
                tri3.append(hap_name)
            else:
                tri2.append(hap_name)

        # find two copy haplotypes
        two_cp_haps = []
        if len(ass_haps) == 3 and len(tri1) == 1 and len(tri2) == 1 and len(tri3) == 1:
            two_cp_haps = list(ass_haps.values())
        elif len(ass_haps) < 6 and len(ass_haps) > 1:
            two_cp_haps = self.compare_depth(haplotypes, ass_haps, loose=True)
            if two_cp_haps == []:
                # check if one haplotype has more reads than others
                two_cp_haps = self.get_cn2_haplotype(
                    read_counts, ass_haps, prob_cutoff=0.15
                )

        for hap in two_cp_haps:
            if hap in tri1:
                tri1.append(hap)
            elif hap in tri2:
                tri2.append(hap)
            elif hap in tri3:
                tri3.append(hap)
        if len(tri1) == 1:
            two_cp_haps.append(tri1[0])
            tri1.append(tri1[0])
        if len(tri3) == 1:
            two_cp_haps.append(tri3[0])
            tri3.append(tri3[0])

        alleles = []
        linked_haps = []
        hap_links = {}
        if two_cp_haps == []:
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

        total_cn = len(ass_haps) + len(two_cp_haps)
        # incorrect phasing suggests haplotypes with cn > 1
        if len(alleles) == 1 and sorted(alleles[0]) == sorted(ass_haps.values()):
            alleles = []
            linked_haps = []
            total_cn = None
        elif len(tri1) > 2 or len(tri3) > 2:
            total_cn = None
            alleles = []

        # if tri2 has links to both haps in tri1 or tri3
        if two_cp_haps == [] and len(tri2) == 1 and len(tri1) == 2 and len(tri3) == 2:
            tri2_links = hap_links.get(tri2[0])
            if tri2_links is not None:
                if (
                    len(set(tri2_links).intersection(set(tri1))) > 1
                    or len(set(tri2_links).intersection(set(tri3))) > 1
                ):
                    total_cn = 6
                    two_cp_haps = [tri2[0]]
                    alleles = []
                    linked_haps = []

        # homozygous case
        if total_cn == 0:
            total_cn = None

        self.close_handle()

        return self.GeneCall(
            total_cn,
            ass_haps,
            two_cp_haps,
            alleles,
            hap_links,
            {"tri1": tri1, "tri2": tri2, "tri3": tri3},
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
        )
