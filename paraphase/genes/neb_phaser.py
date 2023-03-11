# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


from collections import namedtuple
from .rccx_phaser import RccxPhaser
from ..phaser import Phaser


class NebPhaser(Phaser):
    GeneCall = namedtuple(
        "GeneCall",
        "total_cn final_haplotypes two_copy_haplotypes alleles_raw alleles_final \
        repeat_name highest_total_cn assembled_haplotypes sites_for_phasing \
        unique_supporting_reads het_sites_not_used_in_phasing homozygous_sites \
        haplotype_details variant_genotypes nonunique_supporting_reads \
        read_details genome_depth",
    )

    def __init__(self, sample_id, outdir, wgs_depth=None):
        Phaser.__init__(self, sample_id, outdir, wgs_depth)

    def set_parameter(self, config):
        super().set_parameter(config)

    def call(self):
        if self.check_coverage_before_analysis() is False:
            return None
        self.get_homopolymer()
        self.get_candidate_pos()
        self.het_sites = sorted(list(self.candidate_pos))
        self.remove_noisy_sites()

        raw_read_haps = self.get_haplotypes_from_reads(self.het_sites)

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
            tmp.setdefault(hap, f"hap{i+1}")
        ass_haps = tmp

        haplotypes = None
        dvar = None
        if self.het_sites != []:
            haplotypes, dvar = self.output_variants_in_haplotypes(
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
        elif len(ass_haps) < 6:
            two_cp_haps = self.compare_depth(haplotypes, loose=True)
            if two_cp_haps == [] and read_counts is not None:
                # check if one haplotype has more reads than others
                haps = list(read_counts.keys())
                counts = list(read_counts.values())
                max_count = max(counts)
                cp2_hap = haps[counts.index(max_count)]
                others_max = sorted(counts, reverse=True)[1]
                probs = self.depth_prob(max_count, others_max)
                if probs[0] < 0.15 and others_max >= 10:
                    two_cp_haps.append(ass_haps[cp2_hap])

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
        new_alleles = []
        if two_cp_haps == []:
            alleles, links = RccxPhaser.get_alleles(uniquely_supporting_reads)
            for allele in alleles:
                new_allele = []
                for hap in allele:
                    new_allele.append(ass_haps[hap])
                new_alleles.append(new_allele)
        total_cn = len(ass_haps) + len(two_cp_haps)
        # incorract phasing suggests haplotypes with cn > 1
        if len(new_alleles) == 1 and sorted(new_alleles[0]) == sorted(
            ass_haps.values()
        ):
            new_alleles = []
            total_cn = None
        elif len(tri1) > 2 or len(tri3) > 2:
            total_cn = None
            new_alleles = []

        self.close_handle()

        return self.GeneCall(
            total_cn,
            ass_haps,
            two_cp_haps,
            alleles,
            new_alleles,
            {"tri1": tri1, "tri2": tri2, "tri3": tri3},
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
