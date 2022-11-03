from pprint import pprint
import numpy as np
from collections import namedtuple
import re
from .phaser import Phaser


class F8Phaser(Phaser):
    GeneCall = namedtuple(
        "GeneCall",
        "total_cn final_haplotypes two_copy_haplotypes \
        highest_total_cn assembled_haplotypes het_sites \
        unique_supporting_reads het_sites_not_used_in_phasing homozygous_sites \
        haplotype_details variant_genotypes nonunique_supporting_reads \
        read_details genome_depth",
    )

    def __init__(self, sample_id, outdir, wgs_depth=None):
        Phaser.__init__(self, sample_id, outdir, wgs_depth)

    def call(self):
        """
        Main function that calls SMN1/SMN2 copy number and variants
        """
        self.get_homopolymer()
        self.get_candidate_pos()

        var_found = False
        for var in self.candidate_pos:
            pos = int(var.split("_")[0])
            if self.clip_3p_positions[0] < pos < self.clip_3p_positions[1]:
                var_found = True
                break
        if var_found is False:
            self.candidate_pos.add("155386300_A_C")

        var_found = False
        for var in self.candidate_pos:
            pos = int(var.split("_")[0])
            if pos > self.clip_3p_positions[1]:
                var_found = True
                break
        if var_found is False:
            self.candidate_pos.add("155386860_C_G")

        self.het_sites = sorted(list(self.candidate_pos))
        problematic_sites = []
        # for site in self.het_sites:
        #    if 5980880 < int(site.split("_")[0]) < 5980980:
        #        problematic_sites.append(site)
        for site in problematic_sites:
            self.het_sites.remove(site)

        raw_read_haps = self.get_haplotypes_from_reads(self.het_sites, check_clip=True)

        (
            ass_haps,
            original_haps,
            hcn,
            uniquely_supporting_reads,
            nonuniquely_supporting_reads,
            raw_read_haps,
            read_counts,
        ) = self.phase_haps(raw_read_haps)

        total_cn = len(ass_haps)
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

        self.close_handle()

        return self.GeneCall(
            total_cn,
            ass_haps,
            [],
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
