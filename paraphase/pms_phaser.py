from pprint import pprint
import numpy as np
from collections import namedtuple
from .phaser import Phaser


class PmsPhaser(Phaser):
    GeneCall = namedtuple(
        "GeneCall",
        "total_cn final_haplotypes two_copy_haplotypes \
        highest_total_cn assembled_haplotypes sites_for_phasing \
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
        self.het_sites = sorted(list(self.candidate_pos))
        problematic_sites = []
        for site in self.het_sites:
            for region in self.noisy_region:
                if region[0] <= int(site.split("_")[0]) <= region[1]:
                    problematic_sites.append(site)
        for site in problematic_sites:
            self.het_sites.remove(site)
        self.het_sites.append("5989137_G_A")

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
        counter_gene = 1
        counter_pseudo = 1
        for hap in ass_haps:
            if hap[-1] in ["0", "x"]:
                tmp.setdefault(hap, f"pms2cl_hap{counter_pseudo}")
                counter_pseudo += 1
            else:
                tmp.setdefault(hap, f"pms2_hap{counter_gene}")
                counter_gene += 1
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
