import numpy as np
from collections import namedtuple
from .phaser import Phaser


class PmsPhaser(Phaser):
    GeneCall = namedtuple(
        "GeneCall",
        "total_cn pms2_cn final_haplotypes two_copy_haplotypes \
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
        self.remove_noisy_sites()
        # for distinguishing pms2 from pms2cl
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

        tmp = {}
        counter_gene = 0
        counter_pseudo = 0
        for hap in ass_haps:
            if hap[-1] in ["0", "x"]:
                counter_pseudo += 1
                tmp.setdefault(hap, f"pms2cl_hap{counter_pseudo}")
            else:
                counter_gene += 1
                tmp.setdefault(hap, f"pms2_hap{counter_gene}")
        ass_haps = tmp

        haplotypes = None
        dvar = None
        if self.het_sites != []:
            haplotypes, dvar = self.output_variants_in_haplotypes(
                ass_haps,
                uniquely_supporting_reads,
                nonuniquely_supporting_reads,
            )

        # two-cp haplotypes
        two_cp_haps = []
        if len(ass_haps) < 4:
            if counter_gene == 1 and counter_pseudo == 1:
                two_cp_haps = list(ass_haps.values())
            elif len(ass_haps) == 3 and 2 in [counter_gene, counter_pseudo]:
                two_cp_haps = self.compare_depth(haplotypes, loose=True)

        total_cn = len(ass_haps) + len(two_cp_haps)
        pms2_cn = len([a for a in ass_haps.values() if "cl" not in a]) + len(
            [a for a in two_cp_haps if "cl" not in a]
        )
        # bigger cnvs are not handled here yet
        if pms2_cn != 2:
            pms2_cn = None
        self.close_handle()

        return self.GeneCall(
            total_cn,
            pms2_cn,
            ass_haps,
            two_cp_haps,
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
