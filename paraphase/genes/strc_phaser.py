from collections import namedtuple
import numpy as np
import pysam
from ..phaser import Phaser


class StrcPhaser(Phaser):
    GeneCall = namedtuple(
        "GeneCall",
        "total_cn gene_cn final_haplotypes two_copy_haplotypes intergenic_depth \
        highest_total_cn assembled_haplotypes sites_for_phasing \
        unique_supporting_reads het_sites_not_used_in_phasing homozygous_sites \
        haplotype_details variant_genotypes nonunique_supporting_reads \
        read_details genome_depth",
    )

    def __init__(self, sample_id, outdir, wgs_depth=None):
        Phaser.__init__(self, sample_id, outdir, wgs_depth)
        self.del1_reads = set()
        self.del1_reads_partial = set()

    def set_parameter(self, config):
        super().set_parameter(config)
        self.mdepth, self.region_depth = self.mdepth
        self.region_depth = self.region_depth[0]
        self.deletion1_size = config["coordinates"]["hg38"]["deletion1_size"]
        self.del1_3p_pos1 = config["coordinates"]["hg38"]["del1_3p_pos1"]
        self.del1_3p_pos2 = config["coordinates"]["hg38"]["del1_3p_pos2"]
        self.del1_5p_pos1 = config["coordinates"]["hg38"]["del1_5p_pos1"]
        self.del1_5p_pos2 = config["coordinates"]["hg38"]["del1_5p_pos2"]
        self.intergenic = config["coordinates"]["hg38"]["depth_region"]

    def call(self):
        """
        Main function that calls SMN1/SMN2 copy number and variants
        """
        self.get_homopolymer()
        self.del1_reads, self.del1_reads_partial = self.get_long_del_reads(
            self.del1_3p_pos1,
            self.del1_3p_pos2,
            self.del1_5p_pos1,
            self.del1_5p_pos2,
            self.deletion1_size,
        )
        self.get_candidate_pos()
        self.het_sites = sorted(list(self.candidate_pos))
        problematic_sites = []
        for site in self.het_sites:
            for region in self.noisy_region:
                if region[0] <= int(site.split("_")[0]) <= region[1]:
                    problematic_sites.append(site)
        for site in problematic_sites:
            self.het_sites.remove(site)

        raw_read_haps = self.get_haplotypes_from_reads(self.het_sites)
        het_sites = self.het_sites
        if self.del1_reads_partial != set():
            raw_read_haps, het_sites = self.update_reads_for_deletions(
                raw_read_haps,
                het_sites,
                self.del1_3p_pos1,
                self.del1_5p_pos2,
                self.del1_reads_partial,
                "3",
                "43602630_del314",
            )
        self.het_sites = het_sites
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
            if "3" in hap:
                counter_pseudo += 1
                tmp.setdefault(hap, f"strcp1_hap{counter_pseudo}")
            else:
                counter_gene += 1
                tmp.setdefault(hap, f"strc_hap{counter_gene}")
        ass_haps = tmp

        haplotypes = None
        dvar = None
        if self.het_sites != []:
            haplotypes, dvar = self.output_variants_in_haplotypes(
                ass_haps,
                uniquely_supporting_reads,
                nonuniquely_supporting_reads,
            )

        two_cp_haps = []
        if self.region_depth > 5 and len(ass_haps) == 2:
            two_cp_haps = list(ass_haps.values())
        elif counter_gene == 1 or counter_pseudo == 1:
            two_cp_haps = self.compare_depth(haplotypes)
        for hap in two_cp_haps:
            if "strcp1" not in hap:
                counter_gene += 1
            else:
                counter_pseudo += 1

        total_cn = len(ass_haps) + len(two_cp_haps)

        # check depth between STRC and pseudogene
        if self.mdepth is not None:
            prob = self.depth_prob(int(self.region_depth), self.mdepth / 2)
            # print(self.sample_id, self.region_depth, self.mdepth, prob[0])
            if prob[0] < 0.9 and counter_gene == 1:
                counter_gene = None
                total_cn = None
            if prob[0] > 0.95 and counter_gene > 1:
                counter_gene = None
                total_cn = None

        self.close_handle()

        return self.GeneCall(
            total_cn,
            counter_gene,
            ass_haps,
            two_cp_haps,
            self.region_depth,
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
