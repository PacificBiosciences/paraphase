# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>

import copy
from collections import namedtuple
from ..phaser import Phaser


class Ncf1Phaser(Phaser):
    new_fields = copy.deepcopy(Phaser.fields)
    new_fields.insert(4, "pseudo_reads")
    new_fields.insert(4, "gene_reads")
    new_fields.remove("alleles_final")
    new_fields.remove("hap_links")
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

    def call(self):
        if self.check_coverage_before_analysis() is False:
            return self.GeneCall()
        pivot_site = self.pivot_site
        for pileupcolumn in self._bamh.pileup(
            self.nchr, pivot_site, pivot_site + 1, truncate=True
        ):
            bases = [
                a.upper() for a in pileupcolumn.get_query_sequences(add_indels=True)
            ]
        gene_reads = bases.count("G")
        pseudo_reads = bases.count("G-2NN")

        self.get_homopolymer()
        self.get_candidate_pos()
        self.het_sites = sorted(list(self.candidate_pos))
        self.remove_noisy_sites()

        raw_read_haps = self.get_haplotypes_from_reads(add_sites=["74777265_A_T"])

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

        hap_rename = {}
        counter_gene = 0
        counter_pseudo = 0
        # main variant is 74777266_GGT_G
        pivot_var = "74777266_GGT_G"
        for hap in haplotypes:
            var = haplotypes[hap]["variants"]
            if pivot_var not in var:
                counter_gene += 1
                hap_rename.setdefault(hap, f"ncf1_hap{counter_gene}")
            else:
                counter_pseudo += 1
                hap_rename.setdefault(hap, f"pseudo_hap{counter_pseudo}")

        tmp = {}
        for hap, hap_name in ass_haps.items():
            tmp.setdefault(hap, hap_rename[hap_name])
        ass_haps = tmp

        tmp = {}
        for hap_name, hap_info in haplotypes.items():
            tmp.setdefault(hap_rename[hap_name], hap_info)
        haplotypes = tmp

        two_cp_haps = []
        if counter_gene == 1:
            two_cp_hap_candidate = self.compare_depth(haplotypes)
            if "ncf1_hap1" in two_cp_hap_candidate:
                two_cp_haps = two_cp_hap_candidate
                counter_gene += 1
                total_cn += 1

        if self.mdepth is not None:
            prob = self.depth_prob(gene_reads, self.mdepth / 2)
            if prob[0] < 0.9 and counter_gene == 1:
                counter_gene = None
                total_cn = None
            if prob[0] > 0.95 and counter_gene > 1 and two_cp_haps != []:
                counter_gene = None
                total_cn = None
        # scenario where only three haplotypes are found, possibly each at CN2
        if counter_gene == 1 and counter_pseudo == 2 and total_cn == 3:
            counter_gene = None
            total_cn = None

        # homozygous case
        if total_cn == 0:
            total_cn = None

        self.close_handle()

        return self.GeneCall(
            total_cn,
            counter_gene,
            ass_haps,
            two_cp_haps,
            gene_reads,
            pseudo_reads,
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
