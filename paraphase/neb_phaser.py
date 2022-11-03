from pprint import pprint
import numpy as np
import pysam
from collections import namedtuple
from .cyp21_phaser import Cyp21Phaser


class NebPhaser(Cyp21Phaser):
    GeneCall = namedtuple(
        "GeneCall",
        "total_cn final_haplotypes two_copy_haplotypes alleles alleles2 \
        highest_total_cn assembled_haplotypes het_sites \
        unique_supporting_reads het_sites_not_used_in_phasing homozygous_sites \
        haplotype_details variant_genotypes nonunique_supporting_reads \
        read_details genome_depth",
    )

    def __init__(self, sample_id, outdir, wgs_depth=None):
        Cyp21Phaser.__init__(self, sample_id, outdir, wgs_depth)

    def set_parameter(self, config):
        self.homopolymer_file = config["data"]["homopolymer"]
        self.nchr = config["coordinates"]["hg38"]["nchr"]
        self.ref = config["data"]["reference"]
        self._refh = pysam.FastaFile(self.ref)
        self.left_boundary = config["coordinates"]["hg38"]["left_boundary"]
        self.right_boundary = config["coordinates"]["hg38"]["right_boundary"]
        self.pivot_site = config["coordinates"]["hg38"]["pivot_site"]
        self.nchr_old = config["coordinates"]["hg38"]["nchr_old"]
        self.offset = int(self.nchr_old.split("_")[1]) - 1

    def call(self):
        """
        Main function that calls SMN1/SMN2 copy number and variants
        """
        self.get_homopolymer()

        ## get deletion ##

        self.get_candidate_pos()
        # last snp outside of repeat
        # if self.candidate_pos != set():
        #    self.candidate_pos.add("154555882_C_G")

        self.het_sites = sorted(list(self.candidate_pos))
        problematic_sites = []
        # for site in self.het_sites:
        #    if 5980880 < int(site.split("_")[0]) < 5980980:
        #        problematic_sites.append(site)
        for site in problematic_sites:
            self.het_sites.remove(site)

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

        alleles = self.get_alleles(uniquely_supporting_reads)
        new_alleles = []
        for allele in alleles:
            new_allele = []
            for hap in allele:
                new_allele.append(ass_haps[hap])
            new_alleles.append(new_allele)

        self.close_handle()

        return self.GeneCall(
            total_cn,
            ass_haps,
            [],
            alleles,
            new_alleles,
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
