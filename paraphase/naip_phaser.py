from pprint import pprint
import numpy as np
import pysam
from collections import namedtuple
from .cyp21_phaser import Cyp21Phaser


class NaipPhaser(Cyp21Phaser):
    GeneCall = namedtuple(
        "GeneCall",
        "total_cn final_haplotypes two_copy_haplotypes alleles alleles2 \
        del_read_number1 del_read_number2 highest_total_cn assembled_haplotypes het_sites \
        unique_supporting_reads het_sites_not_used_in_phasing homozygous_sites \
        haplotype_details variant_genotypes nonunique_supporting_reads \
        read_details genome_depth",
    )

    def __init__(self, sample_id, outdir, wgs_depth=None):
        Cyp21Phaser.__init__(self, sample_id, outdir, wgs_depth)

    def call(self):
        """
        Main function that calls SMN1/SMN2 copy number and variants
        """
        self.get_homopolymer()

        ## get deletion ##
        self.del1_reads, self.del1_reads_partial = self.get_long_del_reads(
            self.del1_3p_pos1,
            self.del1_3p_pos2,
            self.del1_5p_pos1,
            self.del1_5p_pos2,
            self.deletion1_size,
        )
        self.del2_reads, self.del2_reads_partial = self.get_long_del_reads(
            self.del2_3p_pos1,
            self.del2_3p_pos2,
            self.del2_5p_pos1,
            self.del2_5p_pos2,
            self.deletion2_size,
        )

        regions_to_check = []
        if self.del2_reads_partial != set():
            regions_to_check += [
                [self.del2_3p_pos1, self.del2_3p_pos2],
                [self.del2_5p_pos1, self.del2_5p_pos2],
            ]
        if self.del1_reads_partial != set():
            regions_to_check += [
                [self.del1_3p_pos1, self.del1_3p_pos2],
                [self.del1_5p_pos1, self.del1_5p_pos2],
            ]
        self.get_candidate_pos(min_vaf=0.095, regions_to_check=regions_to_check)

        var_found = False
        for var in self.candidate_pos:
            pos = int(var.split("_")[0])
            if pos < self.clip_5p_positions[0]:
                var_found = True
                break
        if var_found is False:
            self.candidate_pos.add("70960448_T_C")

        het_sites = sorted(list(self.candidate_pos))

        problematic_sites = []
        for site in het_sites:
            for region in self.noisy_region:
                if region[0] <= int(site.split("_")[0]) < region[1]:
                    problematic_sites.append(site)
        for site in problematic_sites:
            het_sites.remove(site)

        raw_read_haps = self.get_haplotypes_from_reads(het_sites, check_clip=True)

        if self.del1_reads_partial != set():
            raw_read_haps, het_sites = self.update_reads_for_deletions(
                raw_read_haps,
                het_sites,
                self.del1_3p_pos1,
                self.del1_5p_pos2,
                self.del1_reads_partial,
                "3",
                "70971390_naip_del1",
            )
        if self.del2_reads_partial != set():
            raw_read_haps, het_sites = self.update_reads_for_deletions(
                raw_read_haps,
                het_sites,
                self.del2_3p_pos1,
                self.del2_5p_pos2,
                self.del2_reads_partial,
                "4",
                "71009640_naip_del2",
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

        total_cn = len(ass_haps)
        tmp = {}
        for i, hap in enumerate(ass_haps):
            hap_name = f"hap{i+1}"
            if i >= 9:
                hap_name = hap_name.replace("1", "I")
            tmp.setdefault(hap, hap_name)
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
            len(self.del1_reads_partial),
            len(self.del2_reads_partial),
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
