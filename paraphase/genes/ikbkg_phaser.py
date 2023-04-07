# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


from collections import namedtuple
from .rccx_phaser import RccxPhaser
from ..phaser import Phaser


class IkbkgPhaser(Phaser):
    GeneCall = namedtuple(
        "GeneCall",
        "total_cn gene_cn final_haplotypes deletion_haplotypes two_copy_haplotypes \
        alleles_raw alleles_final del_read_number highest_total_cn assembled_haplotypes \
        sites_for_phasing unique_supporting_reads het_sites_not_used_in_phasing \
        homozygous_sites haplotype_details variant_genotypes nonunique_supporting_reads \
        read_details genome_depth",
    )

    def __init__(self, sample_id, outdir, wgs_depth=None):
        Phaser.__init__(self, sample_id, outdir, wgs_depth)
        self.del1_reads = set()
        self.del1_reads_partial = set()

    def set_parameter(self, config):
        super().set_parameter(config)
        self.deletion1_size = config["coordinates"]["hg38"]["deletion1_size"]
        self.del1_3p_pos1 = config["coordinates"]["hg38"]["del1_3p_pos1"]
        self.del1_3p_pos2 = config["coordinates"]["hg38"]["del1_3p_pos2"]
        self.del1_5p_pos1 = config["coordinates"]["hg38"]["del1_5p_pos1"]
        self.del1_5p_pos2 = config["coordinates"]["hg38"]["del1_5p_pos2"]

    def call(self):
        if self.check_coverage_before_analysis() is False:
            return None
        self.get_homopolymer()

        ## get deletion ##
        self.del1_reads, self.del1_reads_partial = self.get_long_del_reads(
            self.del1_3p_pos1,
            self.del1_3p_pos2,
            self.del1_5p_pos1,
            self.del1_5p_pos2,
            self.deletion1_size,
        )

        self.get_candidate_pos(min_vaf=0.095)

        # add these sites for duplication/deletion calling
        var_found = False
        for var in self.candidate_pos:
            pos = int(var.split("_")[0])
            if pos > self.clip_3p_positions[0]:
                var_found = True
                break
        if var_found is False:
            self.candidate_pos.add("154569800_T_G")

        var_found = False
        for var in self.candidate_pos:
            pos = int(var.split("_")[0])
            if pos < self.clip_5p_positions[0]:
                var_found = True
                break
        if var_found is False:
            self.candidate_pos.add("154555882_C_G")

        self.het_sites = sorted(list(self.candidate_pos))
        self.remove_noisy_sites()

        raw_read_haps = self.get_haplotypes_from_reads(
            check_clip=True,
            partial_deletion_reads=self.del1_reads_partial,
            kept_sites=["154569800_T_G"],
            add_sites=["154555882_C_G"],
        )
        het_sites = self.het_sites
        if self.del1_reads_partial != set():
            raw_read_haps, het_sites = self.update_reads_for_deletions(
                raw_read_haps,
                het_sites,
                self.del1_3p_pos1,
                self.del1_5p_pos2,
                self.del1_reads_partial,
                "3",
                "154558014_del10806",
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
        gene_counter = 0
        pseudo_counter = 0
        dup_counter = 0
        deletion_haplotypes = []
        pivot_index, index_found = self.get_pivot_site_index()
        if index_found is True:
            for i, hap in enumerate(ass_haps):
                start_seq = hap[: pivot_index + 1]
                if start_seq.startswith("0") is False:
                    if start_seq.count("1") >= start_seq.count("2"):
                        gene_counter += 1
                        hap_name = f"ikbkg_hap{gene_counter}"
                        tmp.setdefault(hap, hap_name)
                        if "3" in hap:
                            deletion_haplotypes.append(hap_name)
                    else:
                        pseudo_counter += 1
                        hap_name = f"pseudo_hap{pseudo_counter}"
                        tmp.setdefault(hap, hap_name)
                        if "3" in hap:
                            deletion_haplotypes.append(hap_name)
                else:
                    dup_counter += 1
                    tmp.setdefault(hap, f"dup_hap{dup_counter}")
        ass_haps = tmp

        haplotypes = None
        dvar = None
        if self.het_sites != []:
            haplotypes, dvar = self.output_variants_in_haplotypes(
                ass_haps,
                uniquely_supporting_reads,
                nonuniquely_supporting_reads,
            )

        # this is on chrX, males have one copy of gene and one copy of pseudogene
        two_cp_haps = []
        if gene_counter == 1 and pseudo_counter > 1:
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
            total_cn += 1
            if "ikbkg" in hap:
                gene_counter += 1
        if gene_counter == 1 and pseudo_counter != 1:
            gene_counter = None
            total_cn = None

        alleles, links = RccxPhaser.get_alleles(uniquely_supporting_reads)
        new_alleles = []
        for allele in alleles:
            new_allele = []
            for hap in allele:
                new_allele.append(ass_haps[hap])
            new_alleles.append(new_allele)

        if gene_counter == 0 or pseudo_counter == 0:
            total_cn = None
            gene_counter = None
        # homozygous case
        if total_cn == 0:
            total_cn = None
            gene_counter = None

        self.close_handle()

        return self.GeneCall(
            total_cn,
            gene_counter,
            ass_haps,
            deletion_haplotypes,
            two_cp_haps,
            alleles,
            new_alleles,
            len(self.del1_reads_partial),
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
