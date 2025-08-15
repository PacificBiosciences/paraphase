# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>

import copy
from collections import namedtuple
from ..phaser import Phaser


class IkbkgPhaser(Phaser):
    new_fields = copy.deepcopy(Phaser.fields)
    new_fields.insert(3, "del_read_number")
    new_fields.insert(3, "deletion_haplotypes")
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
        self.del1_reads = set()
        self.del1_reads_partial = set()

    def set_parameter(self, config):
        super().set_parameter(config)
        self.deletion1_size = config["deletion1_size"]
        self.deletion1_name = config["deletion1_name"]
        self.del1_3p_pos1 = config["del1_3p_pos1"]
        self.del1_3p_pos2 = config["del1_3p_pos2"]
        self.del1_5p_pos1 = config["del1_5p_pos1"]
        self.del1_5p_pos2 = config["del1_5p_pos2"]

    def call(self):
        if self.check_coverage_before_analysis() is False:
            return self.GeneCall(
                genome_depth=self.mdepth,
                region_depth=self.region_avg_depth._asdict(),
                sample_sex=self.sample_sex,
                phase_region=f"{self.genome_build}:{self.nchr}:{self.left_boundary}-{self.right_boundary}",
            )
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
            self.candidate_pos.add(self.add_sites[1])

        self.het_sites = sorted(list(self.candidate_pos))

        kept_sites = [self.add_sites[1]]
        self.het_sites = sorted(self.het_sites)
        self.remove_noisy_sites()
        self.init_het_sites = [a for a in self.het_sites]
        homo_sites_to_add = self.add_homo_sites(min_no_var_region_size=6000)
        # print(homo_sites_to_add)
        kept_sites += homo_sites_to_add

        raw_read_haps = self.get_haplotypes_from_reads(
            check_clip=True,
            partial_deletion_reads=self.del1_reads_partial,
            kept_sites=kept_sites,
            add_sites=[self.add_sites[0]],  # pivot site
            homo_sites=homo_sites_to_add,
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
                self.deletion1_name,
            )
        self.het_sites = het_sites

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

        total_cn = len(ass_haps)
        tmp = {}
        gene_counter = 0
        pseudo_counter = 0
        unknown_counter = 0
        dup_counter = 0
        deletion_haplotypes = []
        pivot_index, index_found = self.get_pivot_site_index()
        if index_found is False:
            return self.GeneCall()
        else:
            for i, hap in enumerate(ass_haps):
                clip_5p = self.get_5pclip_from_hap(hap)
                if clip_5p is None:
                    unknown_counter += 1
                    hap_name = f"{self.gene}_unknownhap{unknown_counter}"
                else:
                    if clip_5p == self.clip_5p_positions[0]:
                        pseudo_counter += 1
                        hap_name = f"{self.gene}_pseudohap{pseudo_counter}"
                    elif clip_5p == self.clip_5p_positions[1]:
                        dup_counter += 1
                        hap_name = f"{self.gene}_duphap{dup_counter}"
                    else:
                        assert clip_5p == 0
                        gene_counter += 1
                        hap_name = f"{self.gene}_ikbkghap{gene_counter}"
                tmp.setdefault(hap, hap_name)
                if "3" in hap:
                    deletion_haplotypes.append(hap_name)

        ass_haps = tmp

        haplotypes = None
        if self.het_sites != []:
            haplotypes = self.output_variants_in_haplotypes(
                ass_haps,
                uniquely_supporting_reads,
                nonuniquely_supporting_reads,
                known_del={
                    "3": self.deletion1_name,
                },
            )

        # this is on chrX, males have one copy of gene and one copy of pseudogene
        if self.sample_sex is not None and self.sample_sex == "male":
            if unknown_counter == 0 and (gene_counter > 1 or pseudo_counter > 1):
                gene_counter = None
                total_cn = None
        two_cp_haps = []
        if (
            unknown_counter == 0
            and self.sample_sex is not None
            and self.sample_sex == "female"
        ):
            if gene_counter == 1 and pseudo_counter == 1:
                two_cp_haps = [a for a in ass_haps.values() if "dup" not in a]
            elif (gene_counter > 1 and pseudo_counter == 1) or (
                gene_counter == 1 and pseudo_counter > 1
            ):
                two_cp_haps = self.compare_depth(haplotypes, ass_haps, loose=True)
                if two_cp_haps == []:
                    # check if one haplotype has more reads than others
                    two_cp_haps = self.get_cn2_haplotype(
                        read_counts, ass_haps, prob_cutoff=0.15
                    )
            for hap in two_cp_haps:
                total_cn += 1
                if "ikbkghap" in hap:
                    gene_counter += 1
            if gene_counter == 1 and pseudo_counter != 1:
                gene_counter = None
                total_cn = None

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
        # all haplotypes are phased into one allele
        if (
            len(alleles) == 1
            and sorted(alleles[0]) == sorted(ass_haps.values())
            and self.sample_sex is not None
            and self.sample_sex == "female"
        ):
            alleles = []
            linked_haps = []

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
            len(self.del1_reads_partial),
            two_cp_haps,
            [],
            hap_links,
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
