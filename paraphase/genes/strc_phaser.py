# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>

import copy
from collections import namedtuple
import pysam
from ..phaser import Phaser
from paraphase.prepare_bam_and_vcf import pysam_handle


class StrcPhaser(Phaser):
    new_fields = copy.deepcopy(Phaser.fields)
    new_fields.insert(4, "intergenic_depth")
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
        self.reference_fasta = None
        if args is not None:
            self.reference_fasta = args.reference

    def set_parameter(self, config):
        super().set_parameter(config)
        self.deletion1_size = config["deletion1_size"]
        self.deletion1_name = config["deletion1_name"]
        self.del1_3p_pos1 = config["del1_3p_pos1"]
        self.del1_3p_pos2 = config["del1_3p_pos2"]
        self.del1_5p_pos1 = config["del1_5p_pos1"]
        self.del1_5p_pos2 = config["del1_5p_pos2"]
        self.depth_region = config["depth_region"]

    def call(self):
        if self.check_coverage_before_analysis() is False:
            return self.GeneCall(
                genome_depth=self.mdepth,
                region_depth=self.region_avg_depth._asdict(),
                sample_sex=self.sample_sex,
                phase_region=f"{self.genome_build}:{self.nchr}:{self.left_boundary}-{self.right_boundary}",
            )
        genome_bamh = pysam_handle(self.genome_bam, self.reference_fasta)
        intergenic_depth = self.get_regional_depth(genome_bamh, self.depth_region)[
            0
        ].median
        genome_bamh.close()
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
        self.remove_noisy_sites()
        self.init_het_sites = [a for a in self.het_sites]
        homo_sites_to_add = self.add_homo_sites()
        raw_read_haps = self.get_haplotypes_from_reads(
            kept_sites=homo_sites_to_add,
            homo_sites=homo_sites_to_add,
            add_sites=self.add_sites,
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

        haplotypes = None
        dvar = None
        two_cp_haps = []
        tmp = {}
        counter_gene = 0
        counter_pseudo = 0
        for hap in ass_haps:
            if "3" in hap:
                counter_pseudo += 1
                tmp.setdefault(hap, f"{self.gene}_strcp1hap{counter_pseudo}")
            else:
                counter_gene += 1
                tmp.setdefault(hap, f"{self.gene}_strchap{counter_gene}")
        ass_haps = tmp

        if self.het_sites != []:
            haplotypes = self.output_variants_in_haplotypes(
                ass_haps,
                uniquely_supporting_reads,
                nonuniquely_supporting_reads,
                known_del={"3": self.deletion1_name},
            )
            # one haplotype, identical on both alleles
            if len(ass_haps) == 1 and self.init_het_sites == []:
                two_cp_haps.append(list(ass_haps.values())[0])
            # identify cn2 haplotypes, stringent
            elif counter_gene == 1 or counter_pseudo == 1:
                two_cp_haps = self.compare_depth(haplotypes, ass_haps)
            # two haplotypes, identical on both alleles
            if (
                intergenic_depth > 5
                and counter_gene == 1
                and counter_pseudo == 1
                and two_cp_haps == []
            ):
                two_cp_haps = list(ass_haps.values())
            # identify cn2 haplotypes, loose
            elif two_cp_haps == [] and counter_gene == 1 and counter_pseudo > 1:
                # check if the strc haplotype has more reads than others
                two_cp_haps = self.get_cn2_haplotype(
                    read_counts, ass_haps, prob_cutoff=0.15
                )
                two_cp_haps = [a for a in two_cp_haps if "strcp1" not in a]

            for hap in two_cp_haps:
                if "strcp1" not in hap:
                    counter_gene += 1
                else:
                    counter_pseudo += 1

            total_cn = len(ass_haps) + len(two_cp_haps)

            # check depth between STRC and pseudogene
            if self.mdepth is not None:
                prob = self.depth_prob(int(intergenic_depth), self.mdepth / 2)
                if prob[0] < 0.9 and counter_gene == 1 and counter_pseudo == 2:
                    counter_gene = None
                    total_cn = None
                if prob[0] > 0.95 and counter_gene > 1 and counter_pseudo > 1:
                    counter_gene = None
                    total_cn = None

        # homozygous case
        else:
            total_cn = 2
            # both copies are pseudogene
            if self.del1_reads_partial != set():
                counter_gene = 0
            else:
                counter_gene = 2

        self.close_handle()

        return self.GeneCall(
            total_cn,
            counter_gene,
            ass_haps,
            two_cp_haps,
            intergenic_depth,
            None,
            None,
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
        )
