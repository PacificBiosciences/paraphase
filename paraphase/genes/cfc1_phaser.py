# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


from ..phaser import Phaser


class Cfc1Phaser(Phaser):
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

    def call(self):
        if self.check_coverage_before_analysis() is False:
            return self.GeneCall(
                genome_depth=self.mdepth,
                region_depth=self.region_avg_depth._asdict(),
                sample_sex=self.sample_sex,
                phase_region=f"{self.genome_build}:{self.nchr}:{self.left_boundary}-{self.right_boundary}",
            )
        self.get_homopolymer()
        self.get_candidate_pos()
        self.het_sites = sorted(list(self.candidate_pos))
        self.remove_noisy_sites()
        self.init_het_sites = [a for a in self.het_sites]
        homo_sites_to_add = self.add_homo_sites()
        raw_read_haps = self.get_haplotypes_from_reads(
            kept_sites=homo_sites_to_add,
            add_sites=self.add_sites,
            homo_sites=homo_sites_to_add,
        )

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

        tmp = {}
        for i, hap in enumerate(ass_haps):
            tmp.setdefault(hap, f"{self.gene}_hap{i+1}")
        ass_haps = tmp

        haplotypes = None
        if self.het_sites != []:
            haplotypes = self.output_variants_in_haplotypes(
                ass_haps,
                uniquely_supporting_reads,
                nonuniquely_supporting_reads,
            )

        two_cp_haps = []
        if len(ass_haps) <= 2:
            two_cp_haps = list(ass_haps.values())
        elif len(ass_haps) == 3:
            two_cp_haps = self.compare_depth(haplotypes, ass_haps, loose=True)
            if two_cp_haps == []:
                # check if one haplotype has more reads than others
                two_cp_haps = self.get_cn2_haplotype(
                    read_counts, ass_haps, prob_cutoff=0.15
                )

        total_cn = len(ass_haps) + len(two_cp_haps)

        if self.init_het_sites == []:
            total_cn = 4

        if total_cn < 4:
            total_cn = None

        self.close_handle()

        return self.GeneCall(
            total_cn,
            None,
            ass_haps,
            two_cp_haps,
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
