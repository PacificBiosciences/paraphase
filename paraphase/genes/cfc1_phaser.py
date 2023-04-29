# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


from ..phaser import Phaser


class Cfc1Phaser(Phaser):
    def __init__(
        self, sample_id, outdir, genome_depth=None, genome_bam=None, sample_sex=None
    ):
        Phaser.__init__(self, sample_id, outdir, genome_depth, genome_bam, sample_sex)

    def call(self):
        if self.check_coverage_before_analysis() is False:
            return self.GeneCall()
        self.get_homopolymer()
        self.get_candidate_pos()
        self.het_sites = sorted(list(self.candidate_pos))
        self.remove_noisy_sites()

        raw_read_haps = self.get_haplotypes_from_reads(add_sites=["130593061_A_G"])

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

        two_cp_haps = []
        if len(ass_haps) == 2:
            two_cp_haps = list(ass_haps.values())
        elif len(ass_haps) == 3:
            two_cp_haps = self.compare_depth(haplotypes, loose=True)
            if two_cp_haps == [] and read_counts is not None:
                # check if one smn1 haplotype has more reads than others
                haps = list(read_counts.keys())
                counts = list(read_counts.values())
                max_count = max(counts)
                cp2_hap = haps[counts.index(max_count)]
                others_max = sorted(counts, reverse=True)[1]
                probs = self.depth_prob(max_count, others_max)
                if probs[0] < 0.15 and others_max >= 10:
                    two_cp_haps.append(ass_haps[cp2_hap])

        total_cn = len(ass_haps) + len(two_cp_haps)
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
            dvar,
            nonuniquely_supporting_reads,
            raw_read_haps,
            self.mdepth,
        )
