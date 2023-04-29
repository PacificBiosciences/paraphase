# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


from ..phaser import Phaser


class Pms2Phaser(Phaser):
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
        # for distinguishing pms2 from pms2cl
        raw_read_haps = self.get_haplotypes_from_reads(
            check_clip=True, add_sites=["5989137_G_A"]
        )

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
                tmp.setdefault(hap, f"pms2clhap{counter_pseudo}")
            else:
                counter_gene += 1
                tmp.setdefault(hap, f"pms2hap{counter_gene}")
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

        # homozygous case
        if total_cn == 0:
            total_cn = None

        self.close_handle()

        return self.GeneCall(
            total_cn,
            pms2_cn,
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
