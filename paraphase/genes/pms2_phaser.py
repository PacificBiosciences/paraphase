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
        self.find_big_deletion(min_size=2900)

        if self.deletion1_size is not None:
            self.del1_reads, self.del1_reads_partial = self.get_long_del_reads(
                self.del1_3p_pos1,
                self.del1_3p_pos2,
                self.del1_5p_pos1,
                self.del1_5p_pos2,
                self.deletion1_size,
            )
        regions_to_check = []
        if self.del1_reads_partial != set():
            regions_to_check += [
                [self.del1_3p_pos1, self.del1_3p_pos2],
                [self.del1_5p_pos1, self.del1_5p_pos2],
            ]
        self.get_candidate_pos()
        self.het_sites = sorted(list(self.candidate_pos))
        self.remove_noisy_sites()
        # for distinguishing pms2 from pms2cl
        raw_read_haps = self.get_haplotypes_from_reads(add_sites=self.add_sites)

        het_sites = self.het_sites
        known_del = {}
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
            known_del.setdefault("3", self.deletion1_name)
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
        counter_unknown = 0
        pivot_index, index_found = self.get_pivot_site_index()
        if index_found is False:
            return self.GeneCall()
        for hap in ass_haps:
            start_seq = hap[pivot_index:]
            if start_seq.count("2") >= 15:
                counter_pseudo += 1
                hap_name = f"pms2clhap{counter_pseudo}"
            elif start_seq.count("2") <= 5:
                counter_gene += 1
                hap_name = f"pms2hap{counter_gene}"
            else:
                counter_unknown += 1
                hap_name = f"pms2_unknown_hap{counter_gene}"
            tmp.setdefault(hap, hap_name)
        ass_haps = tmp

        haplotypes = None
        if self.het_sites != []:
            haplotypes = self.output_variants_in_haplotypes(
                ass_haps,
                uniquely_supporting_reads,
                nonuniquely_supporting_reads,
                known_del=known_del,
            )

        # two-cp haplotypes
        two_cp_haps = []
        if len(ass_haps) < 4 and counter_unknown == 0:
            if counter_gene == 1 and counter_pseudo == 1:
                two_cp_haps = list(ass_haps.values())
            elif len(ass_haps) == 3 and 2 in [counter_gene, counter_pseudo]:
                two_cp_haps = self.compare_depth(haplotypes, loose=True)

        total_cn = len(ass_haps) + len(two_cp_haps)
        pms2_cn = len(
            [a for a in ass_haps.values() if "cl" not in a and "unknown" not in a]
        ) + len([a for a in two_cp_haps if "cl" not in a and "unknown" not in a])
        # bigger cnvs are not handled here yet
        if pms2_cn != 2 or counter_unknown > 0:
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
            nonuniquely_supporting_reads,
            raw_read_haps,
            self.mdepth,
            self.region_avg_depth._asdict(),
            self.sample_sex,
        )
