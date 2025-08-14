# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


from ..phaser import Phaser


class Pms2Phaser(Phaser):
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
        self.init_het_sites = [a for a in self.het_sites]
        homo_sites_to_add = self.add_homo_sites(min_no_var_region_size=6000)
        # for distinguishing pms2 from pms2cl
        raw_read_haps = self.get_haplotypes_from_reads(
            clip_buffer=50,
            check_clip=True,
            kept_sites=homo_sites_to_add,
            add_sites=self.add_sites,
            homo_sites=homo_sites_to_add,
        )

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
        counter_gene = 0
        counter_pseudo = 0
        counter_unknown = 0
        for hap in ass_haps:
            clip_position = self.get_3pclip_from_hap(hap)
            if clip_position is None:
                counter_unknown += 1
                hap_name = f"{self.gene}_unknownhap{counter_unknown}"
            elif clip_position in self.clip_3p_positions:
                counter_pseudo += 1
                hap_name = f"{self.gene}_pms2clhap{counter_pseudo}"
            else:
                assert clip_position == 0
                counter_gene += 1
                hap_name = f"{self.gene}_pms2hap{counter_gene}"
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
            elif len(ass_haps) == 3:
                if counter_gene == 2 and counter_pseudo == 1:
                    two_cp_haps = [a for a in ass_haps.values() if "pms2cl" in a]
                elif counter_gene == 1 and counter_pseudo == 2:
                    two_cp_haps = [a for a in ass_haps.values() if "pms2hap" in a]

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
            self.init_het_sites,
            f"{self.genome_build}:{self.nchr}:{self.left_boundary}-{self.right_boundary}",
        )
