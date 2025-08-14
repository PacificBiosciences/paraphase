# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>

import copy
from collections import namedtuple
import pysam
from ..phaser import Phaser
from paraphase.prepare_bam_and_vcf import pysam_handle


class F8Phaser(Phaser):
    new_fields = copy.deepcopy(Phaser.fields)
    new_fields.remove("gene_cn")
    new_fields.insert(3, "sv_called")
    new_fields.insert(3, "flanking_summary")
    new_fields.insert(3, "exon1_to_exon22_depth")
    GeneCall = namedtuple(
        "GeneCall",
        new_fields,
        defaults=(None,) * len(new_fields),
    )

    def __init__(
        self,
        sample_id,
        outdir,
        args=None,
        genome_depth=None,
        genome_bam=None,
        sample_sex=None,
    ):
        Phaser.__init__(
            self, sample_id, outdir, args, genome_depth, genome_bam, sample_sex
        )
        self.reference_fasta = None
        if args is not None:
            self.reference_fasta = args.reference

    def set_parameter(self, config):
        super().set_parameter(config)
        self.depth_region = config["depth_region"]
        self.extract_region1, self.extract_region2, self.extract_region3 = config[
            "extract_regions"
        ].split()

    def get_read_positions(self, genome_bamh, min_extension=5000):
        """Get mapped region of the part of reads not overlapping repeat"""
        dpos5 = {}
        dpos3 = {}

        for i, extract_region in enumerate(
            [self.extract_region1, self.extract_region2, self.extract_region3]
        ):
            pos1, pos2 = extract_region.split(":")[1].split("-")
            pos1 = int(pos1)
            pos2 = int(pos2)
            pos_name = f"region{i+1}"
            for pileupcolumn in genome_bamh.pileup(
                self.nchr,
                pos1 - 1,
                pos1,
                truncate=True,
            ):
                for read in pileupcolumn.pileups:
                    if read.alignment.reference_start < pos1 - min_extension:
                        read_name = read.alignment.query_name
                        if i in [0, 1]:
                            dpos5.setdefault(read_name, []).append(pos_name)
                        else:
                            dpos3.setdefault(read_name, []).append(pos_name)
            for pileupcolumn in genome_bamh.pileup(
                self.nchr,
                pos2 - 1,
                pos2,
                truncate=True,
            ):
                for read in pileupcolumn.pileups:
                    if read.alignment.reference_end > pos2 + min_extension:
                        read_name = read.alignment.query_name
                        if i in [0, 1]:
                            dpos3.setdefault(read_name, []).append(pos_name)
                        else:
                            dpos5.setdefault(read_name, []).append(pos_name)
        return dpos5, dpos3

    def call(self):
        if self.check_coverage_before_analysis() is False:
            return self.GeneCall(
                genome_depth=self.mdepth,
                region_depth=self.region_avg_depth._asdict(),
                sample_sex=self.sample_sex,
                phase_region=f"{self.genome_build}:{self.nchr}:{self.left_boundary}-{self.right_boundary}",
            )

        genome_bamh = pysam_handle(self.genome_bam, self.reference_fasta)
        e1_e22_depth = self.get_regional_depth(genome_bamh, self.depth_region)[0].median

        dpos5, dpos3 = self.get_read_positions(genome_bamh)
        genome_bamh.close()
        self.get_homopolymer()
        self.get_candidate_pos()

        var_found = False
        for var in self.candidate_pos:
            pos = int(var.split("_")[0])
            if self.clip_3p_positions[0] < pos < self.clip_3p_positions[1]:
                var_found = True
                break
        if var_found is False:
            self.candidate_pos.add(self.add_sites[0])

        var_found = False
        for var in self.candidate_pos:
            pos = int(var.split("_")[0])
            if pos > self.clip_3p_positions[1]:
                var_found = True
                break
        if var_found is False:
            self.candidate_pos.add(self.add_sites[1])

        self.het_sites = sorted(list(self.candidate_pos))
        self.remove_noisy_sites()
        self.init_het_sites = [a for a in self.het_sites]

        raw_read_haps = self.get_haplotypes_from_reads(
            check_clip=True,
            kept_sites=self.add_sites,
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

        total_cn = len(ass_haps)
        tmp = {}
        h1_count = 0
        h2_count = 0
        h3_count = 0
        unknown_count = 0
        inversion_count = 0
        deletion_count = 0
        invdup_count = 0
        for i, hap in enumerate(ass_haps):
            clip_5p = self.get_5pclip_from_hap(hap)
            clip_3p = self.get_3pclip_from_hap(hap)
            if clip_3p is None or clip_5p is None:
                unknown_count += 1
                hap_name = f"{self.gene}_unknownhap{unknown_count}"
            elif clip_5p == 0 and clip_3p == 0:
                h2_count += 1
                hap_name = f"{self.gene}_int22h2hap{h2_count}"
            elif (
                clip_5p == self.clip_5p_positions[0]
                and clip_3p == self.clip_3p_positions[0]
            ):
                h1_count += 1
                hap_name = f"{self.gene}_int22h1hap{h1_count}"
            elif clip_5p == 0 and clip_3p == self.clip_3p_positions[1]:
                h3_count += 1
                hap_name = f"{self.gene}_int22h3hap{h3_count}"
            elif (
                clip_5p == self.clip_5p_positions[0]
                and clip_3p == self.clip_3p_positions[1]
            ):
                inversion_count += 1
                hap_name = f"{self.gene}_int22invhap{inversion_count}"
            elif clip_5p == self.clip_5p_positions[0] and clip_3p == 0:
                deletion_count += 1
                hap_name = f"{self.gene}_int22delhap{deletion_count}"
            elif clip_5p == 0 and clip_3p == self.clip_3p_positions[0]:
                invdup_count += 1
                hap_name = f"{self.gene}_int22invorduphap{invdup_count}"
            else:
                unknown_count += 1
                hap_name = f"{self.gene}_unknownhap{unknown_count}"
            tmp.setdefault(hap, hap_name)
        ass_haps = tmp

        haplotypes = None
        if self.het_sites != []:
            haplotypes = self.output_variants_in_haplotypes(
                ass_haps,
                uniquely_supporting_reads,
                nonuniquely_supporting_reads,
            )

        # check flanking region, call sv
        haplotype_flanking_regions = {}
        # to-do: nonuniquely supporting reads
        for hap, reads in uniquely_supporting_reads.items():
            haplotype_flanking_regions.setdefault(hap, [[], []])
            for read in reads:
                if read in dpos5 and len(dpos5[read]) == 1:
                    haplotype_flanking_regions[hap][0].append(dpos5[read][0])
                if read in dpos3 and len(dpos3[read]) == 1:
                    haplotype_flanking_regions[hap][1].append(dpos3[read][0])
        flanking_sum = {}
        sv_hap = {}
        for hap, regions in haplotype_flanking_regions.items():
            hap_name = ass_haps[hap]
            p5region = list(set(regions[0]))
            p3region = list(set(regions[1]))
            flanking_sum.setdefault(
                hap_name, "-".join(["/".join(p5region), "/".join(p3region)])
            )

        # region2 and region3 are homologous at 5p for another 50kb,
        # so we cannot separate the upstream regions
        # we look for read evidence only at downstream region of region2 and region3
        # so we drop duplication as it involves upstream region of region2 and it's not pathogenic (?)
        for hap in ass_haps.values():
            if "del" in hap:
                if self.sample_sex is not None:
                    if self.sample_sex == "female" and self.mdepth is not None:
                        prob = self.depth_prob(int(e1_e22_depth), self.mdepth / 2)
                        if prob[0] > 0.75:
                            sv_hap.setdefault(hap, "deletion")
                    if self.sample_sex == "male":
                        if e1_e22_depth < 1:
                            sv_hap.setdefault(hap, "deletion")
            if "inv" in hap and "dup" not in hap:
                sv_hap.setdefault(hap, "inversion")
        """
        for hap, links in flanking_sum.items():
            if links == "region1-region2" and "int22h2" in hap:
                if self.sample_sex is not None:
                    if self.sample_sex == "female" and self.mdepth is not None:
                        prob = self.depth_prob(int(e1_e22_depth), self.mdepth / 2)
                        if prob[0] > 0.75:
                            sv_hap.setdefault(hap, "deletion")
                    if self.sample_sex == "male":
                        if e1_e22_depth < 1:
                            sv_hap.setdefault(hap, "deletion")
            elif links == "region1-region3" and "int22h3" in hap:
                sv_hap.setdefault(hap, "inversion")
        """
        if sv_hap == {} and self.sample_sex is not None:
            if self.sample_sex == "female" and total_cn < 6:
                total_cn = None
            if self.sample_sex == "male" and total_cn < 3:
                total_cn = None

        self.close_handle()

        return self.GeneCall(
            total_cn,
            ass_haps,
            [],
            e1_e22_depth,
            flanking_sum,
            sv_hap,
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
