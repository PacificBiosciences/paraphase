# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>

import copy
from collections import namedtuple
import pysam
from ..phaser import Phaser


class F8Phaser(Phaser):
    new_fields = copy.deepcopy(Phaser.fields)
    new_fields.remove("gene_cn")
    new_fields.remove("alleles_final")
    new_fields.remove("hap_links")
    new_fields.insert(3, "sv_called")
    new_fields.insert(3, "flanking_summary")
    new_fields.insert(3, "exon1_to_exon22_depth")
    GeneCall = namedtuple(
        "GeneCall",
        new_fields,
        defaults=(None,) * len(new_fields),
    )

    def __init__(
        self, sample_id, outdir, genome_depth=None, genome_bam=None, sample_sex=None
    ):
        Phaser.__init__(self, sample_id, outdir, genome_depth, genome_bam, sample_sex)

    def set_parameter(self, config):
        super().set_parameter(config)
        self.depth_region = config["depth_region"]
        self.extract_region1, self.extract_region2, self.extract_region3 = config[
            "extract_regions"
        ].split()

    def get_read_positions(self, min_extension=5000):
        """Get mapped region of the part of reads not overlapping repeat"""
        dpos5 = {}
        dpos3 = {}
        genome_bamh = pysam.AlignmentFile(self.genome_bam, "rb")
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
        genome_bamh.close()
        return dpos5, dpos3

    def call(self):
        if self.check_coverage_before_analysis() is False:
            return self.GeneCall()

        genome_bamh = pysam.AlignmentFile(self.genome_bam, "rb")
        e1_e22_depth = self.get_regional_depth(genome_bamh, self.depth_region)[0]
        genome_bamh.close()

        dpos5, dpos3 = self.get_read_positions()
        self.get_homopolymer()
        self.get_candidate_pos()

        var_found = False
        for var in self.candidate_pos:
            pos = int(var.split("_")[0])
            if self.clip_3p_positions[0] < pos < self.clip_3p_positions[1]:
                var_found = True
                break
        if var_found is False:
            self.candidate_pos.add("155386300_A_C")

        var_found = False
        for var in self.candidate_pos:
            pos = int(var.split("_")[0])
            if pos > self.clip_3p_positions[1]:
                var_found = True
                break
        if var_found is False:
            self.candidate_pos.add("155386860_C_G")

        self.het_sites = sorted(list(self.candidate_pos))
        self.remove_noisy_sites()

        raw_read_haps = self.get_haplotypes_from_reads(
            check_clip=True, kept_sites=["155386300_A_C", "155386860_C_G"]
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

        total_cn = len(ass_haps)
        tmp = {}
        for i, hap in enumerate(ass_haps):
            tmp.setdefault(hap, f"hap{i+1}")
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
        for hap, links in flanking_sum.items():
            if links == "region1-region2":
                if self.mdepth is not None and self.sample_sex is not None:
                    if self.sample_sex == "female":
                        prob = self.depth_prob(int(e1_e22_depth), self.mdepth / 2)
                        if prob[0] > 0.75:
                            sv_hap.setdefault(hap, "deletion")
                    if self.sample_sex == "male":
                        if e1_e22_depth < 1:
                            sv_hap.setdefault(hap, "deletion")
            elif links == "region1-region3":
                sv_hap.setdefault(hap, "inversion")

        self.close_handle()

        return self.GeneCall(
            total_cn,
            ass_haps,
            [],
            e1_e22_depth,
            flanking_sum,
            sv_hap,
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
            self.region_avg_depth,
            self.sample_sex,
        )
