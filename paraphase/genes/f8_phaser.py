# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


from collections import namedtuple
import pysam
from ..phaser import Phaser


class F8Phaser(Phaser):
    GeneCall = namedtuple(
        "GeneCall",
        "total_cn final_haplotypes two_copy_haplotypes flanking_summary sv_called \
        highest_total_cn assembled_haplotypes sites_for_phasing \
        unique_supporting_reads het_sites_not_used_in_phasing homozygous_sites \
        haplotype_details variant_genotypes nonunique_supporting_reads \
        read_details genome_depth",
    )

    def __init__(self, sample_id, outdir, wgs_depth=None, genome_bam=None):
        Phaser.__init__(self, sample_id, outdir, wgs_depth, genome_bam)

    def set_parameter(self, config):
        super().set_parameter(config)
        self.extract_region1 = config["coordinates"]["hg38"]["extract_region1"]
        self.extract_region2, self.extract_region3 = config["coordinates"]["hg38"][
            "extract_region2"
        ].split()

    def get_read_positions(self, min_extension=5000):
        """Get mapped region of the part of reads not overlapping repeat"""
        dpos5 = {}
        dpos3 = {}
        genome_bamh = pysam.AlignmentFile(self.genome_bam, "rb")
        for i, extract_region in enumerate(
            [self.extract_region2, self.extract_region1, self.extract_region3]
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
            return None
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
        dvar = None
        if self.het_sites != []:
            haplotypes, dvar = self.output_variants_in_haplotypes(
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
        for hap, links in flanking_sum.items():
            if "region1" in links and links != "region1-region1":
                if links == "region1-region2":
                    sv_hap.setdefault(hap_name, "deletion")
                elif links == "region2-region1":
                    sv_hap.setdefault(hap_name, "duplication")
                elif "region3" in links:
                    sv_hap.setdefault(hap_name, "inversion")

        self.close_handle()

        return self.GeneCall(
            total_cn,
            ass_haps,
            [],
            flanking_sum,
            sv_hap,
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
