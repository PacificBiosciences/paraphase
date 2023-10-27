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
        self, sample_id, outdir, genome_depth=None, genome_bam=None, sample_sex=None
    ):
        Phaser.__init__(self, sample_id, outdir, genome_depth, genome_bam, sample_sex)
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
            return self.GeneCall()
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
        self.remove_noisy_sites()

        raw_read_haps = self.get_haplotypes_from_reads(
            check_clip=True,
            partial_deletion_reads=self.del1_reads_partial,
            kept_sites=[self.add_sites[1]],
            add_sites=[self.add_sites[0]],  # pivot site
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
                start_seq = hap[: pivot_index + 1]
                if start_seq.startswith("0") is False:
                    if len(start_seq) < 15:
                        unknown_counter += 1
                        hap_name = f"ikbkg_unknown_hap{unknown_counter}"
                    elif start_seq.count("2") <= 5:
                        gene_counter += 1
                        hap_name = f"ikbkg_hap{gene_counter}"
                    elif start_seq.count("2") >= 15:
                        pseudo_counter += 1
                        hap_name = f"ikbkg_pseudo_hap{pseudo_counter}"
                    else:
                        unknown_counter += 1
                        hap_name = f"ikbkg_unknown_hap{unknown_counter}"
                    tmp.setdefault(hap, hap_name)
                    if "3" in hap:
                        deletion_haplotypes.append(hap_name)
                else:
                    dup_counter += 1
                    tmp.setdefault(hap, f"ikbkg_dup_hap{dup_counter}")
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
                two_cp_haps = self.compare_depth(haplotypes, loose=True)
                if two_cp_haps == [] and read_counts is not None:
                    # check if one haplotype has more reads than others
                    haps = list(read_counts.keys())
                    counts = list(read_counts.values())
                    max_count = max(counts)
                    cp2_hap = haps[counts.index(max_count)]
                    others_max = sorted(counts, reverse=True)[1]
                    probs = self.depth_prob(max_count, others_max)
                    if probs[0] < 0.15 and others_max >= 10:
                        two_cp_haps.append(ass_haps[cp2_hap])
            for hap in two_cp_haps:
                total_cn += 1
                if "ikbkg_hap" in hap:
                    gene_counter += 1
            if gene_counter == 1 and pseudo_counter != 1:
                gene_counter = None
                total_cn = None

        (
            alleles,
            hap_links,
            _,
            _,
        ) = self.phase_alleles(
            uniquely_supporting_reads,
            nonuniquely_supporting_reads,
            raw_read_haps,
            ass_haps,
            reverse=self.is_reverse,
        )

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
            alleles,
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
        )
