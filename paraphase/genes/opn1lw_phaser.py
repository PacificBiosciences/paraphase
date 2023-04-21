# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


from collections import namedtuple
from ..phaser import Phaser


class Opn1lwPhaser(Phaser):
    GeneCall = namedtuple(
        "GeneCall",
        "total_cn opn1lw_cn opn1mw_cn final_haplotypes two_copy_haplotypes first_copies \
        last_copies middle_copies annotated_haplotypes phasing_success annotated_alleles alleles_final \
        alleles_full hap_links directional_links links_loose alleles_raw highest_total_cn assembled_haplotypes \
        sites_for_phasing unique_supporting_reads het_sites_not_used_in_phasing \
        homozygous_sites haplotype_details variant_genotypes nonunique_supporting_reads \
        read_details genome_depth sample_sex",
    )
    exon3_variants = {
        154152987: (["154152987_A_C"], "L", "M"),
        154153041: (["154153041_G_A", "154153043_G_T"], "I", "V"),
        154153051: (["154153051_C_T"], "V", "A"),
        154153062: (["154153062_A_G"], "V", "I"),
        154153068: (["154153068_G_T"], "S", "A"),
    }
    pathogenic_haps = ["LIAVA", "LVAVA", "LIAVS", "MIAVA", "MVVVA", "MVAVA", "LIAIA"]

    def __init__(self, sample_id, outdir, wgs_depth=None, sample_sex=None):
        Phaser.__init__(
            self, sample_id, outdir, wgs_depth=wgs_depth, sample_sex=sample_sex
        )

    def set_parameter(self, config):
        super().set_parameter(config)

    def call_exon3(self, vars):
        annotated_vars = [None] * 5
        i = 0
        for pos, pos_item in self.exon3_variants.items():
            query_vars, alt, ref = pos_item
            if len(set(vars).intersection(set(query_vars))) == len(query_vars):
                annotated_vars[i] = alt
            else:
                annotated_vars[i] = ref
            i += 1
        if None not in annotated_vars:
            return "".join(annotated_vars)
        return None

    @staticmethod
    def phase_single_allele(
        first_copy,
        all_haps_on_this_allele,
        last_copies,
        middle_copies,
        new_links,
        alleles,
        directional_links,
        directional_links_loose,
    ):
        hap = first_copy
        remaining_haps = [a for a in all_haps_on_this_allele if a != hap]
        middle_copies = [a for a in middle_copies if a in all_haps_on_this_allele]
        # two copies
        if len(remaining_haps) == 1:
            return [hap, remaining_haps[0]]

        # three copies, last one known
        remaining_noending_haps = [a for a in remaining_haps if a not in last_copies]
        if len(remaining_noending_haps) == 1:
            return [hap, remaining_noending_haps[0]]
        # three total copies. we don't know the last copy but we know the middle copy
        if len(remaining_haps) == 2 and len(middle_copies) == 1:
            return [hap, middle_copies[0]]

        # four copies and we know the one before the last copy
        if len(last_copies) == 1:
            last_copy_before = new_links.get(last_copies[0])
            if last_copy_before is not None and len(last_copy_before) == 1:
                remaining_noending_not_before_ending_haps = [
                    a for a in remaining_noending_haps if a not in last_copy_before
                ]
                if len(remaining_noending_not_before_ending_haps) == 1:
                    return [
                        hap,
                        remaining_noending_not_before_ending_haps[0],
                    ]
        # four total copies. we know the directional link between the middle two copies
        if len(remaining_noending_haps) == 2:
            hap1, hap2 = remaining_noending_haps
            if (
                f"{hap1}-{hap2}" in directional_links
                and f"{hap2}-{hap1}" not in directional_links
            ):
                return [hap, hap1]
            if (
                f"{hap2}-{hap1}" in directional_links
                and f"{hap1}-{hap2}" not in directional_links
            ):
                return [hap, hap2]

        # if other haps are phased on this allele
        if len(alleles) == 1 and last_copies[0] in alleles[0]:
            remaining_haps_not_phased = [
                a for a in remaining_noending_haps if a not in alleles[0]
            ]
            if len(remaining_haps_not_phased) == 1:
                return [
                    hap,
                    remaining_haps_not_phased[0],
                ]
            if len(remaining_haps_not_phased) == 2:
                hap1, hap2 = remaining_haps_not_phased
                if (
                    f"{hap1}-{hap2}" in directional_links
                    and f"{hap2}-{hap1}" not in directional_links
                ):
                    return [hap, hap1]
                if (
                    f"{hap2}-{hap1}" in directional_links
                    and f"{hap1}-{hap2}" not in directional_links
                ):
                    return [hap, hap2]

        # loose links of the first copy
        hap_next_loose = directional_links_loose.get(hap)
        if hap_next_loose is not None:
            hap_next_loose_noending = [
                a
                for a in hap_next_loose
                if a not in last_copies and hap_next_loose[a] >= 3
            ]
            if (
                len(hap_next_loose_noending) == 1
                and hap_next_loose_noending[0] in all_haps_on_this_allele
            ):
                return [hap, hap_next_loose_noending[0]]

        return None

    def call(self):
        if self.check_coverage_before_analysis() is False:
            return None
        self.get_homopolymer()
        self.get_candidate_pos()
        self.het_sites = sorted(list(self.candidate_pos))
        self.remove_noisy_sites()

        raw_read_haps = self.get_haplotypes_from_reads(
            check_clip=True, add_sites=["154143410_G_A"]
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
        for i, hap in enumerate(ass_haps):
            tmp.setdefault(hap, f"hap{i+1}")
        ass_haps = tmp

        haplotypes = None
        dvar = None
        two_cp_haps = []
        alleles = []
        new_alleles = []
        new_links = {}
        directional_links = {}
        directional_links_loose = {}
        alleles_1st_2nd = []
        first_copies = []
        last_copies = []
        middle_copies = []
        counter_lw = 0
        counter_mw = 0
        counter_unknown = 0
        annotated_alleles = []
        annotated_haps = {}
        phasing_success = False
        pivot_vars = set(["154156379_A_T", "154156402_A_G"])
        last_copy_vars = set(["154182157_T_G", "154182166_A_C", "154182167_G_A"])

        if self.het_sites != []:
            haplotypes, dvar = self.output_variants_in_haplotypes(
                ass_haps,
                uniquely_supporting_reads,
                nonuniquely_supporting_reads,
            )

            hap_rename = {}
            for hap, hap_name in ass_haps.items():
                var = set(haplotypes[hap_name]["variants"])
                var_intersect = var.intersection(pivot_vars)
                if var_intersect == set():
                    counter_lw += 1
                    renamed_hap = f"opn1lw_hap{counter_lw}"
                elif len(var_intersect) == 2:
                    counter_mw += 1
                    renamed_hap = f"opn1mw_hap{counter_mw}"
                else:
                    counter_unknown += 1
                    renamed_hap = f"opnunknown_hap{counter_unknown}"
                hap_rename.setdefault(hap_name, renamed_hap)
                if hap[0] not in ["x", "0"]:
                    first_copies.append(renamed_hap)
                if var.intersection(last_copy_vars) != set():
                    last_copies.append(renamed_hap)
                elif renamed_hap not in first_copies and "x" not in hap:
                    middle_copies.append(renamed_hap)
                # annotate exon3
                hap_annotated = self.call_exon3(var)
                gene_annotated = renamed_hap.split("_")[0]
                if hap_annotated in self.pathogenic_haps:
                    gene_annotated += "_" + hap_annotated
                annotated_haps.setdefault(renamed_hap, gene_annotated)

            tmp = {}
            for hap, hap_name in ass_haps.items():
                tmp.setdefault(hap, hap_rename[hap_name])
            ass_haps = tmp

            tmp = {}
            for hap_name, hap_info in haplotypes.items():
                tmp.setdefault(hap_rename[hap_name], hap_info)
            haplotypes = tmp

            # find two copy haplotypes
            if len(ass_haps) > 2 and self.sample_sex == "female":
                two_cp_haps = self.compare_depth(haplotypes)
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

            total_cn = len(ass_haps) + len(two_cp_haps)

            (
                new_alleles,
                new_links,
                directional_links,
                directional_links_loose,
            ) = self.get_alleles_directional(
                uniquely_supporting_reads,
                nonuniquely_supporting_reads,
                raw_read_haps,
                ass_haps,
            )
            # from pprint import pprint
            # pprint(directional_links)
            # pprint(directional_links_loose)
            # for allele in alleles:
            #    new_allele = []
            #    for hap in allele:
            #        new_allele.append(ass_haps[hap])
            #    new_alleles.append(new_allele)

            # for hap in links:
            #    hap_links = [ass_haps[a] for a in links[hap]]
            #    new_links.setdefault(ass_haps[hap], hap_links)  # change in rccx?

            # incorrect phasing suggests haplotypes with cn > 1
            if (
                self.sample_sex == "female"
                and len(new_alleles) == 1
                and sorted(new_alleles[0]) == sorted(ass_haps.values())
            ):
                new_alleles = []
                total_cn = None

            # focus on first and second copies
            if self.sample_sex is not None:
                if (self.sample_sex == "male" and len(first_copies) == 1) or (
                    self.sample_sex == "female" and len(first_copies) == 2
                ):
                    for hap in first_copies:
                        hap_next = new_links.get(hap)
                        if hap_next is not None and len(hap_next) == 1:
                            alleles_1st_2nd.append([hap, hap_next[0]])
                        elif hap_next is None and hap in last_copies:
                            alleles_1st_2nd.append([hap])
                        else:
                            hap_next_loose = directional_links_loose.get(hap)
                            if hap_next_loose is not None and len(hap_next_loose) == 1:
                                next_hap_name, read_count = list(
                                    hap_next_loose.items()
                                )[0]
                                if read_count >= 3:
                                    alleles_1st_2nd.append([hap, next_hap_name])

                    # easier phasing for males as there is one allele
                    if self.sample_sex == "male" and alleles_1st_2nd == []:
                        allele_found = self.phase_single_allele(
                            first_copies[0],
                            list(ass_haps.values()),
                            last_copies,
                            middle_copies,
                            new_links,
                            new_alleles,
                            directional_links,
                            directional_links_loose,
                        )
                        if allele_found is not None:
                            alleles_1st_2nd.append(allele_found)
                    if self.sample_sex == "female" and len(alleles_1st_2nd) == 1:
                        complete_allele = []
                        incomplete_allele = []
                        for each_allele in new_alleles:
                            if (
                                set(each_allele).intersection(set(first_copies))
                                != set()
                                and set(each_allele).intersection(set(last_copies))
                                != set()
                            ):
                                complete_allele.append(each_allele)
                            else:
                                incomplete_allele.append(each_allele)
                        if len(complete_allele) == 1:
                            other_allele_copies = [
                                a
                                for a in ass_haps.values()
                                if a not in complete_allele[0]
                            ]
                            other_first_copy = [
                                a for a in first_copies if a not in complete_allele[0]
                            ]
                            other_last_copy = [
                                a for a in last_copies if a not in complete_allele[0]
                            ]
                            if len(other_first_copy) == 1 and len(other_last_copy) == 1:
                                allele_found = self.phase_single_allele(
                                    other_first_copy[0],
                                    other_allele_copies,
                                    other_last_copy,
                                    middle_copies,
                                    new_links,
                                    incomplete_allele,
                                    directional_links,
                                    directional_links_loose,
                                )
                                if allele_found is not None:
                                    alleles_1st_2nd.append(allele_found)

                    if (self.sample_sex == "male" and len(alleles_1st_2nd) == 1) or (
                        self.sample_sex == "female" and len(alleles_1st_2nd) == 2
                    ):
                        phasing_success = True
                    else:
                        for hap in first_copies:
                            if True not in [hap in a for a in alleles_1st_2nd]:
                                alleles_1st_2nd.append([hap, None])

                for each_allele in alleles_1st_2nd:
                    each_allele_annotated = []
                    for each_hap in each_allele:
                        if each_hap is not None:
                            hap_vars = haplotypes[each_hap]["variants"]
                            hap_annotated = self.call_exon3(hap_vars)
                            gene_annotated = each_hap.split("_")[0]
                            if hap_annotated in self.pathogenic_haps:
                                gene_annotated += "_" + hap_annotated
                        else:
                            gene_annotated = None
                        each_allele_annotated.append(gene_annotated)
                    annotated_alleles.append(each_allele_annotated)

        # homozygous cases
        else:
            baseline_cn = None
            phasing_success = True
            if self.sample_sex == "female":
                baseline_cn = 2
            elif self.sample_sex == "male":
                baseline_cn = 1
            hap_annotated = self.call_exon3(self.homo_sites)
            if set(self.homo_sites).intersection(pivot_vars) != set():
                total_cn = baseline_cn
                counter_mw = baseline_cn
                gene_annotated = "opn1mw"
            else:
                total_cn = baseline_cn
                counter_lw = baseline_cn
                gene_annotated = "opn1lw"
            if hap_annotated in self.pathogenic_haps:
                gene_annotated += "_" + hap_annotated
            annotated_alleles.append([gene_annotated])
            if self.sample_sex == "female":
                annotated_alleles.append([gene_annotated])

        self.close_handle()

        return self.GeneCall(
            total_cn,
            counter_lw,
            counter_mw,
            ass_haps,
            two_cp_haps,
            first_copies,
            last_copies,
            middle_copies,
            annotated_haps,
            phasing_success,
            annotated_alleles,
            alleles_1st_2nd,
            new_alleles,
            new_links,
            directional_links,
            directional_links_loose,
            alleles,
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
            self.sample_sex,
        )
