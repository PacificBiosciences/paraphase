from collections import namedtuple
import re
import copy
from .phaser import Phaser


class Cyp21Phaser(Phaser):
    GeneCall = namedtuple(
        "GeneCall",
        "total_cn final_haplotypes two_copy_haplotypes starting_hap ending_hap deletion_hap \
        phasing_success alleles_simple annotated_alleles alleles hap_links \
        hap_variants highest_total_cn assembled_haplotypes gene1_read_number gene2_read_number het_sites \
        unique_supporting_reads het_sites_not_used_in_phasing homozygous_sites \
        haplotype_details variant_genotypes nonunique_supporting_reads \
        read_details genome_depth",
    )

    def __init__(self, sample_id, outdir, wgs_depth=None):
        Phaser.__init__(self, sample_id, outdir, wgs_depth)
        self.has_gene1 = False
        self.has_gene2 = False
        self.gene1_reads = set()
        self.gene2_reads = set()
        self.gene1_read_number = None
        self.gene2_read_number = None
        self.del1_reads = set()
        self.del1_reads_partial = set()
        self.del2_reads = set()
        self.del2_reads_partial = set()

    def set_parameter(self, config):
        super().set_parameter(config)
        self.variant_def = config["data"]["snp_file"]
        self.known_variants = {}
        with open(self.variant_def) as f:
            for line in f:
                split_line = line.split()
                self.known_variants.setdefault(
                    "_".join([split_line[1], split_line[2], split_line[3]]),
                    split_line[-1],
                )
        self.deletion1_size = config["coordinates"]["hg38"]["deletion1_size"]
        self.deletion2_size = config["coordinates"]["hg38"]["deletion2_size"]
        self.del2_3p_pos1 = config["coordinates"]["hg38"]["del2_3p_pos1"]
        self.del2_3p_pos2 = config["coordinates"]["hg38"]["del2_3p_pos2"]
        self.del2_5p_pos1 = config["coordinates"]["hg38"]["del2_5p_pos1"]
        self.del2_5p_pos2 = config["coordinates"]["hg38"]["del2_5p_pos2"]
        self.del1_3p_pos1 = config["coordinates"]["hg38"]["del1_3p_pos1"]
        self.del1_3p_pos2 = config["coordinates"]["hg38"]["del1_3p_pos2"]
        self.del1_5p_pos1 = config["coordinates"]["hg38"]["del1_5p_pos1"]
        self.del1_5p_pos2 = config["coordinates"]["hg38"]["del1_5p_pos2"]

    def allow_del_bases(self, pos):
        """
        During variant calling, allow "bases in deletions" for positions in
        two long deletions
        """
        if (
            self.del2_reads_partial != set()
            and self.del2_3p_pos1 <= pos <= self.del2_5p_pos2
        ):
            return True
        if (
            self.del1_reads_partial != set()
            and self.del1_3p_pos1 <= pos <= self.del1_5p_pos2
        ):
            return True
        return False

    @staticmethod
    def get_alleles(reads):
        """
        Phase haplotypes into alleles using read evidence
        """
        new_reads = {}
        for hap in reads:
            hap_reads = set()
            for read in reads[hap]:
                if "sup" not in read:
                    hap_reads.add(read)
                else:
                    hap_reads.add(read.split("_sup")[0])
            new_reads.setdefault(hap, hap_reads)
        links = {}
        checked = set()
        for hap1 in new_reads:
            r1 = new_reads[hap1]
            for hap2 in new_reads:
                hap_pair = "_".join(sorted([hap1, hap2]))
                if hap_pair not in checked and hap1 != hap2:
                    checked.add(hap_pair)
                    r2 = new_reads[hap2]
                    read_overlap = r1.intersection(r2)
                    # print(hap1, hap2, read_overlap)
                    if len(read_overlap) >= 2:
                        links.setdefault(hap1, []).append(hap2)
                        links.setdefault(hap2, []).append(hap1)
        links = dict(sorted(links.items(), key=lambda item: len(item[1]), reverse=True))
        # print(d)
        alleles = []
        if links != {}:
            alleles = [[list(links.keys())[0]] + list(links.values())[0]]
            for hap1 in links:
                for hap2 in links[hap1]:
                    hap1_in = sum([hap1 in a for a in alleles])
                    hap2_in = sum([hap2 in a for a in alleles])
                    if hap1_in == 0 and hap2_in == 0:
                        alleles.append([hap1, hap2])
                    elif hap1_in == 0:
                        for a in alleles:
                            if hap2 in a:
                                a_index = alleles.index(a)
                                if hap1 not in alleles[a_index]:
                                    alleles[a_index].append(hap1)
                                if hap2 not in alleles[a_index]:
                                    alleles[a_index].append(hap2)
                    else:
                        for a in alleles:
                            if hap1 in a:
                                a_index = alleles.index(a)
                                if hap2 not in alleles[a_index]:
                                    alleles[a_index].append(hap2)
                                if hap1 not in alleles[a_index]:
                                    alleles[a_index].append(hap1)
        return alleles, links

    def check_gene_presence(self):
        """
        check if either copy is present using a SNP site
        """
        bamh = self._bamh
        for pileupcolumn in bamh.pileup(
            self.nchr, self.pivot_site - 1, self.pivot_site, truncate=True
        ):
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    read_seq = pileupread.alignment.query_sequence
                    read_pos = pileupread.query_position
                    base1 = read_seq[read_pos].upper()
                    if base1 == "T":
                        self.gene1_reads.add(pileupread.alignment.query_name)
                    elif base1 == "A":
                        self.gene2_reads.add(pileupread.alignment.query_name)

        self.gene1_read_number = len(self.gene1_reads)
        self.gene2_read_number = len(self.gene2_reads)
        if self.gene1_read_number >= 2:
            self.has_gene1 = True
        if self.gene2_read_number >= 2:
            self.has_gene2 = True

    def output_variants_in_haplotypes(self, haps, reads, nonunique, two_cp_haps=[]):
        """
        Summarize all variants in each haplotype.
        Output all variants and their genotypes.
        Haplotypes are different length, so a range (boundary) is reported
        """
        het_sites = self.het_sites
        haplotype_variants = {}
        haplotype_info = {}
        dvar = {}
        var_no_phasing = copy.deepcopy(self.het_no_phasing)
        for hap, hap_name in haps.items():
            haplotype_variants.setdefault(hap_name, [])
        # het sites not used in phasing
        if reads != {}:
            for var in var_no_phasing:
                genotypes = []
                var_reads = self.check_variants_in_haplotypes(var)
                haps_with_variant = []
                for hap, hap_name in haps.items():
                    hap_reads = reads[hap]
                    hap_reads_nonunique = [a for a in nonunique if hap in nonunique[a]]
                    genotype = self.get_genotype_in_hap(
                        var_reads, hap_reads, hap_reads_nonunique
                    )
                    genotypes.append(genotype)
                    if genotype == "1":
                        haps_with_variant.append(hap_name)
                if haps_with_variant == []:
                    self.het_no_phasing.remove(var)
                else:
                    for hap_name in haps_with_variant:
                        haplotype_variants[hap_name].append(var)
                    dvar.setdefault(var, genotypes)
        # het sites and homo sites
        for hap, hap_name in haps.items():
            for i in range(len(hap)):
                if hap[i] == "2":
                    haplotype_variants[hap_name].append(het_sites[i])
                elif (
                    hap[i] == "3"
                    and "32043718_del120" not in haplotype_variants[hap_name]
                ):
                    haplotype_variants[hap_name].append("32043718_del120")
                elif (
                    hap[i] == "4"
                    and "32017431_del6367" not in haplotype_variants[hap_name]
                ):
                    haplotype_variants[hap_name].append("32017431_del6367")
            if "32017431_del6367" in haplotype_variants[hap_name]:
                pos1 = self.del1_3p_pos1
                pos2 = self.del1_5p_pos2
                for var in self.homo_sites:
                    pos = int(var.split("_")[0])
                    if pos < pos1 or pos > pos2:
                        haplotype_variants[hap_name].append(var)
            elif "32043718_del120" in haplotype_variants[hap_name]:
                pos1 = self.del2_3p_pos1
                pos2 = self.del2_5p_pos2
                for var in self.homo_sites:
                    pos = int(var.split("_")[0])
                    if pos < pos1 or pos > pos2:
                        haplotype_variants[hap_name].append(var)
            else:
                haplotype_variants[hap_name] += self.homo_sites

            var_nstart, var_nend = self.get_hap_variant_ranges(hap)
            var_tmp = haplotype_variants[hap_name]
            var_tmp1 = [
                a for a in var_tmp if var_nstart <= int(a.split("_")[0]) <= var_nend
            ]
            var_tmp1 = list(set(var_tmp1))
            var_tmp2 = sorted(var_tmp1, key=lambda x: int(x.split("_")[0]))
            haplotype_info.setdefault(
                hap_name, {"variants": var_tmp2, "boundary": [var_nstart, var_nend]}
            )

        # summary per variant
        all_haps = haps
        nhap = len(all_haps)
        for var in self.homo_sites:
            dvar.setdefault(var, ["1"] * nhap)
        for i, var in enumerate(het_sites):
            dvar.setdefault(var, [])
            for hap, hap_name in haps.items():
                base_call = "."
                if hap[i] == "2":
                    base_call = "1"
                elif hap[i] == "1":
                    base_call = "0"
                dvar[var].append(base_call)
                if hap_name in two_cp_haps:
                    dvar[var].append(base_call)

        return haplotype_info, {
            var: "|".join(dvar[var]) for var in dict(sorted(dvar.items()))
        }

    def annotate_var(self, allele_var):
        """annotate an allele with variants"""
        annotated_allele = None
        if len(allele_var) == 2:
            if [] in allele_var:
                annotated_allele = "WT"
            else:
                tmp = sorted(allele_var, key=lambda x: len(x))
                annotated_allele = ",".join(tmp[0])
        elif len(allele_var) == 1:
            if allele_var == [[]]:
                annotated_allele = "pseudogene_deletion"
            else:
                annotated_allele = "deletion_" + ",".join(allele_var[0])
        elif len(allele_var) == 3:
            tmp = sorted(allele_var, key=lambda x: len(x))
            if tmp[0] == [] and tmp[1] == []:
                annotated_allele = "gene_duplication"
            elif tmp[0] == []:
                if abs(len(tmp[1]) - len(tmp[2])) <= 1:
                    annotated_allele = "pseudogene_duplication"
                else:
                    annotated_allele = "duplicaton_plus_" + ",".join(tmp[1])
        return annotated_allele

    def call(self):
        self.get_homopolymer()
        self.del2_reads, self.del2_reads_partial = self.get_long_del_reads(
            self.del2_3p_pos1,
            self.del2_3p_pos2,
            self.del2_5p_pos1,
            self.del2_5p_pos2,
            self.deletion2_size,
        )
        self.del1_reads, self.del1_reads_partial = self.get_long_del_reads(
            self.del1_3p_pos1,
            self.del1_3p_pos2,
            self.del1_5p_pos1,
            self.del1_5p_pos2,
            self.deletion1_size,
        )

        self.check_gene_presence()
        # if self.has_gene1 is False and self.has_gene2 is False:
        #    raise Exception("Not enough reads to analyze this region.")

        # scan for polymorphic sites
        regions_to_check = []
        if self.del2_reads_partial != set():
            regions_to_check += [
                [self.del2_3p_pos1, self.del2_3p_pos2],
                [self.del2_5p_pos1, self.del2_5p_pos2],
            ]
        if self.del1_reads_partial != set():
            regions_to_check += [
                [self.del1_3p_pos1, self.del1_3p_pos2],
                [self.del1_5p_pos1, self.del1_5p_pos2],
            ]
        self.get_candidate_pos(regions_to_check=regions_to_check)
        # always add splice site
        if self.candidate_pos != set():
            self.candidate_pos.add("32039816_T_A")

        # add last snp outside of repeat
        var_found = False
        for var in self.candidate_pos:
            pos = int(var.split("_")[0])
            if pos > self.clip_3p_positions[0]:
                var_found = True
                break
        if var_found is False and self.candidate_pos != set():
            self.candidate_pos.add("32046300_G_A")
        var_found = False
        for var in self.candidate_pos:
            pos = int(var.split("_")[0])
            if pos < self.clip_5p_positions[0]:
                var_found = True
                break
        if var_found is False and self.candidate_pos != set():
            self.candidate_pos.add("32013265_A_T")

        het_sites = sorted(list(self.candidate_pos))
        problematic_sites = []
        for site in het_sites:
            for region in self.noisy_region:
                if region[0] <= int(site.split("_")[0]) <= region[1]:
                    problematic_sites.append(site)
        for site in problematic_sites:
            het_sites.remove(site)

        raw_read_haps = self.get_haplotypes_from_reads(
            het_sites, check_clip=True, partial_deletion_reads=self.del1_reads_partial
        )
        if self.del2_reads_partial != set():
            raw_read_haps, het_sites = self.update_reads_for_deletions(
                raw_read_haps,
                het_sites,
                self.del2_3p_pos1,
                self.del2_5p_pos2,
                self.del2_reads_partial,
                "3",
                "32043718_del120",
            )
        if self.del1_reads_partial != set():
            raw_read_haps, het_sites = self.update_reads_for_deletions(
                raw_read_haps,
                het_sites,
                self.del1_3p_pos1,
                self.del1_5p_pos2,
                self.del1_reads_partial,
                "4",
                "32017431_del6367",
            )
        self.het_sites = het_sites

        # assemble haplotypes
        (
            ass_haps,
            original_haps,
            hcn,
            uniquely_supporting_reads,
            nonuniquely_supporting_reads,
            raw_read_haps,
            read_counts,
        ) = self.phase_haps(raw_read_haps)

        tmp1 = {}
        tmp2 = {}
        for i, hap in enumerate(ass_haps):
            hap_name = f"hap{i+1}"
            tmp1.setdefault(hap, hap_name)
            tmp2.setdefault(hap_name, hap)
        final_haps = tmp1
        # get haps that extend into tnxb
        ending_copies = [
            final_haps[a]
            for a in ass_haps
            if a[0] not in ["0", "x"] and a[-1] not in ["0", "x"]
        ]
        starting_copies = [
            final_haps[a] for a in ass_haps if a[0] == "0" and a[-1] == "0"
        ]
        single_copies = [
            final_haps[a] for a in ass_haps if a[0] == "0" and a[-1] not in ["0", "x"]
        ]

        haplotypes = None
        dvar = None
        if self.het_sites != []:
            haplotypes, dvar = self.output_variants_in_haplotypes(
                final_haps,
                uniquely_supporting_reads,
                nonuniquely_supporting_reads,
            )

        successful_phasing = False
        # phase haplotypes into alleles
        alleles, links = self.get_alleles(uniquely_supporting_reads)
        # switch to hap name
        new_alleles = []
        for pair in alleles:
            new_pair = []
            for hap1 in pair:
                new_pair.append(final_haps[hap1])
            new_alleles.append(new_pair)
        new_links = {}
        for hap in links:
            hap_links = [final_haps[a] for a in links[hap]]
            new_links.setdefault(final_haps[hap], []).append(hap_links)
        links = new_links
        # print(alleles)
        # print(links)
        two_cp_haplotypes = []
        # the deletion haplotype will be reported as an allele
        if len(single_copies) == 1 and len(ass_haps) < 5:
            if (
                len(new_alleles) == 1
                and len(new_alleles[0]) == len(ass_haps) - 1
                and single_copies not in new_alleles
            ):
                new_alleles.append(single_copies)
            elif (
                len(new_alleles) == 1
                and len(new_alleles[0]) < len(ass_haps) - 1
                and single_copies not in new_alleles
                and len(starting_copies) == 1
                and len(ending_copies) == 1
            ):
                new_alleles = []
            if new_alleles == []:
                new_alleles.append(single_copies)
                remaining_hap = [
                    a for a in final_haps.values() if a not in single_copies
                ]
                if len(remaining_hap) == len(ass_haps) - 1:
                    new_alleles.append(remaining_hap)
        elif (
            len(single_copies) == 2
            and len(alleles) == 0
            and len(ass_haps) == 2
            and len(starting_copies) == 0
            and len(ending_copies) == 0
        ):
            new_alleles = [[single_copies[0]], [single_copies[1]]]
            successful_phasing = True
        elif single_copies == []:
            # 2 cp, one allele, one haplotype extends to tnxb -> each haplotype has cn 2
            if (
                len(ass_haps) == 2
                and len(ending_copies) == 1
                and len(starting_copies) == 1
                and len(alleles) == 1
            ):
                two_cp_haplotypes = [final_haps[a] for a in ass_haps]
                new_alleles.append(new_alleles[0])
                successful_phasing = True

            # depth-based adjustment when found 3 haplotypes or <2 ending haplotypes
            if haplotypes is not None:
                two_cp_hap_candidate = self.compare_depth(haplotypes, loose=True)
                if len(ending_copies) == 1 and len(starting_copies) == 2:
                    if two_cp_hap_candidate == ending_copies:
                        two_cp_haplotypes = two_cp_hap_candidate
                        if len(ass_haps) == 3:
                            new_alleles = [
                                [starting_copies[0], ending_copies[0]],
                                [starting_copies[1], ending_copies[0]],
                            ]
                            successful_phasing = True
                elif len(starting_copies) == 1 and len(ending_copies) == 2:
                    if two_cp_hap_candidate == starting_copies:
                        two_cp_haplotypes = two_cp_hap_candidate
                        if len(ass_haps) == 3:
                            new_alleles = [
                                [starting_copies[0], ending_copies[0]],
                                [starting_copies[0], ending_copies[1]],
                            ]
                            successful_phasing = True

            # add the missing link in cn=4
            if (
                len(ass_haps) in [3, 4]
                and two_cp_haplotypes == []
                and len(new_alleles) == 1
                and len(new_alleles[0]) == 2
            ):
                remaining_hap = [
                    a for a in final_haps.values() if a not in new_alleles[0]
                ]
                if len(remaining_hap) == len(ass_haps) - 2:
                    new_alleles.append(remaining_hap)
            # add the missing link in cn=5
            if len(ass_haps) == 5 and two_cp_haplotypes == []:
                if (
                    len(new_alleles) == 1
                    and len(new_alleles[0]) == 2
                    and (
                        (
                            new_alleles[0][0] in starting_copies
                            and new_alleles[0][1] in ending_copies
                        )
                        or (
                            new_alleles[0][1] in starting_copies
                            and new_alleles[0][0] in ending_copies
                        )
                    )
                ):
                    remaining_hap = [
                        a for a in final_haps.values() if a not in new_alleles[0]
                    ]
                    if len(remaining_hap) == 3:
                        new_alleles.append(remaining_hap)
                if len(new_alleles) == 2:
                    allele1 = (
                        new_alleles[0][0] in starting_copies
                        and new_alleles[0][1] in ending_copies
                    ) or (
                        new_alleles[0][1] in starting_copies
                        and new_alleles[0][0] in ending_copies
                    )
                    allele2 = (
                        new_alleles[1][0] in starting_copies
                        and new_alleles[1][1] in ending_copies
                    ) or (
                        new_alleles[1][1] in starting_copies
                        and new_alleles[1][0] in ending_copies
                    )
                    if allele1 is True and allele2 is False:
                        remaining_hap = [
                            a for a in final_haps.values() if a not in new_alleles[0]
                        ]
                        if len(remaining_hap) == 3:
                            new_alleles = [new_alleles[0], remaining_hap]
                    elif allele1 is False and allele2 is True:
                        remaining_hap = [
                            a for a in final_haps.values() if a not in new_alleles[1]
                        ]
                        if len(remaining_hap) == 3:
                            new_alleles = [new_alleles[1], remaining_hap]

        if len(new_alleles) == 2:
            if sorted(new_alleles[0] + new_alleles[1]) == sorted(
                list(final_haps.values())
            ):
                successful_phasing = True

        # check wrong phasing
        wrong_allele = False
        for allele in new_alleles:
            hp_set = set(allele)
            for hp in hp_set:
                if (
                    hp in ending_copies
                    and allele.count(hp) > 1
                    and hp not in two_cp_haplotypes
                ):
                    wrong_allele = True
                    break
        if wrong_allele:
            new_alleles = []
            alleles = []

        # annotate haplotypes by checking the diff sites
        # output variants carried by each haplotype
        hap_variants = {}
        if haplotypes is not None:
            for hap, hap_info in haplotypes.items():
                hap_variants.setdefault(hap, [])
                for var in hap_info["variants"]:
                    if var in self.known_variants:
                        hap_variants[hap].append(self.known_variants[var])
                    elif "del120" in var:
                        hap_variants[hap].append(var)

        total_cn = len(ass_haps) + len(two_cp_haplotypes)
        if ass_haps == [] and self.het_sites == []:
            # feed all reads to call variants
            if True in [self.has_gene1, self.has_gene2] and False in [
                self.has_gene1,
                self.has_gene2,
            ]:
                total_cn = 2
        if total_cn == 0:
            total_cn = None

        annotated_alleles = []
        if successful_phasing:
            allele1 = new_alleles[0]
            allele2 = new_alleles[1]
            allele1_var = [hap_variants[a] for a in allele1]
            allele2_var = [hap_variants[a] for a in allele2]
            for allele_var in [allele1_var, allele2_var]:
                annotated_allele = self.annotate_var(allele_var)
                annotated_alleles.append(annotated_allele)
        else:
            if (
                len(ending_copies) == 2
                and len(ass_haps) == 4
                and two_cp_haplotypes == []
            ):
                for allele_var in [hap_variants[a] for a in ending_copies]:
                    annotated_allele = None
                    if allele_var == [[]]:
                        annotated_allele = "WT"
                    else:
                        annotated_allele = ",".join(allele_var)
                annotated_alleles.append(annotated_allele)
        self.close_handle()

        return self.GeneCall(
            total_cn,
            final_haps,
            two_cp_haplotypes,
            starting_copies,
            ending_copies,
            single_copies,
            successful_phasing,
            new_alleles,
            annotated_alleles,
            alleles,
            new_links,
            hap_variants,
            hcn,
            original_haps,
            self.gene1_read_number,
            self.gene2_read_number,
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
