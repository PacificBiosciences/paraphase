from collections import namedtuple
import re
from .phaser import Phaser


class Cyp21Phaser(Phaser):
    GeneCall = namedtuple(
        "GeneCall",
        "total_cn final_haplotypes two_copy_haplotypes ending_hap gene1_read_number gene2_read_number alleles \
        highest_total_cn assembled_haplotypes het_sites \
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

    """
    def get_long_del_reads(
        self,
        p3_pos1,
        p3_pos2,
        p5_pos1,
        p5_pos2,
        del_size,
        min_clip_len=300,
        min_extend=1000,
    ):
        bamh = self._bamh
        p5_reads = set()
        p3_reads = set()
        del_reads = set()
        # SMN1, 3 prime clip
        pos1 = p3_pos1
        pos2 = p3_pos2
        reference_start_cutoff = pos1 - min_extend
        for read in bamh.fetch(self.nchr, pos1, pos2):
            read_name = read.query_name
            if read.is_supplementary:
                read_name = (
                    read_name + f"_sup_{read.reference_start}_{read.reference_length}"
                )
            find_clip_3p = re.findall(self.clip_3p, read.cigarstring)
            if find_clip_3p != [] and pos1 < read.reference_end < pos2:
                if (
                    int(find_clip_3p[0][:-1]) >= min_clip_len
                    and read.reference_start < reference_start_cutoff
                ):
                    p3_reads.add(read_name)
            if self.check_del(read, del_size):
                del_reads.add(read_name)
        # SMN1, 5 prime clip
        pos1 = p5_pos1
        pos2 = p5_pos2
        reference_end_cutoff = pos2 + min_extend
        for read in bamh.fetch(self.nchr, pos1, pos2):
            read_name = read.query_name
            if read.is_supplementary:
                read_name = (
                    read_name + f"_sup_{read.reference_start}_{read.reference_length}"
                )
            find_clip_5p = re.findall(self.clip_5p, read.cigarstring)
            if find_clip_5p != [] and pos1 < read.reference_start < pos2:
                if (
                    int(find_clip_5p[0][:-1]) >= min_clip_len
                    and read.reference_end > reference_end_cutoff
                ):
                    p5_reads.add(read_name)
            if self.check_del(read, del_size):
                del_reads.add(read_name)
        if del_reads != set() or (p3_reads != set() and p5_reads != set()):
            return (
                del_reads.union(p3_reads.intersection(p5_reads)),
                del_reads.union(p3_reads).union(p5_reads),
            )
        return set(), set()

    def get_haplotypes_from_reads(self, het_sites, exclude_reads=[], min_mapq=5):
        read_haps = {}
        nvar = len(het_sites)
        for dsnp_index, allele_site in enumerate(het_sites):
            snp_position_gene1, allele1, allele2, *at = allele_site.split("_")
            snp_position = int(snp_position_gene1)
            for pileupcolumn in self._bamh.pileup(
                self.nchr,
                snp_position - 1,
                snp_position,
                truncate=True,
                min_base_quality=29,  # lowered base quality cutoff
                # min_base_quality=30,
            ):
                for read in pileupcolumn.pileups:
                    read_names = [read.alignment.query_name]
                    if read.alignment.is_supplementary:
                        sup_name = (
                            read.alignment.query_name
                            + f"_sup_{read.alignment.reference_start}_{read.alignment.reference_length}"
                        )
                        read_names = [sup_name]
                        if (
                            sup_name in self.del1_reads_partial
                            and read.alignment.query_name in self.del1_reads_partial
                        ):
                            read_names.append(read.alignment.query_name)
                    for read_name in read_names:
                        if (
                            not read.is_del
                            and not read.is_refskip
                            and read.alignment.is_secondary == 0
                            and read.alignment.is_duplicate == 0
                            and read.alignment.mapping_quality >= min_mapq
                            and read_name not in exclude_reads
                        ):
                            read_seq = read.alignment.query_sequence
                            start_pos = read.query_position
                            end_pos = start_pos + 1
                            if end_pos < len(read_seq):
                                hap = read_seq[start_pos:end_pos]
                                if read_name not in read_haps:
                                    read_haps.setdefault(read_name, ["x"] * nvar)
                                if hap.upper() == allele1.upper():
                                    read_haps[read_name][dsnp_index] = "1"
                                elif hap.upper() == allele2.upper():
                                    read_haps[read_name][dsnp_index] = "2"
        return read_haps
    """

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
        new_reads = {}
        for i, hap in enumerate(reads):
            hap_reads = set()
            for read in reads[hap]:
                if "sup" not in read:
                    hap_reads.add(read)
                else:
                    hap_reads.add(read.split("_sup")[0])
            new_reads.setdefault(hap, hap_reads)
        d = {}
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
                        d.setdefault(hap1, []).append(hap2)
                        d.setdefault(hap2, []).append(hap1)
        d = dict(sorted(d.items(), key=lambda item: len(item[1]), reverse=True))
        # print(d)
        alleles = []
        if d != {}:
            alleles = [[list(d.keys())[0]] + list(d.values())[0]]
            for hap1 in d:
                for hap2 in d[hap1]:
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
        return alleles

    def check_gene_presence(self):
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

        # last snp outside of repeat
        # if self.candidate_pos != set():
        #    self.candidate_pos.add("32046300_G_A")

        var_found = False
        for var in self.candidate_pos:
            pos = int(var.split("_")[0])
            if pos > self.clip_3p_positions[0]:
                var_found = True
                break
        if var_found is False:
            self.candidate_pos.add("32046300_G_A")

        het_sites = sorted(list(self.candidate_pos))
        if "32029159_T_C" in het_sites:
            het_sites.remove("32029159_T_C")
        if "32022483_G_A" in het_sites:
            het_sites.remove("32022483_G_A")

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

        # pprint(read_counts)
        # get haps that extend into tnxb
        last_genes = [a for a in ass_haps if a[-1] != "0"]
        """
        last_genes = []
        ending_haps = {}
        for hap in uniquely_supporting_reads:
            ending_haps.setdefault(
                hap,
                [
                    a
                    for a in uniquely_supporting_reads[hap]
                    if raw_read_haps[a][-1] != "x"
                ],
            )
        for hap in ending_haps:
            if len(ending_haps[hap]) >= 2:
                last_genes.append(hap)
        if len(last_genes) == 0:
            for read in nonuniquely_supporting_reads:
                if raw_read_haps[read][-1] != "x":
                    for hap in nonuniquely_supporting_reads[read]:
                        ending_haps.setdefault(hap, []).append(read)
            # pprint(ending_haps)
            for hap in ending_haps:
                if len(ending_haps[hap]) >= 2 and hap not in last_genes:
                    last_genes.append(hap)
        """

        total_cn = None
        if ass_haps == [] and self.het_sites == []:
            # feed all reads to deepvariant and call in a diploid mode
            if self.has_gene1 is True and self.has_gene2 is False:
                total_cn = 2
            elif self.has_gene1 is False and self.has_gene2 is True:
                total_cn = 0
        else:
            total_cn = len(ass_haps)

        tmp = {}
        for i, hap in enumerate(ass_haps):
            tmp.setdefault(hap, f"hap{i+1}")
        final_haps = tmp

        haplotypes = None
        dvar = None
        if self.het_sites != []:
            haplotypes, dvar = self.output_variants_in_haplotypes(
                final_haps,
                uniquely_supporting_reads,
                nonuniquely_supporting_reads,
            )

        alleles = self.get_alleles(uniquely_supporting_reads)
        # print(alleles)
        # print(last_genes)
        # check wrong phasing
        wrong_allele = False
        for allele in alleles:
            hp_set = set(allele)
            for hp in hp_set:
                if hp in last_genes and allele.count(hp) > 1:
                    wrong_allele = True
                    # alleles = []
                    break
        if wrong_allele:
            alleles = []

        two_cp_haplotypes = []
        # 2 cp, one allele, one haplotype extends to tnxb -> each haplotype has cn 2
        # call deepvariant in diploid mode
        if len(last_genes) == 1 and total_cn == 2 and len(alleles) == 1:
            two_cp_haplotypes = ass_haps

        # 3 or 4 cp, only one haplotype extends to tnxb -> identify if this haplotype has twice coverage
        # check unassigned reads
        if (total_cn == 3 or len(last_genes) < 2) and read_counts is not None:
            # check if one smn1 haplotype has more reads than others
            haps = list(read_counts.keys())
            counts = list(read_counts.values())
            # print(counts)
            max_count = max(counts)
            cp2_hap = haps[counts.index(max_count)]
            if cp2_hap in ass_haps:
                others_max = sorted(counts, reverse=True)[1]
                probs = self.depth_prob(max_count, others_max)
                if probs[0] < 0.0005 and others_max >= 10:
                    two_cp_haplotypes.append(cp2_hap)

        # annotate haplotypes by checking the diff sites
        # output variants carried by each haplotype

        self.close_handle()

        return self.GeneCall(
            total_cn,
            final_haps,
            two_cp_haplotypes,
            last_genes,
            self.gene1_read_number,
            self.gene2_read_number,
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
        )
