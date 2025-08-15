# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


import os
import pysam
import subprocess
import random
import re
import math
from scipy.stats import poisson
from collections import Counter
from .haplotype_assembler import VariantGraph
import paraphase


class BamRealigner:
    """
    Extract and realign reads to region of interest
    """

    min_mapq = 50
    min_aln = 800
    max_mismatch = 0.05
    deletion = r"\d+D"
    insertion = r"\d+I"

    def __init__(self, bam, outdir, config, prog_cmd, sample_id, genome_reference):
        self.bam = bam
        self.outdir = outdir
        self.prog_cmd = prog_cmd
        self.sample_id = sample_id
        self.genome_reference = genome_reference
        self.gene = config["gene"]
        self.ref = config["data"]["reference"]
        self.nchr_old = config["nchr_old"]
        self.nchr = config["nchr"]
        self.offset = int(self.nchr_old.split("_")[1]) - 1
        self.nchr_length = config["nchr_length"]
        self.extract_regions = config["extract_regions"]
        self.samtools = config["tools"]["samtools"]
        self.minimap2 = config["tools"]["minimap2"]
        self.methyl = {}
        if "check_nm" in config:
            self.max_mismatch = config["check_nm"]
        self.use_r2k = ""
        if "use_r2k" in config:
            self.use_r2k = "-r2k"
        self.realign_bam = os.path.join(
            self.outdir, self.sample_id + f"_{self.gene}_realigned_old.bam"
        )
        self.realign_out_bam = os.path.join(
            self.outdir, self.sample_id + f"_{self.gene}_realigned.bam"
        )
        self.gene2_region = config.get("gene2_region")
        if self.gene2_region is not None:
            self.nchr_gene2 = config["nchr_gene2"]
            self.offset_gene2 = int(self.gene2_region.split(":")[1].split("-")[0]) - 1
            self.nchr_length_gene2 = config["nchr_length_gene2"]
            self.gene2_ref = config["data"].get("reference_gene2")
            self.realign_bam_gene2 = os.path.join(
                self.outdir, self.sample_id + f"_{self.gene}_gene2_realigned_old.bam"
            )
            self.realign_out_bam_gene2 = os.path.join(
                self.outdir, self.sample_id + f"_{self.gene}_gene2_realigned.bam"
            )

    def get_nm(self, read, min_size=300):
        """Get number of mismatches excluding big deletions or insertions"""
        cigar = read.cigarstring
        deletions_on_read = [int(a[:-1]) for a in re.findall(self.deletion, cigar)]
        large_deletions = [a for a in deletions_on_read if a > min_size]
        insertions_on_read = [int(a[:-1]) for a in re.findall(self.insertion, cigar)]
        large_insertions = [a for a in insertions_on_read if a > min_size]
        return read.get_tag("NM") - sum(large_deletions + large_insertions)

    def write_realign_bam(self, gene2=False):
        """
        Realign reads to region of interest and output a tagged bam for visualization
        """
        if gene2 is False:
            realign_ref = self.ref
            realign_out_tmp = self.realign_bam
            realign_out = self.realign_out_bam
            ref_name_tmp = self.nchr_old
            offset = self.offset
            nchr = self.nchr
            nchr_length = self.nchr_length
        else:
            realign_ref = self.gene2_ref
            realign_out_tmp = self.realign_bam_gene2
            realign_out = self.realign_out_bam_gene2
            ref_name_tmp = self.gene2_region.replace(":", "_").replace("-", "_")
            offset = self.offset_gene2
            nchr = self.nchr_gene2
            nchr_length = self.nchr_length_gene2

        wgs_bamh = pysam_handle(self.bam, reference_fasta=self.genome_reference)
        has_rq = False
        read = None
        region = self.extract_regions.split()[0]
        region_split = region.split(":")
        for read in wgs_bamh.fetch(
            region_split[0],
            int(region_split[1].split("-")[0]),
            int(region_split[1].split("-")[1]),
        ):
            if read.is_supplementary is False and read.is_secondary is False:
                break
        if read is not None:
            if read.has_tag("rq"):
                has_rq = True
        wgs_bamh.close()
        if has_rq:
            realign_cmd = (
                f'echo {self.extract_regions} | sed -e "s/ /\\n/g"  | sed -e "s/:/\\t/" | sed -e "s/-/\\t/" | '
                + f"{self.samtools} view --region-file - -u -F 0x100 -F 0x200 -F 0x800 -e '[rq]>=0.99' -T {self.genome_reference} {self.bam} | {self.samtools} fastq -T MM,ML - | {self.minimap2} {self.use_r2k} -ay -x map-pb {realign_ref} - | {self.samtools} view -bh | {self.samtools} sort -n | {self.samtools} sort > {realign_out_tmp}"
            )
        else:
            realign_cmd = (
                f'echo {self.extract_regions} | sed -e "s/ /\\n/g"  | sed -e "s/:/\\t/" | sed -e "s/-/\\t/" | '
                + f"{self.samtools} view --region-file - -u -F 0x100 -F 0x200 -F 0x800 -T {self.genome_reference} {self.bam} | {self.samtools} fastq -T MM,ML - | {self.minimap2} {self.use_r2k} -ay -x map-pb {realign_ref} - | {self.samtools} view -bh | {self.samtools} sort -n | {self.samtools} sort > {realign_out_tmp}"
            )
        result = subprocess.run(realign_cmd, capture_output=True, text=True, shell=True)
        result.check_returncode()
        if os.path.exists(realign_out_tmp) is False:
            raise Exception("Realigned bam does not exist.")

        pysam.index(realign_out_tmp)
        realign_bamh = pysam.AlignmentFile(realign_out_tmp, "rb")
        realign_bamh_header = realign_bamh.header
        realign_bamh_header = realign_bamh_header.to_dict()
        pg_lines = []
        pg_lines_original = realign_bamh_header.get("PG")
        if pg_lines_original is not None:
            pg_lines = [
                a for a in pg_lines_original if "ID" in a and a["ID"] == "minimap2"
            ]
        pg_lines.append(
            {
                "PN": "paraphase",
                "ID": "paraphase",
                "VN": paraphase.__version__,
                "CL": f"paraphase {self.prog_cmd}",
            }
        )
        new_header = {
            "SQ": [{"SN": nchr, "LN": nchr_length}],
            "PG": pg_lines,
        }
        realign_out_bamh = pysam.AlignmentFile(
            realign_out,
            "wb",
            header=new_header,
        )
        unique_reads_added = set()
        for read in realign_bamh.fetch(ref_name_tmp):
            num_mismatch = self.get_nm(read)
            if (
                read.mapping_quality >= self.min_mapq
                and read.query_alignment_length >= self.min_aln
                and num_mismatch < read.reference_length * self.max_mismatch
            ):
                if (
                    read.qname,
                    read.cigarstring,
                ) not in unique_reads_added:
                    unique_reads_added.add((read.qname, read.cigarstring))
                    read.reference_start += offset
                    ltags = read.tags
                    new_ltags = []
                    for tag in ltags:
                        if tag[0] != "SA":
                            new_ltags.append(tag)
                        else:
                            # update the SA tag
                            new_sa = []
                            for sa_item in tag[1].split(";"):
                                if sa_item != "":
                                    at = sa_item.split(",")
                                    at[0] = nchr
                                    tmp_pos = int(at[1])
                                    at[1] = str(tmp_pos + offset)
                                    mapq = int(at[4])
                                    if mapq >= self.min_mapq:
                                        new_at = ",".join(at)
                                        new_sa.append(new_at)
                            if new_sa != []:
                                new_sa.append("")
                                new_ltags.append(("SA", ";".join(new_sa)))
                    read.tags = new_ltags
                    realign_out_bamh.write(read)
        realign_bamh.close()
        realign_out_bamh.close()
        pysam.index(realign_out)
        os.remove(realign_out_tmp)
        os.remove(realign_out_tmp + ".bai")


class BamTagger:
    """
    Add HP tags to bam
    """

    read_color = "166,206,227"
    read_color_allele1 = "178,223,138"
    read_color_allele2 = "177,156,217"

    def __init__(self, sample_id, outdir, config, call_sum):
        self.sample_id = sample_id
        self.outdir = outdir
        self.call_sum = call_sum
        self.gene = config["gene"]
        self.bam = os.path.join(
            self.outdir, self.sample_id + f"_{self.gene}_realigned.bam"
        )
        self.nchr = config["nchr"]
        self.tmp_bam = os.path.join(
            self.outdir, self.sample_id + f"_{self.gene}_tmp.bam"
        )
        self.tagged_realigned_bam = os.path.join(
            self.outdir, self.sample_id + f"_{self.gene}_realigned_tagged.bam"
        )
        self.use_supplementary = False
        if "use_supplementary" in config or "is_tandem" in config:
            self.use_supplementary = True
        self.gene2_region = config.get("gene2_region")
        if self.gene2_region is not None:
            self.nchr_gene2 = config["nchr_gene2"]
            self.bam_gene2 = os.path.join(
                self.outdir, self.sample_id + f"_{self.gene}_gene2_realigned.bam"
            )
            self.tmp_bam_gene2 = os.path.join(
                self.outdir, self.sample_id + f"_{self.gene}_gene2_tmp.bam"
            )
            self.tagged_realigned_bam_gene2 = os.path.join(
                self.outdir, self.sample_id + f"_{self.gene}_gene2_realigned_tagged.bam"
            )
        random.seed(0)

    def add_tag_to_read(
        self,
        read,
        hp_keys,
        reads_to_tag,
        nonunique,
        read_details,
        alleles=[],
        random_assign=False,
        gene2=False,
    ):
        """
        Add HP tag to each read.
        Option to randomly assign a read to a possible haplotype when
        nonuniquelly matching. The randomly assigned reads are in gray
        while unique reads are in blue.
        """
        read.set_tag("RN", self.gene, "Z")
        hp_found = False
        read_name = read.qname
        if (
            read.is_supplementary
            and self.use_supplementary
            and gene2 is False
            # coordinates are all changed after realigning to gene2
        ):
            read_name = (
                read_name + f"_sup_{read.reference_start}_{read.reference_length}"
            )
        for hap, hap_name in hp_keys.items():
            if read_name in reads_to_tag[hap]:
                read.set_tag("HP", hap_name, "Z")
                read.set_tag("YC", self.read_color, "Z")
                if alleles != []:
                    if hap_name in alleles[0]:
                        read.set_tag("YC", self.read_color_allele1, "Z")
                    elif len(alleles) > 1 and hap_name in alleles[1]:
                        read.set_tag("YC", self.read_color_allele2, "Z")
                hp_found = True
        if hp_found is False:
            # find closest match
            if read_name in read_details:
                read_seq = read_details[read_name]
                keys = []
                mismatches = []
                for ass_hap in hp_keys.keys():
                    match, mismatch, extend = VariantGraph.compare_two_haps(
                        read_seq, ass_hap
                    )
                    # if match > 0:
                    keys.append(ass_hap)
                    mismatches.append(mismatch)
                mismatches_sorted = sorted(mismatches)
                if (
                    len(mismatches_sorted) > 1
                    and 0 < mismatches_sorted[0] <= 2
                    and mismatches_sorted[1] >= mismatches_sorted[0] + 2
                ):
                    best_match = keys[mismatches.index(mismatches_sorted[0])]
                    selected_hap = hp_keys[best_match]
                    read.set_tag("HP", selected_hap, "Z")
                    hp_found = True
                    read.set_tag("YC", self.read_color, "Z")
                    if alleles != []:
                        if selected_hap in alleles[0]:
                            read.set_tag("YC", self.read_color_allele1, "Z")
                        elif len(alleles) > 1 and selected_hap in alleles[1]:
                            read.set_tag("YC", self.read_color_allele2, "Z")
        if hp_found is False:
            if random_assign is False:
                read.set_tag("HP", "Unassigned", "Z")
            else:
                if read_name in nonunique:
                    possible_haps = nonunique[read_name]
                    random_hap = possible_haps[
                        random.randint(0, len(possible_haps) - 1)
                    ]
                    if random_hap in hp_keys:
                        read.set_tag("HP", hp_keys[random_hap], "Z")
                else:
                    read.set_tag("HP", "Unassigned", "Z")
        return read

    def write_bam(self, random_assign=False, gene2=False):
        """
        Write HP tags in output bams.
        """
        if gene2 is False:
            inbam = self.bam
            tmpbam = self.tmp_bam
            outbam = self.tagged_realigned_bam
            nchr = self.nchr
        else:
            inbam = self.bam_gene2
            tmpbam = self.tmp_bam_gene2
            outbam = self.tagged_realigned_bam_gene2
            nchr = self.nchr_gene2

        call_sum = self.call_sum
        hp_keys = call_sum.get("final_haplotypes")
        self._bamh = pysam.AlignmentFile(inbam, "rb")
        out_bamh = pysam.AlignmentFile(tmpbam, "wb", template=self._bamh)
        unique_reads = call_sum.get("unique_supporting_reads")
        nonunique_reads = call_sum.get("nonunique_supporting_reads")
        if nonunique_reads is None:
            nonunique_reads = {}
        read_details = call_sum.get("read_details")
        alleles = call_sum.get("raw_alleles")
        if alleles is None:
            alleles = []

        for read in self._bamh.fetch(nchr):
            if read.is_secondary is False:
                read.set_tag("RN", self.gene, "Z")
                if hp_keys is not None:
                    read = self.add_tag_to_read(
                        read,
                        hp_keys,
                        unique_reads,
                        nonunique_reads,
                        read_details,
                        alleles,
                        random_assign=random_assign,
                        gene2=gene2,
                    )
                out_bamh.write(read)
        out_bamh.close()
        self._bamh.close()
        pysam.sort("-o", outbam, tmpbam)
        pysam.index(outbam)
        os.remove(tmpbam)
        os.remove(inbam)
        os.remove(inbam + ".bai")


class VcfGenerater:
    """
    Call variants and generate individual/merged vcfs
    """

    search_range = 200
    min_base_quality_for_variant_calling = 25

    def __init__(self, sample_id, outdir, call_sum, args=None):
        self.sample_id = sample_id
        self.outdir = outdir
        self.call_sum = call_sum
        self.lowqual = False
        if args is not None:
            self.lowqual = args.write_nocalls_in_vcf
        self.match = {}

    def set_parameter(self, config, tmpdir=None, prog_cmd=None):
        self.gene = config["gene"]
        self.nchr = config["nchr"]
        self.nchr_old = config["nchr_old"]
        self.offset = int(self.nchr_old.split("_")[1]) - 1
        self.nchr_length = config["nchr_length"]
        self.ref = config["data"]["reference"]
        self.left_boundary = config.get("left_boundary")
        self.right_boundary = config.get("right_boundary")
        self.deletion1_in_gene1 = config.get("deletion1_in_gene1")
        self.deletion1_in_gene2 = config.get("deletion1_in_gene2")
        self.extract_regions = config.get("extract_regions")
        if self.left_boundary is None:
            self.left_boundary = int(self.nchr_old.split("_")[1])
        if self.right_boundary is None:
            self.right_boundary = int(self.nchr_old.split("_")[2])
        self.use_supplementary = False
        if "use_supplementary" in config or "is_tandem" in config:
            self.use_supplementary = True

        self.prog_cmd = prog_cmd
        self.tmpdir = tmpdir
        if self.tmpdir is None:
            self.tmpdir = self.outdir
        self.bam = os.path.join(
            tmpdir, self.sample_id + f"_{self.gene}_realigned_tagged.bam"
        )
        self.vcf_dir = os.path.join(self.outdir, f"{self.sample_id}_paraphase_vcfs")

    def get_range_in_other_gene(self, pos):
        """
        Find the correponding coordinates in the other gene
        """
        if pos in self.match:
            return self.match[pos]
        for i in range(self.search_range):
            new_pos = pos + i
            if new_pos in self.match:
                return self.match[new_pos]
        return None

    def write_header(self, fout, has_sv=False):
        """Write VCF header"""
        fout.write("##fileformat=VCFv4.2\n")
        fout.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        fout.write('##FILTER=<ID=LowQual,Description="Nonpassing variant">\n')
        fout.write(
            '##INFO=<ID=HPBOUND,Number=.,Type=String,Description="Boundary coordinates of the phased haplotype">\n'
        )
        alleles = self.call_sum.get("alleles_final")
        if alleles is not None and alleles != []:
            fout.write(
                '##INFO=<ID=ALLELE,Number=.,Type=String,Description="Haplotypes phased into alleles">\n'
            )
        if has_sv:
            fout.write(
                '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">\n'
            )
            fout.write(
                '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of the SV">\n'
            )
            fout.write(
                '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">\n'
            )
            fout.write('##ALT=<ID=DEL,Description="Deletion">\n')
            fout.write('##ALT=<ID=DUP,Description="Duplication">\n')
            fout.write('##ALT=<ID=INV,Description="Inversion">\n')
        fout.write(
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype per haplotype">\n'
        )
        fout.write(
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth per haplotype">\n'
        )
        fout.write(
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for REF and ALT per haplotype">\n'
        )
        fout.write(f"##contig=<ID={self.nchr},length={self.nchr_length}>\n")
        fout.write(f"##paraphase_version={paraphase.__version__}\n")
        fout.write(f"##paraphase_command=paraphase {self.prog_cmd}\n")

    @staticmethod
    def modify_hapbound(bound1, bound2, truncated):
        """Get haplotype boundaries to appear in vcf"""
        hap_bound = f"{bound1}-{bound2}"
        if truncated == ["5p"]:
            hap_bound = f"{bound1}truncated-{bound2}"
        elif truncated == ["3p"]:
            hap_bound = f"{bound1}-{bound2}truncated"
        elif truncated == ["5p", "3p"]:
            hap_bound = f"{bound1}truncated-{bound2}truncated"
        return hap_bound

    @staticmethod
    def convert_alt_record(ref, alt):
        """Convert variant ALT allele to the pileup format"""
        if len(alt) > 1:
            ins_len = len(alt) - 1
            return f"{ref}+{ins_len}{alt[1:]}"
        if len(ref) > 1:
            del_len = len(ref) - 1
            del_seq_n = "N" * del_len
            return f"{ref[0]}-{del_len}{del_seq_n}"
        return alt

    @staticmethod
    def check_sv_call_presence(vars_list):
        """Check if variant calls contain any SV"""
        for variants_info, _ in vars_list:
            for pos in variants_info:
                for var in variants_info[pos]:
                    if var is not None:
                        _, _, alt = var[0].split("_")
                        if alt.isdigit():
                            return True
        return False

    def merge_vcf(self, vars_list):
        """
        Merge vcfs from multiple haplotypes.
        """
        os.makedirs(self.vcf_dir, exist_ok=True)
        has_sv = self.check_sv_call_presence(vars_list)
        merged_vcf = os.path.join(self.vcf_dir, self.sample_id + f"_{self.gene}.vcf")
        with open(merged_vcf, "w") as fout:
            self.write_header(fout, has_sv=has_sv)
            assert len(vars_list) <= 2
            haps_ids = []
            haps_ids1 = []
            haps_ids2 = []
            haps_bounds = []
            for list_counter, (variants_info, haps_info) in enumerate(vars_list):
                for hap_name, bound1, bound2, truncated in haps_info:
                    haps_ids.append(hap_name)
                    hap_bound = self.modify_hapbound(bound1, bound2, truncated)
                    haps_bounds.append(hap_bound)
                    if list_counter == 0:
                        haps_ids1.append(hap_name)
                    else:
                        haps_ids2.append(hap_name)
            header = [
                "#CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "INFO",
                "FORMAT",
            ] + haps_ids
            fout.write("\t".join(header) + "\n")

            for list_counter, (variants_info, haps_info) in enumerate(vars_list):
                variants_info = dict(sorted(variants_info.items()))
                for pos in variants_info:
                    call_info = variants_info[pos]
                    # unique variants at this site
                    variant_observed = set([a[0] for a in call_info if a is not None])
                    ref_only = False
                    if len(variant_observed) == 1:
                        _, ref, alt = list(variant_observed)[0].split("_")
                        if ref == alt:
                            ref_only = True
                    if len(variant_observed) == 2:
                        has_ref = False
                        has_del = False
                        for var in variant_observed:
                            _, a, b = var.split("_")
                            if a == b:
                                has_ref = True
                            elif b == "*":
                                has_del = True
                        if has_ref and has_del:
                            ref_only = True

                    for variant in variant_observed:
                        _, ref, alt = variant.split("_")
                        valid_gts = []
                        merge_gt = []
                        merge_ad = []
                        merge_dp = []
                        for each_call_index, each_call in enumerate(call_info):
                            if each_call is None:
                                merge_gt.append(".")
                                merge_ad.append(".")
                                merge_dp.append(".")
                                hap_info = haps_info[each_call_index]
                                (
                                    _,
                                    hap_info_bound1,
                                    hap_info_bound2,
                                    hap_info_truncated,
                                ) = hap_info
                                if hap_info_truncated == []:
                                    valid_gts.append(".")
                                elif hap_info_truncated == ["5p"]:
                                    if pos > hap_info_bound1:
                                        valid_gts.append(".")
                                elif hap_info_truncated == ["3p"]:
                                    if pos < hap_info_bound2:
                                        valid_gts.append(".")
                                elif hap_info_truncated == ["5p", "3p"]:
                                    if hap_info_bound1 < pos < hap_info_bound2:
                                        valid_gts.append(".")
                            else:
                                var_name, dp, ad, var_filter, gt, counter = each_call
                                if counter is None:
                                    if ref != alt:
                                        this_ad = ",".join([str(a) for a in ad])
                                    else:
                                        this_ad = ",".join([str(a) for a in [ad[0], 0]])
                                else:
                                    if ref != alt:
                                        alt_converted = self.convert_alt_record(
                                            ref, alt
                                        )
                                        this_ad = ",".join(
                                            [
                                                str(a)
                                                for a in [
                                                    ad[0],
                                                    counter[alt_converted],
                                                ]
                                            ]
                                        )
                                    else:
                                        this_ad = ",".join(
                                            [str(a) for a in [ad[0], dp - ad[0]]]
                                        )
                                if var_filter != []:
                                    gt = "."
                                merge_dp.append(str(dp))
                                if gt == "0":
                                    merge_gt.append(gt)
                                    valid_gts.append(gt)
                                    merge_ad.append(this_ad)
                                elif var_name == variant:
                                    merge_gt.append(gt)
                                    valid_gts.append(gt)
                                    merge_ad.append(this_ad)
                                else:
                                    merge_gt.append(".")
                                    valid_gts.append(".")
                                    merge_ad.append(this_ad)
                        write_variant = False
                        if self.lowqual is True:
                            if (
                                (alt != ref or ref_only)
                                and alt not in [".", "*"]
                                and ("1" in merge_gt or "." in valid_gts)
                            ):
                                write_variant = True
                        elif alt != ref and alt not in [".", "*"] and "1" in merge_gt:
                            write_variant = True
                        final_qual = "."
                        if write_variant:
                            if list_counter == 0 and haps_ids != haps_ids1:
                                for _ in range(len(haps_ids2)):
                                    merge_gt.append(".")
                                    merge_ad.append(".")
                                    merge_dp.append(".")
                            elif list_counter > 0:
                                for _ in range(len(haps_ids1)):
                                    merge_gt.insert(0, ".")
                                    merge_ad.insert(0, ".")
                                    merge_dp.insert(0, ".")

                            if "1" in merge_gt:
                                variant_filter = "PASS"
                            else:
                                variant_filter = "LowQual"
                            info_field = "HPBOUND=" + ",".join(haps_bounds)
                            alleles = self.call_sum.get("alleles_final")
                            if alleles is not None and alleles != []:
                                alleles_rename = []
                                for allele in alleles:
                                    allele_rename = []
                                    for a in allele:
                                        if a is not None:
                                            allele_rename.append(a)
                                        else:
                                            allele_rename.append("unknown")
                                    alleles_rename.append(allele_rename)
                                info_field += (
                                    ";"
                                    + "ALLELE="
                                    + ",".join(["+".join(a) for a in alleles_rename])
                                )
                            if alt.isdigit() is False:
                                if alt == ref:
                                    alt = "."
                                merged_entry = [
                                    self.nchr,
                                    str(pos),
                                    ".",
                                    ref,
                                    alt,
                                    final_qual,
                                    variant_filter,
                                    info_field,
                                    "GT:DP:AD",
                                ] + [
                                    ":".join([merge_gt[j], merge_dp[j], merge_ad[j]])
                                    for j in range(len(haps_ids))
                                ]

                            else:
                                nstart, var_type, nend = variant.split("_")
                                nstart = int(nstart)
                                nend = int(nend)
                                var_size = nend - nstart
                                info_field = (
                                    f"SVTYPE={var_type};END={nend};SVLEN={var_size};"
                                    + info_field
                                )
                                merged_entry = [
                                    self.nchr,
                                    str(pos),
                                    ".",
                                    "N",
                                    f"<{var_type}>",
                                    final_qual,
                                    variant_filter,
                                    info_field,
                                    "GT:DP:AD",
                                ] + [
                                    ":".join([merge_gt[j], merge_dp[j], merge_ad[j]])
                                    for j in range(len(haps_ids))
                                ]
                            fout.write("\t".join(merged_entry) + "\n")

    @staticmethod
    def refine_indels(ref_seq, var_seq, pos, refh, ref_name):
        """process indels"""
        if "+" in var_seq:
            ins_base = var_seq.split(re.findall(r"\+\d+", var_seq)[0])[1]
            var_seq = ref_seq + ins_base
        elif "-" in var_seq:
            del_len = int(re.findall(r"\-\d+", var_seq)[0][1:])
            var_seq = ref_seq
            ref_seq = refh.fetch(
                ref_name,
                pos - 1,
                pos + del_len,
            ).upper()
        return ref_seq, var_seq

    @staticmethod
    def get_var_call_qual(dp, ad, gt, is_snp):
        "Calculate variant call quality. Not used"
        noise_ref = 0.03
        noise_alt = 0.05
        if is_snp is False:
            noise_alt = 0.1
        min_qual = 0
        max_qual = 100
        if gt == "0":
            # non-allele reads given noise
            h0 = poisson.pmf(dp - ad, dp * noise_ref)
            # allele reads given noise
            h1 = poisson.pmf(ad, dp * noise_alt)
        if gt == "1":
            # non-allele reads given noise
            h0 = poisson.pmf(dp - ad, dp * noise_alt)
            # allele reads given noise
            h1 = poisson.pmf(ad, dp * noise_alt)
        qual = math.floor(10 * math.log(h0 / h1, 10))
        return max(min(qual, max_qual), min_qual)

    @staticmethod
    def get_var(all_bases, ref_seq):
        """Get most supported base as variant"""
        dp = len(all_bases)
        gt = "."
        qual = "."
        ad = len([a for a in all_bases if a != ref_seq])
        ad_ref = len([a for a in all_bases if a == ref_seq])
        var_seq = ref_seq
        counter = None
        if all_bases != []:
            counter = Counter(all_bases)
            most_common_base = counter.most_common(2)
            var_seq = most_common_base[0][0]
            ad = most_common_base[0][1]
            is_snp = False
            if len(var_seq) == 1 and len(ref_seq) == 1 and var_seq != "*":
                is_snp = True
            if var_seq == ref_seq:
                gt = "0"
            else:
                gt = "1"
            # qual = VcfGenerater.get_var_call_qual(dp, ad, gt, is_snp)
        return [var_seq, dp, (ad_ref, ad), gt, qual, counter]

    def get_hap_bound(self, hap_name):
        """Get haplotype boundaries"""
        hap_bound = list(self.call_sum["haplotype_details"][hap_name]["boundary"])
        is_truncated = self.call_sum["haplotype_details"][hap_name]["is_truncated"]
        # find the positions next to the existing boundaries
        confident_position = hap_bound[0]
        if confident_position == self.left_boundary or "5p" in is_truncated:
            hap_bound.append(confident_position)
        else:
            for var in self.call_sum["sites_for_phasing"]:
                confident_position = int(var.split("_")[0])
                if confident_position > hap_bound[0]:
                    break
            hap_bound.append(confident_position)
        confident_position = hap_bound[1]
        if confident_position == self.right_boundary or "3p" in is_truncated:
            hap_bound.append(confident_position)
        else:
            for var in reversed(self.call_sum["sites_for_phasing"]):
                confident_position = int(var.split("_")[0])
                if confident_position < hap_bound[1]:
                    break
            hap_bound.append(confident_position)
        return hap_bound

    def pileup_to_variant(
        self,
        pileups_raw,
        read_names,
        uniq_reads,
        refh,
        offset,
        hap_bound,
        min_depth=4,
        min_qual=25,
        variants_to_add={},
    ):
        """
        Filter pileups and make variant calls.
        """
        variants = []
        for pos in variants_to_add:
            if pos not in pileups_raw:
                var_name = variants_to_add[pos]
                nstart, var_type, nend = var_name.split("_")
                nstart = int(nstart)
                nend = int(nend)
                variants.append([nstart, var_name, ".", ".", [], "1", None])

        ref_name = refh.references[0]
        sv_pos = []
        for pos in pileups_raw:
            if pos in variants_to_add:
                var_name = variants_to_add[pos]
                nstart, var_type, nend = var_name.split("_")
                nstart = int(nstart)
                nend = int(nend)
                variants.append([nstart, var_name, ".", ".", [], "1", None])
                sv_pos = [nstart, nend]

            all_bases = pileups_raw[pos]
            if offset < 0:
                true_pos = pos
                refh_pos = pos + offset
            else:
                true_pos = pos + offset
                refh_pos = pos
            ref_seq = refh.fetch(ref_name, refh_pos - 1, refh_pos).upper()
            alt_all_reads = self.get_var(all_bases, ref_seq)
            var_seq, dp, ad, gt, qual, counter = alt_all_reads
            if (
                hap_bound == []
                or (None not in hap_bound and hap_bound[0] < true_pos < hap_bound[1])
            ) and (
                sv_pos == []
                or (
                    self.gene == "ikbkg"
                    and (true_pos < sv_pos[0] or true_pos > sv_pos[1])
                )
                or (self.gene != "ikbkg" and true_pos != sv_pos[0])
            ):
                # use only unique reads for positions at the edge
                # or no-call when using all reads
                if (
                    hap_bound == []
                    or true_pos < hap_bound[2]
                    or true_pos > hap_bound[3]
                    or ad[1] < dp * 0.7
                ):
                    bases_uniq_reads = []
                    for i, read_base in enumerate(all_bases):
                        if uniq_reads is None or read_names[pos][i] in uniq_reads:
                            bases_uniq_reads.append(read_base)
                    alt_uniq_reads = self.get_var(bases_uniq_reads, ref_seq)
                    var_seq, dp, ad, gt, qual, counter = alt_uniq_reads
                    if dp < min_depth or ad[1] < dp * 0.7:
                        var_seq, dp, ad, gt, qual, counter = alt_all_reads

                ref_seq, var_seq = self.refine_indels(
                    ref_seq, var_seq, refh_pos, refh, ref_name
                )
                var = f"{true_pos}_{ref_seq}_{var_seq}"
                var_filter = []
                if dp < min_depth:
                    var_filter.append("LowDP")
                if (gt == "1" and ad[1] < dp * 0.7) or (gt == "0" and ad[0] < dp * 0.7):
                    var_filter.append("LowQual")
                if var_filter != []:
                    gt = "."
                variants.append([true_pos, var, dp, ad, var_filter, gt, counter])
        return variants

    def get_sv(self, gene2=False):
        """Get SVs to report for certain genes"""
        call_sum = self.call_sum
        special_variants = {}
        if self.gene == "ikbkg":
            del_haps = call_sum.get("deletion_haplotypes")
            if del_haps is not None and del_haps != []:
                for del_hap in del_haps:
                    del_name = self.deletion1_in_gene1
                    if "pseudo" in del_hap and gene2 is True:
                        del_name = self.deletion1_in_gene2
                    special_variants.setdefault(del_hap, del_name)
        if self.gene == "f8":
            sv_called = call_sum.get("sv_called")
            if sv_called is not None and sv_called != {}:
                (
                    extract_region1,
                    extract_region2,
                    extract_region3,
                ) = self.extract_regions.split()
                extract_region1_end = extract_region1.split("-")[1]
                extract_region2_start = extract_region2.split(":")[1].split("-")[0]
                extract_region3_start = extract_region3.split(":")[1].split("-")[0]
                for sv_hap, sv in sv_called.items():
                    if sv == "inversion":
                        sv_name = f"{extract_region1_end}_INV_{extract_region3_start}"
                        special_variants.setdefault(sv_hap, sv_name)
                    elif sv == "deletion":
                        sv_name = f"{extract_region1_end}_DEL_{extract_region2_start}"
                        special_variants.setdefault(sv_hap, sv_name)
        """
        # add fusions, disabled for now
        fusion_calls = call_sum.get("fusions_called")
        if fusion_calls is not None and fusion_calls != {}:
            for sv_hap, sv_info in fusion_calls.items():
                sv_type = sv_info["type"]
                bp1 = sv_info["breakpoint"][0][0]
                bp2 = sv_info["breakpoint"][1][0]
                if sv_type == "deletion":
                    sv_name = f"{bp1}_DEL_{bp2}"
                    special_variants.setdefault(sv_hap, sv_name)
                elif sv_type == "duplication":
                    sv_name = f"{bp1}_DUP_{bp2}"
                    special_variants.setdefault(sv_hap, sv_name)
        """
        return special_variants

    def run_without_realign(
        self,
        gene2=False,
        final_haps={},
        match_range=False,
    ):
        """
        Make vcf from existing alignment.
        This works for gene/pseudogene scenarios,
        i.e. no need to realign to pseudogene reference.
        """
        call_sum = self.call_sum
        if gene2 is False:
            final_haps = call_sum.get("final_haplotypes")
        # special SV type variants to report for certain genes
        special_variants = self.get_sv(gene2)

        uniq_reads = []
        for read_set in self.call_sum["unique_supporting_reads"].values():
            for read_name in read_set:
                read_name_split = read_name.split("_sup")
                # supplementary alignments
                if self.use_supplementary and len(read_name_split) > 1:
                    uniq_reads.append(read_name_split[0])
                else:
                    uniq_reads.append(read_name)
        variants_info = {}
        two_cp_haplotypes = self.call_sum.get("two_copy_haplotypes")
        nhap = len(final_haps)
        if two_cp_haplotypes is not None:
            nhap += len([a for a in two_cp_haplotypes if a in final_haps.values()])
        hap_info = []

        # gene1only, or two-gene mode but gene1 side
        if gene2 is False or match_range is False:
            bamh = pysam.AlignmentFile(self.bam, "rb")
            refh = pysam.FastaFile(self.ref)
            offset = self.offset
            nchr = self.nchr
        else:
            bamh = pysam.AlignmentFile(self.bam_gene2, "rb")
            refh = pysam.FastaFile(self.ref_gene2)
            offset = self.offset_gene2
            nchr = self.nchr_gene2

        if (
            gene2 is False
            and final_haps == {}
            and "heterozygous_sites" in call_sum
            and call_sum["heterozygous_sites"] == []
        ):
            hap_name = f"{self.gene}_homozygous_hap1"
            pileups_raw = {}
            read_names = {}

            for pos in range(self.left_boundary, self.right_boundary):
                pileups_raw.setdefault(pos, [])
                read_names.setdefault(pos, [])

            for pileupcolumn in bamh.pileup(
                nchr,
                truncate=True,
                min_base_quality=self.min_base_quality_for_variant_calling,
            ):
                pos = pileupcolumn.pos + 1
                this_pos_bases = [
                    a.upper() for a in pileupcolumn.get_query_sequences(add_indels=True)
                ]
                pileups_raw[pos] = this_pos_bases
                read_names[pos] = pileupcolumn.get_query_names()
            variants_called = self.pileup_to_variant(
                pileups_raw,
                read_names,
                None,
                refh,
                0 - offset,
                [],
            )
            hap_info.append([hap_name, self.left_boundary, self.right_boundary, None])
            hap_info.append(
                [hap_name + "_cp2", self.left_boundary, self.right_boundary, None]
            )

            for pos, var_name, dp, ad, var_filter, gt, counter in variants_called:
                if pos > self.left_boundary and pos < self.right_boundary:
                    variants_info.setdefault(
                        pos,
                        [
                            [var_name, dp, ad, var_filter, gt, counter],
                            [var_name, dp, ad, var_filter, gt, counter],
                        ],
                    )

        i = 0
        for hap_name in final_haps.values():
            variants_to_add = {}
            if hap_name in special_variants:
                variant_to_add = special_variants[hap_name]
                variants_to_add.setdefault(
                    int(variant_to_add.split("_")[0]), variant_to_add
                )

            hap_bound = self.get_hap_bound(hap_name)
            # convert to positions in the other gene
            if match_range:
                if self.match == {}:
                    hap_bound = []
                else:
                    n1 = self.get_range_in_other_gene(hap_bound[0])
                    n2 = self.get_range_in_other_gene(hap_bound[1])
                    n3 = self.get_range_in_other_gene(hap_bound[2])
                    n4 = self.get_range_in_other_gene(hap_bound[3])
                    if None in [n1, n2, n3, n4]:
                        hap_bound = []
                    else:
                        hap_bound = [
                            min(n1, n2),
                            max(n1, n2),
                            min(n3, n4),
                            max(n3, n4),
                        ]
            if gene2 is False or match_range is False:
                this_hap_info = [
                    hap_name,
                    hap_bound[0],
                    hap_bound[1],
                    self.call_sum["haplotype_details"][hap_name]["is_truncated"],
                ]
                hap_info.append(this_hap_info)
            else:
                this_hap_info = [
                    hap_name,
                    hap_bound[0],
                    hap_bound[1],
                    None,
                ]
                hap_info.append(this_hap_info)

            # by HP tag
            pileups_raw = {}
            read_names = {}

            if hap_bound != []:
                for pos in range(hap_bound[0], hap_bound[1]):
                    pileups_raw.setdefault(pos, [])
                    read_names.setdefault(pos, [])

            for pileupcolumn in bamh.pileup(
                nchr,
                truncate=True,
                min_base_quality=self.min_base_quality_for_variant_calling,
            ):
                pos = pileupcolumn.pos + 1
                this_pos_bases = [
                    a.upper() for a in pileupcolumn.get_query_sequences(add_indels=True)
                ]
                this_position_hps = []
                this_pos_read_names = []
                this_pos_read_names_sup = []

                for read in pileupcolumn.pileups:
                    read_tag = read.alignment.get_tag("HP")
                    this_position_hps.append(read_tag)
                    read_name = read.alignment.query_name
                    new_read_name = read_name
                    if self.use_supplementary and read.alignment.is_supplementary:
                        new_read_name = (
                            read_name + f"_sup_{read.alignment.reference_start}"
                        )
                    this_pos_read_names.append(read_name)
                    this_pos_read_names_sup.append(new_read_name)

                for read_num, read_hap in enumerate(this_position_hps):
                    if read_hap == hap_name:
                        pileups_raw.setdefault(pos, []).append(this_pos_bases[read_num])
                        if self.use_supplementary:
                            read_names.setdefault(pos, []).append(
                                this_pos_read_names_sup[read_num]
                            )
                        else:
                            read_names.setdefault(pos, []).append(
                                this_pos_read_names[read_num]
                            )

            variants_called = self.pileup_to_variant(
                pileups_raw,
                read_names,
                uniq_reads,
                refh,
                0 - offset,
                hap_bound,
                variants_to_add=variants_to_add,
            )

            for pos, var_name, dp, ad, var_filter, gt, counter in variants_called:
                variants_info.setdefault(pos, [None] * nhap)
                variants_info[pos][i] = [var_name, dp, ad, var_filter, gt, counter]
                if hap_name in two_cp_haplotypes:
                    variants_info[pos][i + 1] = [
                        var_name,
                        dp,
                        ad,
                        var_filter,
                        gt,
                        counter,
                    ]
            if hap_name in two_cp_haplotypes:
                i += 1
                this_hap_info_cp2 = [this_hap_info[0] + "_cp2"] + this_hap_info[1:]
                hap_info.append(this_hap_info_cp2)
            i += 1

        bamh.close()
        refh.close()
        return variants_info, hap_info

    def run(self):
        final_haps = self.call_sum.get("final_haplotypes")
        true_het_sites = self.call_sum.get("heterozygous_sites")
        if final_haps is None:
            return
        if final_haps == {} and true_het_sites is not None and true_het_sites != []:
            return
        variants_info, hap_info = self.run_without_realign()
        self.merge_vcf([(variants_info, hap_info)])


class TwoGeneVcfGenerater(VcfGenerater):
    """
    Make vcf for two-gene scenario
    """

    def __init__(self, sample_id, outdir, call_sum, args):
        VcfGenerater.__init__(self, sample_id, outdir, call_sum, args)

    def set_parameter(self, config, tmpdir=None, prog_cmd=None):
        super().set_parameter(config, tmpdir, prog_cmd)
        self.gene2_region = config["gene2_region"]
        self.nchr_gene2 = config["nchr_gene2"]
        self.offset_gene2 = int(self.gene2_region.split(":")[1].split("-")[0]) - 1
        self.ref_gene2 = config["data"]["reference_gene2"]
        self.bam_gene2 = os.path.join(
            tmpdir, self.sample_id + f"_{self.gene}_gene2_realigned_tagged.bam"
        )
        self.position_match = config["data"].get("gene_position_match")
        if self.position_match is not None:
            with open(self.position_match) as f:
                for line in f:
                    at = line.split()
                    self.match.setdefault(int(at[0]), int(at[1]))

    def separate_two_genes(self):
        """Get haplotypes for gene1 and gene2"""
        all_haps = self.call_sum.get("final_haplotypes")
        gene1_haps = {}
        gene2_haps = {}
        if self.gene == "smn1":
            return self.call_sum["smn1_haplotypes"], {
                **self.call_sum["smn2_haplotypes"],
                **self.call_sum["smn_del78_haplotypes"],
            }
        elif self.gene == "pms2":
            for hap, hap_name in all_haps.items():
                if "pms2cl" not in hap_name:
                    gene1_haps.setdefault(hap, hap_name)
                else:
                    gene2_haps.setdefault(hap, hap_name)
        elif self.gene == "strc":
            for hap, hap_name in all_haps.items():
                if "strcp1" not in hap_name:
                    gene1_haps.setdefault(hap, hap_name)
                else:
                    gene2_haps.setdefault(hap, hap_name)
        elif self.gene == "ncf1":
            for hap, hap_name in all_haps.items():
                if "pseudo" not in hap_name:
                    gene1_haps.setdefault(hap, hap_name)
                else:
                    gene2_haps.setdefault(hap, hap_name)
        elif self.gene == "ikbkg":
            for hap, hap_name in all_haps.items():
                if "pseudo" not in hap_name:
                    gene1_haps.setdefault(hap, hap_name)
                else:
                    gene2_haps.setdefault(hap, hap_name)
        return gene1_haps, gene2_haps

    def run(self):
        """
        Process haplotypes one by one. Realign to different ref sequence
        in this two-gene scenario
        """
        call_sum = self.call_sum
        final_haps = self.call_sum.get("final_haplotypes")
        true_het_sites = self.call_sum.get("heterozygous_sites")
        if final_haps is None:
            return
        if final_haps == {} and true_het_sites is not None and true_het_sites != []:
            return
        gene1_haps, gene2_haps = self.separate_two_genes()
        vars_gene1, gene1_hap_info = self.run_without_realign(
            gene2=True,
            final_haps=gene1_haps,
        )
        vars_gene2, gene2_hap_info = self.run_without_realign(
            gene2=True,
            final_haps=gene2_haps,
            match_range=True,
        )
        if (
            vars_gene1 != {}
            and vars_gene2 != {}
            and list(vars_gene1.keys())[0] < list(vars_gene2.keys())[0]
        ):
            self.merge_vcf(
                [
                    (vars_gene1, gene1_hap_info),
                    (vars_gene2, gene2_hap_info),
                ]
            )
        else:
            self.merge_vcf(
                [
                    (vars_gene2, gene2_hap_info),
                    (vars_gene1, gene1_hap_info),
                ]
            )


def pysam_handle(input_file, reference_fasta=None):
    """Get pysam handle for bam/cram files"""
    if input_file.lower().endswith("cram"):
        return pysam.AlignmentFile(input_file, "rc", reference_filename=reference_fasta)
    return pysam.AlignmentFile(input_file, "rb")
