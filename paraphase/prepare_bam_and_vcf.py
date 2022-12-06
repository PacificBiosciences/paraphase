# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


import os
import pysam
import subprocess
import random
import re
from collections import Counter


class BamRealigner:
    """
    Extract and realign reads to region of interest
    """

    min_mapq = 50
    min_aln = 400

    def __init__(self, bam, outdir, config):
        self.bam = bam
        self.outdir = outdir
        self.nchr = config["coordinates"]["hg38"]["nchr"]
        self.ref = config["data"]["reference"]
        self.nchr_old = config["coordinates"]["hg38"]["nchr_old"]
        self.offset = int(self.nchr_old.split("_")[1]) - 1
        self.nchr_length = config["coordinates"]["hg38"]["nchr_length"]
        self.extract_region1 = config["coordinates"]["hg38"]["extract_region1"]
        self.extract_region2 = config["coordinates"]["hg38"]["extract_region2"]
        self.samtools = config["tools"]["samtools"]
        self.minimap2 = config["tools"]["minimap2"]
        self._bamh = pysam.AlignmentFile(bam, "rb")
        self.sample_id = bam.split("/")[-1].split(".")[0]
        self.realign_bam = os.path.join(
            self.outdir, self.sample_id + "_realigned_old.bam"
        )
        self.realign_out_bam = os.path.join(
            self.outdir, self.sample_id + "_realigned.bam"
        )

    def write_realign_bam(self):
        """
        Realign reads to region of interest and output a tagged bam for visualization
        """
        realign_cmd = (
            f"{self.samtools} view -F 0x100 -F 0x200 -F 0x800 {self.bam} {self.extract_region1} {self.extract_region2} | "
            + f'awk \'BEGIN {{FS="\\t"}} {{print "@" $1 "\\n" $10 "\\n+\\n" $11}}\''
            + f" | {self.minimap2} -a -x map-pb {self.ref} - | {self.samtools} view -b | {self.samtools} sort > {self.realign_bam}"
        )
        result = subprocess.run(realign_cmd, capture_output=True, text=True, shell=True)
        result.check_returncode()
        if os.path.exists(self.realign_bam) is False:
            raise Exception("Realigned bam does not exist.")

        pysam.index(self.realign_bam)
        realign_bamh = pysam.AlignmentFile(self.realign_bam, "rb")
        realign_out_bamh = pysam.AlignmentFile(
            self.realign_out_bam,
            "wb",
            reference_names=[self.nchr],
            reference_lengths=[self.nchr_length],
        )
        for read in realign_bamh.fetch(self.nchr_old):
            if (
                read.mapping_quality >= self.min_mapq
                and read.query_alignment_length >= self.min_aln
            ):
                read.reference_start += self.offset
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
                                at[0] = self.nchr
                                tmp_pos = int(at[1])
                                at[1] = str(tmp_pos + self.offset)
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
        pysam.index(self.realign_out_bam)
        os.remove(self.realign_bam)
        os.remove(self.realign_bam + ".bai")


class BamTagger:
    """
    Add HP tags to bam
    """

    read_color = "166,206,227"

    def __init__(self, sample_id, outdir, config, call_sum):
        self.sample_id = sample_id
        self.outdir = outdir
        self.call_sum = call_sum
        self.bam = os.path.join(outdir, self.sample_id + "_realigned.bam")
        self.nchr = config["coordinates"]["hg38"]["nchr"]
        self._bamh = pysam.AlignmentFile(self.bam, "rb")
        self.tmp_bam = os.path.join(self.outdir, self.sample_id + "_tmp.bam")
        self.tagged_realigned_bam = os.path.join(
            self.outdir, self.sample_id + "_realigned_tagged.bam"
        )
        random.seed(0)

    def add_tag_to_read(
        self,
        read,
        hp_keys,
        reads_to_tag,
        nonunique,
        random_assign=False,
    ):
        """
        Add HP tag to each read.
        Option to randomly assign a read to a possible haplotype when
        nonuniquelly matching. The randomly assigned reads are in gray
        while unique reads are in blue.
        """
        hp_found = False
        for hap, hap_name in hp_keys.items():
            if read.qname in reads_to_tag[hap]:
                read.set_tag("HP", hap_name, "Z")
                if random_assign:
                    read.set_tag("YC", self.read_color, "Z")
                hp_found = True
        if hp_found is False:
            if random_assign is False:
                read.set_tag("HP", "Unassigned", "Z")
            else:
                if read.qname in nonunique:
                    possible_haps = nonunique[read.qname]
                    random_hap = possible_haps[
                        random.randint(0, len(possible_haps) - 1)
                    ]
                    if random_hap in hp_keys:
                        read.set_tag("HP", hp_keys[random_hap], "Z")
                else:
                    read.set_tag("HP", "Unassigned", "Z")
        return read

    def write_bam(self, random_assign=False):
        """
        Write HP tags in output bams.
        """
        call_sum = self.call_sum
        out_bamh = pysam.AlignmentFile(self.tmp_bam, "wb", template=self._bamh)
        unique_reads = call_sum.get("unique_supporting_reads")
        nonunique_reads = call_sum.get("nonunique_supporting_reads")
        if nonunique_reads is None:
            nonunique_reads = {}
        hp_keys = call_sum.get("final_haplotypes")

        for read in self._bamh.fetch(self.nchr):
            if read.is_secondary is False:  # and read.is_supplementary is False:
                read = self.add_tag_to_read(
                    read,
                    hp_keys,
                    unique_reads,
                    nonunique_reads,
                    random_assign=random_assign,
                )
                out_bamh.write(read)
        out_bamh.close()
        self._bamh.close()
        pysam.sort("-o", self.tagged_realigned_bam, self.tmp_bam)
        pysam.index(self.tagged_realigned_bam)
        os.remove(self.tmp_bam)
        os.remove(self.bam)
        os.remove(self.bam + ".bai")


class VcfGenerater:
    """
    Call variants and generate individual/merged vcfs
    """

    search_range = 200

    def __init__(self, sample_id, outdir, config, call_sum):
        self.sample_id = sample_id
        self.outdir = outdir
        self.call_sum = call_sum
        self.bam = os.path.join(outdir, self.sample_id + "_realigned_tagged.bam")
        self.nchr = config["coordinates"]["hg38"]["nchr"]
        self.nchr_old = config["coordinates"]["hg38"]["nchr_old"]
        self.offset = int(self.nchr_old.split("_")[1]) - 1
        self.nchr_length = config["coordinates"]["hg38"]["nchr_length"]
        self.ref = config["data"]["reference"]
        self.samtools = config["tools"]["samtools"]
        self.minimap2 = config["tools"]["minimap2"]
        self.match = {}

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

    def write_header(self, fout):
        """Write VCF header"""
        fout.write("##fileformat=VCFv4.2\n")
        fout.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        fout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fout.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
        fout.write(
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">\n'
        )
        fout.write(f"##contig=<ID={self.nchr},length={self.nchr_length}>\n")
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
            "default",
        ]
        fout.write("\t".join(header) + "\n")

    def merge_vcf(self, vars_list):
        """
        Merge vcfs from multiple haplotypes.
        """
        merged_vcf = os.path.join(self.outdir, self.sample_id + f"_variants.vcf")
        with open(merged_vcf, "w") as fout:
            self.write_header(fout)
            for vars in vars_list:
                vars = dict(sorted(vars.items()))
                for pos in vars:
                    call_info = vars[pos]
                    variant_observed = set([a[0] for a in call_info if a is not None])
                    for variant in variant_observed:
                        _, ref, alt = variant.split("_")
                        merge_gt = []
                        merge_ad = []
                        merge_dp = []
                        for each_call in call_info:
                            if each_call is None:
                                merge_gt.append(".")
                                merge_ad.append(".")
                                merge_dp.append(".")
                            else:
                                var_name, dp, ad, qual, gt = each_call
                                merge_dp.append(str(dp))
                                if gt == "0":
                                    merge_gt.append(gt)
                                    merge_ad.append(str(dp - ad))
                                elif var_name == variant:
                                    merge_gt.append(gt)
                                    merge_ad.append(str(ad))
                                else:
                                    merge_gt.append(".")
                                    if gt == ".":
                                        merge_ad.append(str(ad))
                                    else:
                                        merge_ad.append(".")
                        final_qual = "."
                        if "1" in merge_gt:
                            merged_entry = [
                                self.nchr,
                                str(pos),
                                ".",
                                ref,
                                alt,
                                final_qual,
                                "PASS",
                                ".",
                                "GT:DP:AD",
                                "/".join(merge_gt)
                                + ":"
                                + ",".join(merge_dp)
                                + ":"
                                + ",".join(merge_ad),
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
            )
        return ref_seq, var_seq

    @staticmethod
    def get_var(all_bases, ref_seq):
        """Get most supported base as variant"""
        dp = len(all_bases)
        gt = "."
        ad = len([a for a in all_bases if a != ref_seq])
        var_seq = ref_seq
        if all_bases != []:
            counter = Counter(all_bases)
            most_common_base = counter.most_common(2)
            var_seq = most_common_base[0][0]
            ad = most_common_base[0][1]
            if (
                dp >= 3
                and ad >= 2
                and (len(most_common_base) == 1 or ad > most_common_base[1][1])
            ):
                if var_seq == ref_seq or var_seq == "*":
                    gt = "0"
                else:
                    gt = "1"
        return [var_seq, dp, ad, gt]

    def read_pileup(self, hap_bam, hap_vcf_out, hap_bound, offset, ref, uniq_reads):
        """
        Call variants from bam. Take the most supported base at each position.
        """
        bamh = pysam.AlignmentFile(hap_bam, "rb")
        refh = pysam.FastaFile(ref)
        ref_name = refh.references[0]
        vcf_out = open(hap_vcf_out, "w")
        self.write_header(vcf_out)
        pileups_raw = {}
        read_names = {}
        variants = []
        for pileupcolumn in bamh.pileup(
            ref_name,
            truncate=True,
            min_base_quality=30,
        ):
            pos = pileupcolumn.pos + 1
            pileups_raw.setdefault(
                pos,
                [a.upper() for a in pileupcolumn.get_query_sequences(add_indels=True)],
            )
            read_names.setdefault(
                pos,
                pileupcolumn.get_query_names(),
            )
        for pos in pileups_raw:
            all_bases = pileups_raw[pos]
            true_pos = pos + offset
            ref_seq = refh.fetch(ref_name, pos - 1, pos)
            alt_all_reads = self.get_var(all_bases, ref_seq)
            if None not in hap_bound and hap_bound[0] < true_pos < hap_bound[1]:
                # use only unique reads for positions at the edge
                if true_pos < hap_bound[2] or true_pos > hap_bound[3]:
                    bases_uniq_reads = []
                    for i, read_base in enumerate(all_bases):
                        if read_names[pos][i] in uniq_reads:
                            bases_uniq_reads.append(read_base)
                    alt_uniq_reads = self.get_var(bases_uniq_reads, ref_seq)
                    if alt_uniq_reads[-1] != ".":
                        var_seq, dp, ad, gt = alt_uniq_reads
                    else:
                        var_seq, dp, ad, gt = alt_all_reads
                        gt = "."
                else:
                    var_seq, dp, ad, gt = alt_all_reads

                ref_seq, var_seq = self.refine_indels(
                    ref_seq, var_seq, pos, refh, ref_name
                )
                var = f"{true_pos}_{ref_seq}_{var_seq}"
                qual = "."
                variants.append([true_pos, var, dp, ad, qual, gt])
                if gt == "1":
                    vcf_out_line = [
                        self.nchr,
                        str(true_pos),
                        ".",
                        ref_seq,
                        var_seq,
                        str(qual),
                        "PASS",
                        ".",
                        "GT:DP:AD",
                        f"1:{dp}:{ad}",
                    ]
                    vcf_out.write("\t".join(vcf_out_line) + "\n")
        vcf_out.close()
        bamh.close()
        refh.close()
        return variants

    def run_step(
        self,
        final_haps,
        haplotype_details,
        ref_seq,
        offset,
        match_range=False,
    ):
        """
        Process haplotypes
        """
        uniq_reads = []
        for read_set in self.call_sum["unique_supporting_reads"].values():
            uniq_reads += read_set
        vars = {}
        two_cp_haplotypes = self.call_sum.get("two_copy_haplotypes")
        nhap = len(final_haps) + len(
            [a for a in two_cp_haplotypes if a in final_haps.values()]
        )
        i = 0
        for hap_name in final_haps.values():
            hap_bound = list(haplotype_details[hap_name]["boundary"])
            # find the positions next to the existing boundaries
            for var in self.call_sum["sites_for_phasing"]:
                confident_position = int(var.split("_")[0])
                if confident_position > hap_bound[0]:
                    break
            hap_bound.append(confident_position)
            for var in reversed(self.call_sum["sites_for_phasing"]):
                confident_position = int(var.split("_")[0])
                if confident_position < hap_bound[1]:
                    break
            hap_bound.append(confident_position)
            # convert to positions in the other gene
            if match_range:
                hap_bound = [
                    self.get_range_in_other_gene(hap_bound[0]),
                    self.get_range_in_other_gene(hap_bound[1]),
                    self.get_range_in_other_gene(hap_bound[2]),
                    self.get_range_in_other_gene(hap_bound[3]),
                ]
            hap_bam = os.path.join(self.outdir, self.sample_id + f"_{hap_name}.bam")
            realign_cmd = (
                f"{self.samtools} view -d HP:{hap_name} {self.bam} |"
                + f'awk \'BEGIN {{FS="\\t"}} {{print "@" $1 "\\n" $10 "\\n+\\n" $11}}\''
                + f" | {self.minimap2} -a -x map-pb {ref_seq} - | {self.samtools} view -b | {self.samtools} sort > {hap_bam}"
            )
            result = subprocess.run(
                realign_cmd, capture_output=True, text=True, shell=True
            )
            result.check_returncode()
            pysam.index(hap_bam)

            # call variants
            hap_vcf_out = os.path.join(self.outdir, self.sample_id + f"_{hap_name}.vcf")
            variants_called = self.read_pileup(
                hap_bam, hap_vcf_out, hap_bound, offset, ref_seq, uniq_reads
            )
            os.remove(hap_bam)
            os.remove(hap_bam + ".bai")

            for pos, var_name, dp, ad, qual, gt in variants_called:
                vars.setdefault(pos, [None] * nhap)
                vars[pos][i] = [var_name, dp, ad, qual, gt]
                if hap_name in two_cp_haplotypes:
                    vars[pos][i + 1] = [var_name, dp, ad, qual, gt]
            if hap_name in two_cp_haplotypes:
                i += 1
            i += 1
        return vars

    def run(self):
        """Process haplotypes one by one"""
        call_sum = self.call_sum
        vars = self.run_step(
            call_sum["final_haplotypes"],
            call_sum["haplotype_details"],
            self.ref,
            self.offset,
        )
        self.merge_vcf([vars])


class TwoGeneVcfGenerater(VcfGenerater):
    """
    Make vcf for two-gene scenario
    """

    def __init__(self, sample_id, outdir, config, call_sum):
        VcfGenerater.__init__(self, sample_id, outdir, config, call_sum)
        self.nchr_old_gene2 = config["coordinates"]["hg38"]["nchr_old_smn2"]
        self.offset_gene2 = int(self.nchr_old_gene2.split("_")[1]) - 1
        self.ref_gene2 = config["data"]["reference_smn2"]
        self.position_match = config["data"]["smn_match"]
        with open(self.position_match) as f:
            for line in f:
                at = line.split()
                self.match.setdefault(int(at[1]), int(at[3]))

    def run(self):
        """
        Process haplotypes one by one. Realign to different ref sequence
        in this two-gene scenario
        """
        call_sum = self.call_sum
        vars_smn1 = self.run_step(
            call_sum["smn1_haplotypes"],
            call_sum["haplotype_details"],
            self.ref,
            self.offset,
        )
        vars_smn2 = self.run_step(
            call_sum["smn2_haplotypes"],
            call_sum["haplotype_details"],
            self.ref_gene2,
            self.offset_gene2,
            match_range=True,
        )
        self.merge_vcf([vars_smn2, vars_smn1])
