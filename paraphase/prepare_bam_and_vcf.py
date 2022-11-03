# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


import os
import pysam
import subprocess
import random
import gzip
from .haplotype_assembler import VariantGraph


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
        self.use_supplementary = False
        if "use_supplementary" in config:
            self.use_supplementary = config["use_supplementary"]
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
                # and read.get_tag("NM") < read.reference_length * 0.1
                and (self.use_supplementary is True or read.is_supplementary is False)
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
    read_color_allele1 = "178,223,138"
    read_color_allele2 = "177,156,217"

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
        self.use_supplementary = False
        if "use_supplementary" in config:
            self.use_supplementary = config["use_supplementary"]

    def add_tag_to_read(
        self,
        read,
        hp_keys,
        reads_to_tag,
        nonunique,
        read_details,
        alleles=[],
        random_assign=False,
    ):
        """
        Add HP tag to each read.
        Option to randomly assign a read to a possible haplotype when
        nonuniquelly matching. The randomly assigned reads are in gray
        while unique reads are in blue.
        """
        hp_found = False
        read_name = read.qname
        if read.is_supplementary and self.use_supplementary:
            read_name = (
                read_name + f"_sup_{read.reference_start}_{read.reference_length}"
            )
        for hap, hap_name in hp_keys.items():
            if read_name in reads_to_tag[hap]:
                read.set_tag("HP", hap_name, "Z")
                if random_assign:
                    read.set_tag("YC", self.read_color, "Z")
                if alleles != []:
                    if hap in alleles[0]:
                        read.set_tag("YC", self.read_color_allele1, "Z")
                    elif len(alleles) > 1 and hap in alleles[1]:
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
                    keys.append(ass_hap)
                    mismatches.append(mismatch)
                mismatches_sorted = sorted(mismatches)
                if (
                    len(mismatches_sorted) > 1
                    and 0 < mismatches_sorted[0] <= 2
                    and mismatches_sorted[1] >= mismatches_sorted[0] + 2
                ):
                    best_match = keys[mismatches.index(mismatches_sorted[0])]
                    read.set_tag("HP", hp_keys[best_match], "Z")
                    hp_found = True
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
        read_details = call_sum.get("read_details")
        alleles = call_sum.get("alleles")
        if alleles is None:
            alleles = []

        for read in self._bamh.fetch(self.nchr):
            if read.is_secondary is False and (
                self.use_supplementary or read.is_supplementary is False
            ):
                read = self.add_tag_to_read(
                    read,
                    hp_keys,
                    unique_reads,
                    nonunique_reads,
                    read_details,
                    alleles,
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
    Run DeepVariant for each haplotype and generate individual/merged vcfs
    """

    search_range = 200

    def __init__(self, sample_id, outdir, config, call_sum):
        self.sample_id = sample_id
        self.outdir = outdir
        self.call_sum = call_sum
        self.bam = os.path.join(outdir, self.sample_id + "_realigned_tagged.bam")
        self.nchr = config["coordinates"]["hg38"]["nchr"]
        self.nchr_old = config["coordinates"]["hg38"]["nchr_old"]
        nstart = int(self.nchr_old.split("_")[1])
        nend = int(self.nchr_old.split("_")[2])
        self.offset = nstart - 1
        self.vcf_region = f"1-{nend - nstart}"
        if "vcf_region" in config["coordinates"]["hg38"]:
            self.vcf_region = config["coordinates"]["hg38"]["vcf_region"]
        self.nchr_length = config["coordinates"]["hg38"]["nchr_length"]
        self.ref = config["data"]["reference"]
        self.samtools = config["tools"]["samtools"]
        self.minimap2 = config["tools"]["minimap2"]
        self.singularity = config["tools"]["singularity"]
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

    def filter_gt12(self, at):
        """
        Filter variants in the scenario of GT 1/2
        """
        alt = at[4].split(",")
        info = at[-1].split(":")
        vaf = info[4].split(",")
        vaf_float = [float(a) for a in vaf]
        variant_picked_index = vaf_float.index(max(vaf_float))
        new_alt = alt[variant_picked_index]
        new_info = [
            "1/1",
            info[1],
            info[2],
            ",".join(["0", info[3].split(",")[variant_picked_index + 1]]),
            vaf[variant_picked_index],
            ",".join(info[-1].split(",")[:2] + ["0"]),
        ]
        at[4] = new_alt
        at[-1] = ":".join(new_info)
        return at

    def filter_vcf(self, hap_vcf, hap_vcf_out, hap_range, offset):
        """
        Filter vcf and return to standard coordinate
        """
        variants = []
        with gzip.open(hap_vcf) as f, open(hap_vcf_out, "w") as fout:
            for line in f:
                line = line.decode("utf8")
                if line[0] == "#":
                    if line.startswith("##contig="):
                        fout.write(
                            f"##contig=<ID={self.nchr},length={self.nchr_length}>\n"
                        )
                    else:
                        fout.write(line)
                else:
                    at = line.split()
                    gt = at[-1].split(":")[0]
                    pos = int(at[1]) + offset
                    if (
                        at[6] == "PASS"
                        and hap_range[0] < pos < hap_range[1]
                        and "0" not in gt
                    ):
                        if gt != "1/1":
                            at = self.filter_gt12(at)
                        at[0] = self.nchr
                        at[1] = str(pos)
                        fout.write("\t".join(at) + "\n")

                        ad = at[-1].split(":")[3].split(",")[1]
                        qual = float(at[5])
                        var_name = f"{pos}_{at[3]}_{at[4]}"
                        variants.append([var_name, ad, qual])
        return variants

    def merge_vcf(self, vars_ranges):
        """
        Merge vcfs from multiple haplotypes. Very rough merging.
        """
        merged_vcf = os.path.join(self.outdir, self.sample_id + f"_variants.vcf")
        with open(merged_vcf, "w") as fout:
            fout.write("##fileformat=VCFv4.2\n")
            fout.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
            fout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
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
            for vars, ranges in vars_ranges:
                vars = dict(sorted(vars.items()))
                for var in vars:
                    pos, ref, alt = var.split("_")
                    pos = int(pos)
                    call_info = vars[var]
                    merge_gt = []
                    merge_ad = []
                    merge_qual = []
                    for i, each_call in enumerate(call_info):
                        hap_range = ranges[i]
                        if each_call is None:
                            if hap_range[0] < pos < hap_range[1] or None in hap_range:
                                merge_gt.append("0")
                                merge_ad.append("0")
                            else:
                                merge_gt.append(".")
                                merge_ad.append(".")
                        else:
                            merge_gt.append("1")
                            merge_ad.append(str(each_call[0]))
                            merge_qual.append(each_call[1])

                    merged_entry = [
                        self.nchr,
                        str(pos),
                        ".",
                        ref,
                        alt,
                        str(min(merge_qual)),
                        "PASS",
                        ".",
                        "GT:AD",
                        "/".join(merge_gt) + ":" + ",".join(merge_ad),
                    ]
                    fout.write("\t".join(merged_entry) + "\n")

    def run_step(
        self,
        final_haps,
        haplotype_details,
        ref_seq,
        ref_name,
        offset,
        match_range=False,
    ):
        """
        Process one haplotype
        """
        vars = {}
        ranges = []
        two_cp_haplotypes = self.call_sum.get("two_copy_haplotypes")
        nhap = len(final_haps)
        for i, hap_name in enumerate(final_haps.values()):
            hap_range = haplotype_details[hap_name]["boundary"]
            if match_range:
                hap_range = [
                    self.get_range_in_other_gene(hap_range[0]),
                    self.get_range_in_other_gene(hap_range[1]),
                ]
            ranges.append(hap_range)
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

            hap_vcf = os.path.join(
                self.outdir, self.sample_id + f"_{hap_name}_prefilter.vcf.gz"
            )
            vcf_cmd = (
                f"{self.singularity} exec --bind /usr/lib/locale/ docker://google/deepvariant:1.3.0 "
                + f"/opt/deepvariant/bin/run_deepvariant --model_type PACBIO --ref {ref_seq}  --reads {hap_bam} "
                + f"--output_vcf {hap_vcf} --num_shards 1 --regions {ref_name}:{self.vcf_region}"
            )
            result = subprocess.run(vcf_cmd, capture_output=True, text=True, shell=True)
            result.check_returncode()

            # fileter vcf
            hap_vcf_out = os.path.join(self.outdir, self.sample_id + f"_{hap_name}.vcf")
            variants_called = self.filter_vcf(hap_vcf, hap_vcf_out, hap_range, offset)

            for var_name, ad, qual in variants_called:
                vars.setdefault(var_name, [None] * nhap)
                vars[var_name][i] = [ad, qual]
        # insert info for two-copy haplotype
        if two_cp_haplotypes is not None and two_cp_haplotypes != []:
            for i, hap_name in enumerate(final_haps.values()):
                if hap_name in two_cp_haplotypes:
                    for var in vars:
                        var_call = vars[var][i]
                        vars[var].insert(i, var_call)
                    range_call = ranges[i]
                    ranges.insert(i, range_call)
        return vars, ranges

    def run(self):
        """Process haplotypes one by one"""
        call_sum = self.call_sum
        vars, ranges = self.run_step(
            call_sum["final_haplotypes"],
            call_sum["haplotype_details"],
            self.ref,
            self.nchr_old,
            self.offset,
        )
        self.merge_vcf([(vars, ranges)])


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
        vars_smn1, ranges_smn1 = self.run_step(
            call_sum["smn1_haplotypes"],
            call_sum["haplotype_details"],
            self.ref,
            self.nchr_old,
            self.offset,
        )
        vars_smn2, ranges_smn2 = self.run_step(
            call_sum["smn2_haplotypes"],
            call_sum["haplotype_details"],
            self.ref_gene2,
            self.nchr_old_gene2,
            self.offset_gene2,
            match_range=True,
        )
        self.merge_vcf([(vars_smn2, ranges_smn2), (vars_smn1, ranges_smn1)])
