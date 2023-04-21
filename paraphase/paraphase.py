# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


import os
import argparse
import json
import yaml
import logging
import datetime
import shutil
import pysam
import subprocess
import re
from glob import glob
import multiprocessing as mp
from argparse import RawTextHelpFormatter
from functools import partial
from paraphase.genome_depth import GenomeDepth
from paraphase.prepare_bam_and_vcf import (
    BamRealigner,
    BamTagger,
    VcfGenerater,
    TwoGeneVcfGenerater,
)
from paraphase import genes
from .phaser import Phaser


class Paraphase:
    def __init__(self):
        data_path = os.path.join(os.path.dirname(__file__), "data")
        gene_config = os.path.join(data_path, "genes.yaml")
        self.accepted_genes = []
        for config_file in glob(os.path.join(data_path, "*/*.yaml")):
            gene_name = config_file.split("/")[-2]
            self.accepted_genes.append(gene_name)
        with open(gene_config, "r") as f:
            try:
                gene_config_parsed = yaml.safe_load(f)
            except yaml.YAMLError as yaml_error:
                raise Exception(f"Error reading {gene_config}: \n\n{yaml_error}")
        if gene_config_parsed.get("accepted_genes") is not None:
            self.accepted_genes = gene_config_parsed["accepted_genes"]
        self.genome_depth_genes = gene_config_parsed.get("genome_depth_genes")
        self.no_vcf_genes = gene_config_parsed.get("no_vcf_genes")
        self.check_sex_genes = gene_config_parsed.get("check_sex_genes")
        self.samtools = None
        self.minimap2 = None

    def process_sample(self, bamlist, outdir, configs, dcov={}, novcf=False):
        """Main workflow"""
        for bam in bamlist:
            sample_id = bam.split("/")[-1].split(".")[0]
            logging.info(
                f"Processing sample {sample_id} at {datetime.datetime.now()}..."
            )
            sample_out = {}
            gdepth = None
            sample_sex = None
            query_genes = set(configs.keys())
            if query_genes.intersection(set(self.genome_depth_genes)) != set():
                logging.info(f"Getting genome depth at {datetime.datetime.now()}...")
                if sample_id in dcov:
                    gdepth = dcov[sample_id]
                if gdepth is None:
                    depth = GenomeDepth(
                        bam,
                        os.path.join(
                            os.path.dirname(__file__), "data", "genome_region.bed"
                        ),
                    )
                    gdepth, gmad = depth.call()
                    if gdepth < 10 or gmad > 0.25:
                        logging.warning(
                            "Due to low or highly variable genome coverage, genome coverage is not used for depth correction."
                        )
                        gdepth = None
            if query_genes.intersection(set(self.check_sex_genes)) != set():
                logging.info(f"Checking sample sex at {datetime.datetime.now()}...")
                depth = GenomeDepth(
                    bam,
                    os.path.join(os.path.dirname(__file__), "data", "sex_region.bed"),
                )
                sample_sex = depth.check_sex()

            for gene in configs:
                config = configs[gene]
                logging.info(
                    f"Running analysis for {gene} at {datetime.datetime.now()}..."
                )
                logging.info(
                    f"Realigning reads for {gene} at {datetime.datetime.now()}..."
                )
                bam_realigner = BamRealigner(bam, outdir, config)
                bam_realigner.write_realign_bam()

                logging.info(
                    f"Phasing haplotypes for {gene} at {datetime.datetime.now()}..."
                )

                phasers = {
                    "smn1": genes.Smn1Phaser(sample_id, outdir, gdepth),
                    "rccx": genes.RccxPhaser(sample_id, outdir),
                    "pms2": genes.Pms2Phaser(sample_id, outdir),
                    "strc": genes.StrcPhaser(sample_id, outdir, gdepth, genome_bam=bam),
                    "ncf1": genes.Ncf1Phaser(sample_id, outdir, gdepth),
                    "cfc1": genes.Cfc1Phaser(sample_id, outdir),
                    "neb": genes.NebPhaser(sample_id, outdir),
                    "ikbkg": genes.IkbkgPhaser(sample_id, outdir),
                    "f8": genes.F8Phaser(sample_id, outdir, genome_bam=bam),
                    "opn1lw": genes.Opn1lwPhaser(
                        sample_id, outdir, sample_sex=sample_sex
                    ),
                }
                phaser = phasers.get(gene)
                # not pre-included gene
                if phaser is None:
                    phaser = Phaser(sample_id, outdir)

                phaser.set_parameter(config)

                phaser_call = phaser.call()
                if phaser_call is not None:
                    phaser_call = phaser_call._asdict()

                    logging.info(
                        f"Tagging reads for {gene} at {datetime.datetime.now()}..."
                    )
                    bam_tagger = BamTagger(sample_id, outdir, config, phaser_call)
                    bam_tagger.write_bam(random_assign=True)

                    if novcf is False and gene not in self.no_vcf_genes:
                        logging.info(
                            f"Generating VCFs for {gene} at {datetime.datetime.now()}..."
                        )
                        vcf_dir = os.path.join(outdir, f"{sample_id}_vcfs")
                        if os.path.exists(vcf_dir) is False:
                            os.makedirs(vcf_dir)
                        if gene == "smn1":
                            vcf_generater = TwoGeneVcfGenerater(
                                sample_id, outdir, config, phaser_call
                            )
                            vcf_generater.run()
                        else:
                            vcf_generater = VcfGenerater(
                                sample_id, outdir, config, phaser_call
                            )
                            vcf_generater.run_without_realign()

                sample_out.setdefault(gene, phaser_call)

            logging.info(f"Merging all bams at {datetime.datetime.now()}...")
            self.merge_bams(query_genes, outdir, sample_id)

            logging.info(f"Writing to json at {datetime.datetime.now()}...")
            out_json = os.path.join(outdir, sample_id + ".json")
            with open(out_json, "w") as json_output:
                json.dump(sample_out, json_output, indent=4)

    @staticmethod
    def merge_bams(query_genes, outdir, sample_id):
        """Merge realigned tagged bams for each gene into one bam"""
        bams = []
        for gene in query_genes:
            gene_bam = os.path.join(outdir, f"{sample_id}_{gene}_realigned_tagged.bam")
            if os.path.exists(gene_bam) is False:
                gene_bam = os.path.join(outdir, f"{sample_id}_{gene}_realigned.bam")
            if os.path.exists(gene_bam):
                bams.append(gene_bam)
        bam_list_file = os.path.join(outdir, f"{sample_id}_bam_list.txt")
        with open(bam_list_file, "w") as fout:
            for bam in bams:
                fout.write(bam + "\n")
        merged_bam = os.path.join(outdir, f"{sample_id}_realigned_tagged.bam")
        tmp_bam = os.path.join(outdir, f"{sample_id}_merged.bam")
        pysam.merge("-f", "-o", tmp_bam, "-b", bam_list_file)
        pysam.sort("-o", merged_bam, tmp_bam)
        pysam.index(merged_bam)
        os.remove(tmp_bam)
        os.remove(bam_list_file)
        for bam in bams:
            os.remove(bam)
            os.remove(bam + ".bai")

    def get_samtools_minimap2_path(self, args):
        # check third-party tools
        samtools_check = [
            a for a in [args.samtools, shutil.which("samtools")] if a is not None
        ]
        samtools_check2 = [a for a in samtools_check if os.path.exists(a)]
        if samtools_check2 == []:
            raise Exception("samtools is not found")
        else:
            self.samtools = samtools_check2[0]

        minimap2_check = [
            a for a in [args.minimap2, shutil.which("minimap2")] if a is not None
        ]
        minimap2_check2 = [a for a in minimap2_check if os.path.exists(a)]
        if minimap2_check2 == []:
            raise Exception("minimap2 is not found")
        else:
            self.minimap2 = minimap2_check2[0]

    def update_config(self, gene_list, args):
        """Get config info for each gene"""
        data_path = os.path.join(os.path.dirname(__file__), "data")
        configs = {}
        for gene in gene_list:
            config_file = os.path.join(data_path, gene, f"{gene}_config.yaml")
            # parse config file
            with open(config_file, "r") as f:
                try:
                    config = yaml.safe_load(f)
                except yaml.YAMLError as yaml_error:
                    raise Exception(f"Error reading {config_file}: \n\n{yaml_error}")
            configs.setdefault(gene, config)

            # check third-party tools
            configs[gene].setdefault("tools", {})
            configs[gene]["tools"].setdefault("samtools", self.samtools)
            configs[gene]["tools"].setdefault("minimap2", self.minimap2)

            # update paths
            data_paths = configs[gene].get("data")
            for data_entry in data_paths:
                old_data_file = data_paths[data_entry]
                if data_entry != "depth_region":
                    new_data_file = os.path.join(data_path, gene, old_data_file)
                else:
                    new_data_file = os.path.join(data_path, old_data_file)
                data_paths[data_entry] = new_data_file

            for data_file in list(data_paths.values()):
                if os.path.exists(data_file) is False:
                    raise Exception(f"File {data_file} not found.")
        return configs

    def load_parameters(self):
        parser = argparse.ArgumentParser(
            description="paraphase: HiFi-based caller for highly homologous genes",
            formatter_class=RawTextHelpFormatter,
        )
        all_genes_joined = ",".join(self.accepted_genes)
        inputp = parser.add_argument_group("Input Options")
        outputp = parser.add_argument_group("Output Options")
        newgenep = parser.add_argument_group(
            "Options to provide a user-specified region of interest"
        )
        input_group = inputp.add_mutually_exclusive_group(required=True)
        input_group.add_argument(
            "-b",
            "--bam",
            help="Input BAM file, mutually exclusive with -l",
        )
        input_group.add_argument(
            "-l",
            "--list",
            help="File listing absolute paths to multiple input BAM files, one per line",
        )
        outputp.add_argument(
            "-o",
            "--out",
            help="Output directory",
            required=True,
        )
        parser.add_argument(
            "-g",
            "--gene",
            help="Optionally specify which gene(s) to run (separated by comma).\n"
            + "Will run all genes if not specified.\n"
            + f"Currently supported genes are {all_genes_joined}.\n"
            + "Genes to be analyzed can be defined in paraphase/data/genes.yaml",
            required=False,
        )
        parser.add_argument(
            "-d",
            "--depth",
            help="Optional path to a file listing average depth for each sample",
            required=False,
        )
        parser.add_argument(
            "--novcf",
            help="Optional. If specified, paraphase will not write vcfs",
            required=False,
            action="store_true",
        )
        parser.add_argument(
            "-t",
            "--threads",
            help="Optional number of threads. Only used together with -l",
            required=False,
            type=int,
            default=1,
        )
        parser.add_argument(
            "--samtools",
            help="Optional path to samtools",
            required=False,
        )
        parser.add_argument(
            "--minimap2",
            help="Optional path to minimap2",
            required=False,
        )
        # ask Paraphase to analyze a new region
        newgenep.add_argument(
            "-n",
            "--newgene",
            help=f"New gene name",
            required=False,
        )
        newgenep.add_argument(
            "-r1",
            "--region1",
            help=f"Region of interest, format chr1:1234567-1234589",
            required=False,
        )
        newgenep.add_argument(
            "-r2",
            "--region2",
            help=f"Additional region to extract reads, format chr1:1234567-1234589",
            required=False,
        )
        newgenep.add_argument(
            "--genome",
            help=f"Path to reference fasta file. Required if providing user specified region of interest",
            required=False,
        )
        return parser

    def prepare_data_for_new_gene(self, gene, region1, region2, genome):
        """make data dir and files for new region"""
        gene_data_dir = os.path.join(os.path.dirname(__file__), "data", gene)
        if os.path.exists(gene_data_dir) is False:
            os.makedirs(gene_data_dir)
        nchr = region1.split(":")[0]
        gene_data_config = os.path.join(gene_data_dir, f"{gene}_config.yaml")
        # make ref
        genome_ref = pysam.FastaFile(genome)
        ref_file = os.path.join(gene_data_dir, f"{gene}_ref.fa")
        make_ref_cmd = (
            f"{self.samtools} faidx {genome} {region1} | sed -e "
            + f'"s/-/_/"'
            + " | sed -e "
            + f'"s/:/_/"'
            + f" > {ref_file}"
        )
        result = subprocess.run(
            make_ref_cmd, capture_output=True, text=True, shell=True
        )
        result.check_returncode()
        pysam.faidx(ref_file)

        chr_length = genome_ref.get_reference_length(nchr)
        genome_ref.close()
        # make config file
        new_ref_name = region1.replace(":", "_").replace("-", "_")
        with open(gene_data_config, "w") as fout:
            fout.write(f"gene: {gene}\n")
            fout.write("data:\n")
            fout.write(f" reference: {gene}_ref.fa\n")
            fout.write("coordinates:\n")
            fout.write(" hg38:\n")
            fout.write(f"  nchr: {nchr}\n")
            fout.write(f"  nchr_old: {new_ref_name}\n")
            fout.write(f"  nchr_length: {chr_length}\n")
            fout.write(f"  extract_region1: {region1}\n")
            fout.write(f"  extract_region2: {region2}\n")

    def run(self):
        parser = self.load_parameters()
        args = parser.parse_args()
        outdir = args.out
        logging.basicConfig(level=logging.DEBUG)
        self.get_samtools_minimap2_path(args)
        if os.path.exists(outdir) is False:
            os.makedirs(outdir)

        gene_list = []
        if None not in [args.newgene, args.region1, args.region2, args.genome]:
            # check region format
            # if (
            #    re.findall(r"chr\d+:\d+-\d+", args.region1) == []
            #    or re.findall(r"chr\d+:\d+-\d+", args.region2) == []
            # ):
            #    raise Exception(
            #        "Region should be given in the following format chr1:1234567-1234589"
            #    )
            self.prepare_data_for_new_gene(
                args.newgene, args.region1, args.region2, args.genome
            )
            gene_list = [args.newgene]
        if gene_list == []:
            all_genes_joined = ",".join(self.accepted_genes)
            gene_list = args.gene
            if gene_list is None:
                gene_list = self.accepted_genes
            else:
                gene_list = [
                    a for a in gene_list.split(",") if a in self.accepted_genes
                ]
                if gene_list == []:
                    raise Exception(
                        f"Gene names not recognized. Currently accepted genes are {all_genes_joined}"
                    )

        configs = self.update_config(gene_list, args)

        # parse depth file
        dcov = {}
        if args.depth is not None:
            with open(args.depth) as f:
                for line in f:
                    split_line = line.split()
                    dcov.setdefault(split_line[0], float(split_line[-1]))

        # process sample(s)
        bamlist = []
        if args.bam is not None:
            if os.path.exists(args.bam) and os.path.exists(args.bam + ".bai"):
                bamlist = [args.bam]
                self.process_sample(
                    bamlist,
                    outdir,
                    configs,
                    dcov,
                    args.novcf,
                )
            else:
                print(f"{args.bam} bam or bai file doesn't exist")
        elif args.list is not None:
            with open(args.list) as f:
                for line in f:
                    bam = line[:-1]
                    if os.path.exists(bam) and os.path.exists(bam + ".bai"):
                        bamlist.append(bam)
                    else:
                        print(f"{bam} bam or bai file doesn't exist")

            nCores = args.threads
            process_sample_partial = partial(
                self.process_sample,
                outdir=outdir,
                configs=configs,
                dcov=dcov,
                novcf=args.novcf,
            )
            bam_groups = [bamlist[i::nCores] for i in range(nCores)]
            pool = mp.Pool(nCores)
            pool.map(process_sample_partial, bam_groups)
            pool.close()
            pool.join()
        else:
            raise Exception("No input file is given")
