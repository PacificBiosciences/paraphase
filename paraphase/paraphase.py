# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


import os
import sys
import argparse
import json
import yaml
import logging
import datetime
import shutil
import pysam
import subprocess
import traceback
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
import paraphase
from paraphase import genes
from .phaser import Phaser

logging.basicConfig(level=logging.DEBUG)


class Paraphase:
    def __init__(self):
        self.samtools = None
        self.minimap2 = None

    def parse_configs(self, region_config=None):
        """
        Parse config files
        region_config: parameters for each region of interest
        gene_config: specifies some additional analyses triggered with some genes
        """
        data_path = os.path.join(os.path.dirname(__file__), "data")
        gene_config = os.path.join(data_path, "genes.yaml")
        if region_config is None:
            region_config = os.path.join(data_path, "config.yaml")
        self.region_config_parsed = None
        with open(region_config, "r") as f:
            self.region_config_parsed = yaml.safe_load(f)
        self.accepted_genes = list(self.region_config_parsed.keys())
        with open(gene_config, "r") as f:
            gene_config_parsed = yaml.safe_load(f)
        self.genes_to_call = gene_config_parsed.get("genes_to_call")
        self.genome_depth_genes = gene_config_parsed.get("genome_depth_genes")
        self.no_vcf_genes = gene_config_parsed.get("no_vcf_genes")
        self.check_sex_genes = gene_config_parsed.get("check_sex_genes")

    def process_gene(
        self,
        gene_list,
        configs,
        outdir,
        tmpdir,
        gdepth,
        bam,
        sample_sex,
        novcf,
        prog_cmd,
    ):
        """Workflow for each region"""
        phaser_calls = {}
        for gene in gene_list:
            try:
                sample_id = bam.split("/")[-1].split(".")[0]
                if gene == "smn1":
                    phaser = genes.Smn1Phaser(
                        sample_id, tmpdir, gdepth, bam, sample_sex
                    )
                elif gene == "rccx":
                    phaser = genes.RccxPhaser(
                        sample_id, tmpdir, gdepth, bam, sample_sex
                    )
                elif gene == "pms2":
                    phaser = genes.Pms2Phaser(
                        sample_id, tmpdir, gdepth, bam, sample_sex
                    )
                elif gene == "strc":
                    phaser = genes.StrcPhaser(
                        sample_id, tmpdir, gdepth, bam, sample_sex
                    )
                elif gene == "ncf1":
                    phaser = genes.Ncf1Phaser(
                        sample_id, tmpdir, gdepth, bam, sample_sex
                    )
                elif gene == "cfc1":
                    phaser = genes.Cfc1Phaser(
                        sample_id, tmpdir, gdepth, bam, sample_sex
                    )
                elif gene == "neb":
                    phaser = genes.NebPhaser(sample_id, tmpdir, gdepth, bam, sample_sex)
                elif gene == "ikbkg":
                    phaser = genes.IkbkgPhaser(
                        sample_id, tmpdir, gdepth, bam, sample_sex
                    )
                elif gene == "f8":
                    phaser = genes.F8Phaser(sample_id, tmpdir, gdepth, bam, sample_sex)
                else:
                    phaser = Phaser(sample_id, tmpdir)

                config = configs[gene]
                logging.info(
                    f"Running analysis for {gene} for sample {sample_id} at {datetime.datetime.now()}..."
                )
                logging.info(
                    f"Realigning reads for {gene} for sample {sample_id} at {datetime.datetime.now()}..."
                )
                bam_realigner = BamRealigner(bam, tmpdir, config, prog_cmd)
                bam_realigner.write_realign_bam()

                logging.info(
                    f"Phasing haplotypes for {gene} for sample {sample_id} at {datetime.datetime.now()}..."
                )

                phaser.set_parameter(config)
                phaser_call = phaser.call()._asdict()
                phaser_calls.setdefault(gene, phaser_call)

                logging.info(
                    f"Tagging reads for {gene} for sample {sample_id} at {datetime.datetime.now()}..."
                )
                bam_tagger = BamTagger(sample_id, tmpdir, config, phaser_call)
                bam_tagger.write_bam(random_assign=True)

                if novcf is False and gene not in self.no_vcf_genes:
                    logging.info(
                        f"Generating VCFs for {gene} for sample {sample_id} at {datetime.datetime.now()}..."
                    )
                    if gene == "smn1":
                        vcf_generater = TwoGeneVcfGenerater(
                            sample_id,
                            outdir,
                            phaser_call,
                        )
                        vcf_generater.set_parameter(
                            config, tmpdir=tmpdir, prog_cmd=prog_cmd
                        )
                        vcf_generater.run()
                    else:
                        vcf_generater = VcfGenerater(
                            sample_id,
                            outdir,
                            phaser_call,
                        )
                        vcf_generater.set_parameter(
                            config, tmpdir=tmpdir, prog_cmd=prog_cmd
                        )
                        vcf_generater.run_without_realign()
            except Exception:
                logging.error(
                    f"Error running {gene} for sample {sample_id}...See error message below"
                )
                traceback.print_exc()
                phaser_calls.setdefault(gene, None)
        return phaser_calls

    def process_sample(
        self,
        bamlist,
        outdir,
        configs,
        tmpdir,
        prog_cmd,
        num_threads=1,
        dcov={},
        novcf=False,
    ):
        """Main workflow"""
        for bam in bamlist:
            try:
                sample_id = bam.split("/")[-1].split(".")[0]
                logging.info(
                    f"Processing sample {sample_id} at {datetime.datetime.now()}..."
                )
                sample_out = {}
                gdepth = None
                sample_sex = None
                query_genes = list(configs.keys())
                if set(query_genes).intersection(set(self.genome_depth_genes)) != set():
                    logging.info(
                        f"Getting genome depth for sample {sample_id} at {datetime.datetime.now()}..."
                    )
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
                                f"For sample {sample_id}, due to low or highly variable genome coverage, genome coverage is not used for depth correction."
                            )
                            gdepth = None

                if num_threads == 1:
                    sample_out = self.process_gene(
                        query_genes,
                        configs,
                        outdir,
                        tmpdir,
                        gdepth,
                        bam,
                        sample_sex,
                        novcf,
                        prog_cmd,
                    )
                else:
                    process_gene_partial = partial(
                        self.process_gene,
                        configs=configs,
                        outdir=outdir,
                        tmpdir=tmpdir,
                        gdepth=gdepth,
                        bam=bam,
                        sample_sex=sample_sex,
                        novcf=novcf,
                        prog_cmd=prog_cmd,
                    )
                    gene_groups = [
                        query_genes[i::num_threads] for i in range(num_threads)
                    ]
                    pool = mp.Pool(num_threads)
                    phaser_calls = pool.map(process_gene_partial, gene_groups)
                    pool.close()
                    pool.join()
                    for phaser_call_set in phaser_calls:
                        sample_out.update(phaser_call_set)
                sample_out = dict(sorted(sample_out.items()))

                logging.info(
                    f"Merging all bams for sample {sample_id} at {datetime.datetime.now()}..."
                )
                self.merge_bams(query_genes, outdir, tmpdir, sample_id)

                logging.info(
                    f"Writing to json for sample {sample_id} at {datetime.datetime.now()}..."
                )
                out_json = os.path.join(outdir, sample_id + ".json")
                with open(out_json, "w") as json_output:
                    json.dump(sample_out, json_output, indent=4)
            except Exception:
                logging.error(
                    f"Error running sample {sample_id}...See error message below"
                )
                traceback.print_exc()

    @staticmethod
    def merge_bams(query_genes, outdir, tmpdir, sample_id):
        """Merge realigned tagged bams for each gene into one bam"""
        bams = []
        for gene in query_genes:
            gene_bam = os.path.join(tmpdir, f"{sample_id}_{gene}_realigned_tagged.bam")
            if os.path.exists(gene_bam) is False:
                gene_bam = os.path.join(tmpdir, f"{sample_id}_{gene}_realigned.bam")
            if os.path.exists(gene_bam):
                bams.append(gene_bam)
        bam_list_file = os.path.join(tmpdir, f"{sample_id}_bam_list.txt")
        with open(bam_list_file, "w") as fout:
            for bam in bams:
                fout.write(bam + "\n")
        merged_bam = os.path.join(outdir, f"{sample_id}_realigned_tagged.bam")
        tmp_bam = os.path.join(tmpdir, f"{sample_id}_merged.bam")
        pysam.merge("-f", "-o", tmp_bam, "-b", bam_list_file)
        pysam.sort("-o", merged_bam, tmp_bam)
        pysam.index(merged_bam)

    def get_samtools_minimap2_path(self, args):
        """Get and check path to third-party tools"""
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

    def update_config(self, gene_list, ref_dir, genome):
        """Get config info for each gene"""
        genome_ref = pysam.FastaFile(genome)
        data_path = os.path.join(os.path.dirname(__file__), "data")
        configs = {}
        for gene in gene_list:
            assert gene in self.region_config_parsed
            configs.setdefault(gene, self.region_config_parsed[gene])
            configs[gene].setdefault("gene", gene)
            # check third-party tools
            configs[gene].setdefault("tools", {})
            configs[gene]["tools"].setdefault("samtools", self.samtools)
            configs[gene]["tools"].setdefault("minimap2", self.minimap2)
            # update paths
            configs[gene].setdefault("data", {})
            data_paths = configs[gene]["data"]
            for data_entry in data_paths:
                old_data_file = data_paths[data_entry]
                new_data_file = os.path.join(data_path, gene, old_data_file)
                data_paths[data_entry] = new_data_file
            # add reference fasta
            ref_file = os.path.join(ref_dir, f"{gene}_ref.fa")
            data_paths.setdefault("reference", ref_file)
            realign_region = configs[gene].get("realign_region")
            self.make_ref_fasta(ref_file, realign_region, genome)
            # check files exist
            for data_file in list(data_paths.values()):
                if os.path.exists(data_file) is False:
                    raise Exception(f"File {data_file} not found.")
            nchr = realign_region.split(":")[0]
            nchr_length = genome_ref.get_reference_length(nchr)
            nchr_old = realign_region.replace(":", "_").replace("-", "_")
            configs[gene].setdefault("nchr", nchr)
            configs[gene].setdefault("nchr_length", nchr_length)
            configs[gene].setdefault("nchr_old", nchr_old)
        genome_ref.close()
        return configs

    def make_ref_fasta(self, ref_file, realign_region, genome):
        """Create a reference file for each region of interest for realignment"""
        make_ref_cmd = (
            f"{self.samtools} faidx {genome} {realign_region} | sed -e "
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

    def get_gene_list(self, gene_input):
        """Get a list of genes to analyze"""
        if gene_input is None:
            if self.genes_to_call is not None and self.genes_to_call != []:
                return self.genes_to_call
        else:
            gene_list = set()
            for gene in gene_input.split(","):
                if gene in self.accepted_genes:
                    gene_list.add(gene)
                else:
                    logging.warning(
                        f"{gene} is not a valid gene name and will be skipped..."
                    )
            return list(gene_list)
        return self.accepted_genes

    def load_parameters(self):
        parser = argparse.ArgumentParser(
            description="Paraphase: HiFi-based caller for highly homologous genes",
            formatter_class=RawTextHelpFormatter,
        )
        inputp = parser.add_argument_group("Input Options")
        outputp = parser.add_argument_group("Output Options")
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
        inputp.add_argument(
            "-r",
            "--reference",
            help=f"Path to reference genome fasta file",
            required=True,
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
            + "The full set of accepted genes are defined in the config file.\n",
            required=False,
        )
        parser.add_argument(
            "-c",
            "--config",
            help="Optional path to a user-defined config file listing the full set of regions to analyze.\n"
            + "By default Paraphase uses the config file in paraphase/data/config.yaml",
            required=False,
        )
        parser.add_argument(
            "-d",
            "--depth",
            help="Optional path to a file listing average depth for each sample",
            required=False,
        )
        parser.add_argument(
            "-t",
            "--threads",
            help="Optional number of threads",
            required=False,
            type=int,
            default=1,
        )
        parser.add_argument(
            "--novcf",
            help="Optional. If specified, paraphase will not write vcfs",
            required=False,
            action="store_true",
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
        parser.add_argument(
            "-v", "--version", action="version", version=f"{paraphase.__version__}"
        )
        return parser

    def run(self):
        parser = self.load_parameters()
        args = parser.parse_args()
        num_threads = args.threads
        outdir = args.out
        os.makedirs(outdir, exist_ok=True)
        tmpdir = os.path.join(
            outdir, f"tmp_{datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S-%f')}"
        )
        os.makedirs(tmpdir, exist_ok=True)

        try:
            prog_cmd = " ".join(sys.argv[1:])
            self.parse_configs(region_config=args.config)
            self.get_samtools_minimap2_path(args)
            gene_list = self.get_gene_list(args.gene)
            if gene_list == []:
                all_genes_joined = ",".join(self.accepted_genes)
                raise Exception(
                    f"Please provide valid gene name(s). Accepted genes are {all_genes_joined}."
                )
            configs = self.update_config(gene_list, tmpdir, args.reference)

            # parse depth file
            dcov = {}
            if args.depth is not None:
                with open(args.depth) as f:
                    for line in f:
                        split_line = line.split()
                        dcov.setdefault(split_line[0], float(split_line[-1]))

            # process sample(s)
            bamlist = []
            # one bam, multiprocess by gene
            if args.bam is not None:
                if os.path.exists(args.bam) and os.path.exists(args.bam + ".bai"):
                    bamlist = [args.bam]
                    self.process_sample(
                        bamlist,
                        outdir,
                        configs,
                        tmpdir,
                        prog_cmd,
                        num_threads,
                        dcov,
                        args.novcf,
                    )
                else:
                    logging.warning(f"{args.bam} bam or bai file doesn't exist")
            # multiple bams, multiprocess by sample
            elif args.list is not None:
                with open(args.list) as f:
                    for line in f:
                        bam = line[:-1]
                        if os.path.exists(bam) and os.path.exists(bam + ".bai"):
                            bamlist.append(bam)
                        else:
                            logging.warning(f"{bam} bam or bai file doesn't exist")

                process_sample_partial = partial(
                    self.process_sample,
                    outdir=outdir,
                    configs=configs,
                    tmpdir=tmpdir,
                    prog_cmd=prog_cmd,
                    dcov=dcov,
                    novcf=args.novcf,
                )
                bam_groups = [bamlist[i::num_threads] for i in range(num_threads)]
                pool = mp.Pool(num_threads)
                pool.map(process_sample_partial, bam_groups)
                pool.close()
                pool.join()
            else:
                raise Exception("No input file is given")
        except Exception:
            logging.error("Error running the program...See error message below")
            traceback.print_exc()
        finally:
            # remove temp dir
            if os.path.exists(tmpdir):
                shutil.rmtree(tmpdir)
            logging.info(
                f"Completed Paraphase analysis at {datetime.datetime.now()}..."
            )
