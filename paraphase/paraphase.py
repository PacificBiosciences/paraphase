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
import multiprocessing as mp
from functools import partial
from paraphase.genome_depth import GenomeDepth
from paraphase.prepare_bam_and_vcf import (
    BamRealigner,
    BamTagger,
    VcfGenerater,
    TwoGeneVcfGenerater,
)
from paraphase.genes.smn1_phaser import Smn1Phaser
from paraphase.genes.pms2_phaser import Pms2Phaser
from paraphase.genes.rccx_phaser import RccxPhaser
from paraphase.genes.strc_phaser import StrcPhaser
from paraphase.genes.ncf1_phaser import Ncf1Phaser
from paraphase.genes.cfc1_phaser import Cfc1Phaser
from paraphase.genes.neb_phaser import NebPhaser
from paraphase.genes.ikbkg_phaser import IkbkgPhaser
from paraphase.genes.f8_phaser import F8Phaser

ACCEPTED_GENES = [
    "smn1",
    "rccx",
    "pms2",
    "strc",
    "ncf1",
    "cfc1",
    "neb",
    "ikbkg",
    "f8",
]
GENOME_DEPTH_GENES = ["smn1", "strc", "ncf1"]
NO_VCF_GENES = ["f8"]


def process_sample(bamlist, outdir, configs, dcov={}):
    """Main workflow"""
    for bam in bamlist:
        sample_id = bam.split("/")[-1].split(".")[0]
        logging.info(f"Processing sample {sample_id} at {datetime.datetime.now()}...")
        sample_out = {}
        gdepth = None
        query_genes = set(configs.keys())
        if query_genes.intersection(set(GENOME_DEPTH_GENES)) != set():
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

        for gene in configs:
            config = configs[gene]
            logging.info(f"Running analysis for {gene} at {datetime.datetime.now()}...")
            logging.info(f"Realigning reads for {gene} at {datetime.datetime.now()}...")
            bam_realigner = BamRealigner(bam, outdir, config)
            bam_realigner.write_realign_bam()

            logging.info(
                f"Phasing haplotypes for {gene} at {datetime.datetime.now()}..."
            )

            phasers = {
                "smn1": Smn1Phaser(sample_id, outdir, gdepth),
                "rccx": RccxPhaser(sample_id, outdir),
                "pms2": Pms2Phaser(sample_id, outdir),
                "strc": StrcPhaser(sample_id, outdir, gdepth, genome_bam=bam),
                "ncf1": Ncf1Phaser(sample_id, outdir, gdepth),
                "cfc1": Cfc1Phaser(sample_id, outdir),
                "neb": NebPhaser(sample_id, outdir),
                "ikbkg": IkbkgPhaser(sample_id, outdir),
                "f8": F8Phaser(sample_id, outdir, genome_bam=bam),
            }
            phaser = phasers.get(gene)
            phaser.set_parameter(config)

            phaser_call = phaser.call()
            if phaser_call is not None:
                phaser_call = phaser_call._asdict()

                logging.info(
                    f"Tagging reads for {gene} at {datetime.datetime.now()}..."
                )
                bam_tagger = BamTagger(sample_id, outdir, config, phaser_call)
                bam_tagger.write_bam(random_assign=True)

                if 0:  # gene not in NO_VCF_GENES:
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
        merge_bams(query_genes, outdir, sample_id)

        logging.info(f"Writing to json at {datetime.datetime.now()}...")
        out_json = os.path.join(outdir, sample_id + ".json")
        with open(out_json, "w") as json_output:
            json.dump(sample_out, json_output, indent=4)


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


def update_config(gene_list, args):
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
        samtools_check = [
            a for a in [args.samtools, shutil.which("samtools")] if a is not None
        ]
        samtools_check2 = [a for a in samtools_check if os.path.exists(a)]
        if samtools_check2 == []:
            raise Exception("samtools is not found")
        else:
            configs[gene]["tools"].setdefault("samtools", samtools_check2[0])

        minimap2_check = [
            a for a in [args.minimap2, shutil.which("minimap2")] if a is not None
        ]
        minimap2_check2 = [a for a in minimap2_check if os.path.exists(a)]
        if minimap2_check2 == []:
            raise Exception("minimap2 is not found")
        else:
            configs[gene]["tools"].setdefault("minimap2", minimap2_check2[0])

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


def main():
    parser = argparse.ArgumentParser(
        description="paraphase: HiFi-based SMN1/SMN2 variant caller."
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
    outputp.add_argument(
        "-o",
        "--out",
        help="Output directory",
        required=True,
    )
    parser.add_argument(
        "-g",
        "--gene",
        help="Optionally specify which gene(s) to run (separated by comma). Will run all genes if not specified",
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

    args = parser.parse_args()
    outdir = args.out
    logging.basicConfig(level=logging.DEBUG)

    if os.path.exists(outdir) is False:
        os.makedirs(outdir)

    gene_list = args.gene
    if gene_list is None:
        gene_list = ACCEPTED_GENES
    else:
        gene_list = [a for a in gene_list.split(",") if a in ACCEPTED_GENES]
        if gene_list == []:
            joined_list = ",".join(ACCEPTED_GENES)
            raise Exception(
                f"Gene names not recognized. Currently accepted genes are {joined_list}"
            )

    configs = update_config(gene_list, args)

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
            process_sample(
                bamlist,
                outdir,
                configs,
                dcov,
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
            process_sample,
            outdir=outdir,
            configs=configs,
            dcov=dcov,
        )
        bam_groups = [bamlist[i::nCores] for i in range(nCores)]
        pool = mp.Pool(nCores)
        pool.map(process_sample_partial, bam_groups)
        pool.close()
        pool.join()
    else:
        raise Exception("No input file is given")


if __name__ == "__main__":
    main()
