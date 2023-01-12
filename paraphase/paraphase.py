# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


import os
import argparse
import json
import yaml
import logging
import datetime
import shutil
import multiprocessing as mp
from functools import partial
from paraphase.genome_depth import GenomeDepth
from paraphase.prepare_bam_and_vcf import (
    BamRealigner,
    BamTagger,
    TwoGeneVcfGenerater,
)
from paraphase.smn_phaser import SmnPhaser


def process_sample(bamlist, outdir, config, dcov={}, vcf=False):
    for bam in bamlist:
        sample_id = bam.split("/")[-1].split(".")[0]
        logging.info(f"Processing sample {sample_id} at {datetime.datetime.now()}...")

        logging.info(f"Getting genome depth at {datetime.datetime.now()}...")
        gdepth = None
        if sample_id in dcov:
            gdepth = dcov[sample_id]
        if gdepth is None:
            depth = GenomeDepth(bam, config)
            gdepth, gmad = depth.call()
            if gdepth < 10 or gmad > 0.25:
                logging.warning(
                    "Due to low or highly variable genome coverage, genome coverage is not used for depth correction."
                )
                gdepth = None

        logging.info(f"Realigning reads at {datetime.datetime.now()}...")
        bam_realigner = BamRealigner(bam, outdir, config)
        bam_realigner.write_realign_bam()

        logging.info(f"Phasing haplotypes at {datetime.datetime.now()}...")
        smn_phaser = SmnPhaser(sample_id, outdir, config, gdepth)
        smn_phaser_call = smn_phaser.call()._asdict()

        logging.info(f"Tagging reads at {datetime.datetime.now()}...")
        bam_tagger = BamTagger(sample_id, outdir, config, smn_phaser_call)
        bam_tagger.write_bam(random_assign=True)

        if vcf:
            logging.info(f"Generating VCFs at {datetime.datetime.now()}...")
            vcf_generater = TwoGeneVcfGenerater(
                sample_id, outdir, config, smn_phaser_call
            )
            vcf_generater.run()

        sample_out = smn_phaser_call
        logging.info(f"Writing to json at {datetime.datetime.now()}...")
        out_json = os.path.join(outdir, sample_id + ".json")
        with open(out_json, "w") as json_output:
            json.dump(sample_out, json_output, indent=4)


def main():
    parser = argparse.ArgumentParser(
        description="paraphase: HiFi-based SMN1/SMN2 variant caller."
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
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
    parser.add_argument(
        "-o",
        "--out",
        help="Output directory",
        required=True,
    )
    parser.add_argument(
        "-v",
        "--vcf",
        help="Optional. If specified, paraphase will produce vcf for each haplotype",
        required=False,
        action="store_true",
    )
    parser.add_argument(
        "-c",
        "--config",
        help="Optional path to config yaml file",
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
    config_file = args.config
    logging.basicConfig(level=logging.DEBUG)

    if os.path.exists(outdir) is False:
        os.makedirs(outdir)

    data_path = os.path.join(os.path.dirname(__file__), "data")
    if config_file is None:
        config_file = os.path.join(data_path, "smn1", "config.yaml")

    # parse config file
    with open(config_file, "r") as f:
        try:
            config = yaml.safe_load(f)
        except yaml.YAMLError as yaml_error:
            raise Exception(f"Error reading {config_file}: \n\n{yaml_error}")

    # check third-party tools
    tools = config.get("tools")
    samtools_check = [
        a
        for a in [tools.get("samtools"), args.samtools, shutil.which("samtools")]
        if a is not None
    ]
    samtools_check2 = [a for a in samtools_check if os.path.exists(a)]
    if samtools_check2 == []:
        raise Exception("samtools is not found")
    else:
        config["tools"]["samtools"] = samtools_check2[0]

    minimap2_check = [
        a
        for a in [tools.get("minimap2"), args.minimap2, shutil.which("minimap2")]
        if a is not None
    ]
    minimap2_check2 = [a for a in minimap2_check if os.path.exists(a)]
    if minimap2_check2 == []:
        raise Exception("minimap2 is not found")
    else:
        config["tools"]["minimap2"] = minimap2_check2[0]

    # update paths
    gene = config.get("gene")
    data_paths = config.get("data")
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
                config,
                dcov,
                args.vcf,
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
            config=config,
            dcov=dcov,
            vcf=args.vcf,
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
