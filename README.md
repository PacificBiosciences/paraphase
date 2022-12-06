# Paraphase: HiFi-based SMN1/SMN2 variant caller

SMN1, the gene that causes spinal muscular atrophy, is considered a 'dark' region of the genome due to high sequence similarity with its paralog SMN2. Paraphase is a Python tool that takes HiFi BAMs as input (whole-genome or enrichment), phases complete SMN1 and SMN2 haplotypes, determines copy numbers and makes phased variant calls for both genes. It also categorizes the haplotypes, enabling future haplotype-based screening of silent carriers (2+0). Please check out our [preprint](https://www.biorxiv.org/content/10.1101/2022.10.19.512930) for more details about the method and our population-wide haplotype analysis.   

For whole-genome sequencing (WGS) data, we recommend >20X, ideally 30X, genome coverage. Low coverage or short read length could result in less accurate phasing, especially when haplotypes are highly similar to each other in Exons 1-6. For hybrid capture-based enrichment data, a higher read depth (>50X) is recommended as the read length is generally shorter than WGS.

## Contact

If you need assistance or have suggestions, please don't hesitate to reach out by email or open a GitHub issue.

Xiao Chen: xchen@pacificbiosciences.com

## Dependencies

- [samtools](http://www.htslib.org/)
- [minimap2](https://github.com/lh3/minimap2)

## Installation

Paraphase can be installed through pip:
```bash
pip install paraphase
```

Alternatively, Paraphase can be installed from GitHub.
```bash
git clone https://github.com/PacificBiosciences/paraphase
cd paraphase
python setup.py install
```

## Running the program

```bash
paraphase -b input.bam -o output_directory
```

Alternatively when you have a list of bam files
```bash
paraphase -l list.txt -o output_directory
```

Required parameters:
- `-b`: Input BAM file or `-l`: List of BAM files (one per line)
- `-o`: Output directory

Optional parameters:
- `-v`: If specified, Paraphase will produce VCFs for each haplotype.
- `-c`: Config file, default config file is `paraphase/data/smn1/config.yaml`.
- `-t`: Number of threads, used when `-l` is specified.
- `-d`: File listing average genome depth per sample, with two columns, sample ID and depth values, separated by tab or space. This saves run time by skipping the step to calculate genome depth.
- `--samtools`
- `--minimap2`
- `--singularity`

The paths to samtools and minimap2 can be provided through the `--samtools` and `--minimap2` parameters or by modifying the `tools` section of the [config](paraphase/data/smn1/config.yaml) file.

Note that currently only GRCh38 is supported. We will support GRCh37 in the future if there is request.

## Interpreting the output

Paraphase produces a few output files in the directory specified by `-o`, with the sample ID as the prefix.
- `_realigned_tagged.bam`: This BAM file can be loaded into IGV for visualization of haplotypes, see [haplotype visualization](docs/visualization.md).  
- If `-v` is specified, Paraphase will generate VCF files. A VCF file is written for each haplotype, and there is also a `_variants.vcf` file containing merged variants from all haplotypes.
- `.json`: Main output file, summerizes haplotypes and variant calls for each sample. Details of the fields are explained below:
  - `smn1_cn`: copy number of SMN1, a `null` call indicates that Paraphase finds only one haplotype but depth does not unambiguously support a copy number of one or two.
  - `smn2_cn`: copy number of SMN2, a `null` call indicates that Paraphase finds only one haplotype but depth does not unambiguously support a copy number of one or two.
  - `smn2_del78_cn`: copy number of SMN2Δ7–8 (SMN2 with a deletion of Exon7-8)
  - `smn1_read_number`: number of reads containing c.840C
  - `smn2_read_number`: number of reads containing c.840T
  - `smn2_del78_read_number`: number of reads containing the known deletion of Exon7-8 on SMN2
  - `smn1_haplotypes`: phased SMN1 haplotypes
  - `smn2_haplotypes`: phased SMN2 haplotypes
  - `smn2_del78_haplotypes`: phased SMN2Δ7–8 haplotypes
  - `two_copy_haplotypes`: haplotypes that are present in two copies based on depth. This happens when (in a small number of cases) two haplotypes are identical and we infer that there exist two of them instead of one by checking the read depth.
  - `haplotype_details`: lists information about each haplotype 
    - `variants`: The variants contained in the haplotype, excluding those in homopolymer regions. For a complete set of variant calls, please use the `-v` option.
    - `boundary`: The boundary of the region that is resolved on the haplotype. This is useful when a haplotype is only partially phased.
    - `haplogroup`: The haplogroup that the haplotype is assigned to

