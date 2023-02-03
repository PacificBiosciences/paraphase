# Paraphase: HiFi-based caller for highly homologous genes

Many medically relevant genes fall into 'dark' regions where variat calling is limited due to high sequence homology with paralogs or pseudogenes. Paraphase is a Python tool that takes HiFi BAMs as input (whole-genome or enrichment), phases complete haplotypes for genes of the same family, determines copy numbers and makes phased variant calls. 

This branch is work in progress as we adds more genes into Paraphase. Currently it supports SMN1/SMN2, the RCCX module (CYP21A2, TNXB, C4A/C4B) and PMS2.

Please check out our paper on its application to the gene SMN1 for more details about Paraphase.   
Chen X, Harting J, Farrow E, et al. Comprehensive SMN1 and SMN2 profiling for spinal muscular atrophy analysis using long-read PacBio HiFi sequencing. The American Journal of Human Genetics. 2023;0(0). doi:10.1016/j.ajhg.2023.01.001

For whole-genome sequencing (WGS) data, we recommend >20X, ideally 30X, genome coverage. Low coverage or short read length could result in less accurate phasing, especially when haplotypes are highly similar to each other in Exons 1-6. For hybrid capture-based enrichment data, a higher read depth (>50X) is recommended as the read length is generally shorter than WGS.

Currently Paraphase only works on GRCh38. Support for GRCh37 will be adde in the future.

## Contact

If you need assistance or have suggestions, please don't hesitate to reach out by email or open a GitHub issue.

Xiao Chen: xchen@pacificbiosciences.com

## Dependencies

- [samtools](http://www.htslib.org/)
- [minimap2](https://github.com/lh3/minimap2)

## Installation

Paraphase can be installed through pip or conda:
```bash
pip install paraphase
# or
conda install -c conda-forge -c bioconda paraphase
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
- `-g`: Gene(s) to analyze. All supported genes will be analyzed if not specified.
- `-t`: Number of threads, used when `-l` is specified.
- `-d`: File listing average genome depth per sample, with two columns, sample ID and depth values, separated by tab or space. This saves run time by skipping the step to calculate genome depth.
- `--samtools`
- `--minimap2`

The paths to samtools and minimap2 can be provided through the `--samtools` and `--minimap2` parameters.

## Interpreting the output

Paraphase produces a few output files in the directory specified by `-o`, with the sample ID and the gene name as the prefix.
- `_realigned_tagged.bam`: This BAM file can be loaded into IGV for visualization of haplotypes, see [haplotype visualization](docs/visualization.md).  
- `.vcf`: A VCF file is written for each haplotype, and there is also a `_variants.vcf` file containing merged variants from all haplotypes.
- `.json`: Main output file, summerizes haplotypes and variant calls for each sample. Details of the fields are explained below for each gene.

### SMN1

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
  - `variants`: The variants contained in the haplotype, excluding those in homopolymer regions.
  - `boundary`: The boundary of the region that is resolved on the haplotype. This is useful when a haplotype is only partially phased.
  - `haplogroup`: The haplogroup that the haplotype is assigned to

### RCCX & PMS2

- `total_cn`: total copy number of the family (sum of gene and paralog/pseudogene)
- `final_haplotypes`: phased haplotypes
- `two_copy_haplotypes`: haplotypes that are present in two copies based on depth. This happens when (in a small number of cases) two haplotypes are identical and we infer that there exist two of them instead of one by checking the read depth.
- `phasing_success`: whether haplotypes are phased into alleles in the case of RCCX
- `alleles_final`: haplotypes phased into alleles in the case of RCCX
- `annotated_alleles`: CYP21A2 allele annotation in the case of RCCX. This is only based on common gene-pseudogene conversions. Please refer to the vcfs for most thorough variant calling and annotation.





