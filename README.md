# Paraphase: HiFi-based caller for highly homologous genes

Many medically relevant genes fall into 'dark' regions where variat calling is limited due to high sequence homology with paralogs or pseudogenes. Paraphase is a Python tool that takes HiFi BAMs as input (whole-genome or enrichment), phases complete haplotypes for genes of the same family, determines copy numbers and makes phased variant calls. 

![Paraphase diagram](docs/figures/paraphase_diagram.png)
Paraphase takes all reads from a gene family, realigns to just the gene of interest and then phases them into haplotypes. This solves the problem of alignment difficulty due to sequence homology and allows us to examine all copies of genes in a gene family and call copy number changes and other variants.

Paraphase supports 161 segmental duplication regions in GRCh38, listed in the [config](paraphase/data/38/config.yaml) file. Among these, there are 11 medically relevant regions that are also supported in GRCh37/hg19:
- SMN1/SMN2 (spinal muscular atrophy)
- RCCX module
  - CYP21A2 (21-Hydroxylase-Deficient Congenital Adrenal Hyperplasia)
  - TNXB (Ehlers-Danlos syndrome)
  - C4A/C4B (relevant in autoimmune diseases)
- PMS2 (Lynch Syndrome)
- STRC (hereditary hearing loss and deafness)
- IKBKG (Incontinentia Pigmenti)
- NCF1 (chronic granulomatous disease; Williams syndrome)
- NEB (Nemaline myopathy)
- F8 (intron 22 inversion, Hemophilia A)
- CFC1 (heterotaxy syndrome)
- OPN1LW/OPN1MW (color vision deficiencies)
- HBA1/HBA2 (Alpha-Thalassemia)

Please check out our paper on its application to the gene SMN1 for more details about Paraphase.   
Chen X, Harting J, Farrow E, et al. Comprehensive SMN1 and SMN2 profiling for spinal muscular atrophy analysis using long-read PacBio HiFi sequencing. The American Journal of Human Genetics. 2023;0(0). doi:10.1016/j.ajhg.2023.01.001

For whole-genome sequencing (WGS) data, we recommend >20X, ideally 30X, genome coverage. Low coverage or short read length could result in less accurate phasing, especially when haplotypes are highly similar to each other in Exons 1-6. For hybrid capture-based enrichment data, a higher read depth (>50X) is recommended as the read length is generally shorter than WGS.

## Contact

If you have suggestions or need assistance, please don't hesitate to reach out by email or open a GitHub issue.

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
paraphase -b input.bam -o output_directory -r genome_fasta
```

Alternatively when you have a list of bam files
```bash
paraphase -l list.txt -o output_directory -r genome_fasta
```

Required parameters:
- `-b`: Input BAM file or `-l`: text file listing BAM files one per line (a BAI file needs to exist in the same directory)
- `-o`: Output directory
- `-r`: Path to the reference genome fasta file

Please note that the input BAM should be one that's aligned to the ENTIRE reference genome (either GRCh38 or GRCh37/hg19), and this reference should NOT include ALT contigs. The fasta file of this reference genome should be provided to Paraphase with `-r`. 

Optional parameters:
- `-g`: Gene(s) to analyze. All supported genes will be analyzed if not specified.
- `-t`: Number of threads.
- `-d`: File listing average genome depth per sample, with two columns, sample ID and depth values, separated by tab or space. This saves run time by skipping the step to calculate genome depth.
- `--genome`: Genome reference build. Default is `38`. If `37` or `19` is specified, Paraphase will run the analysis for GRCh37/hg19 (note that only 11 medically relevant regions are supported now).
- `gene1only`: If specified, variants calls will be made against the main gene only for SMN1, PMS2, STRC, NCF1 and IKBKG, see [below](#interpreting-the-output).
- `--novcf`: no vcf output if specified.
- `--samtools`: path to samtools. If the paths to samtools or minimap2 are not in the PATH environment variable, they can be provided through the `--samtools` and `--minimap2` parameters.
- `--minimap2`: path to minimap2

## Interpreting the output

Paraphase produces a few output files in the directory specified by `-o`, with the sample ID as the prefix.

1. `.vcf` in `sampleID_vcfs` folder. A VCF file is written for each haplotype per gene family. There is also a `_variants.vcf` file containing merged variants from all haplotypes. Note that this is not a diploid vcf as there are usually more than 2 copies of each gene family in a sample.

As genes of the same family can be highly similar to each other in sequence and not easy to differentiate, variant calls are made against one selected gene from the gene family (e.g. the gene is selected when the family has a gene and a pseudogene). Exceptions are SMN1 (paralog SMN2), PMS2 (pseudogene PMS2CL), STRC (pseudogene STRCP1), NCF1 (pseudogenes NCF1B and NCF1C) and IKBKG (pseudogene IKBKGP1). In these families, haplotypes are assigned to each gene in the family, i.e. gene or paralog/pseudogene, and variants are called against the gene (or paralog/pseudogene) for the gene (or paralog/pseudogene) haplotypes, respectively. Variants calls can be made against the gene only for these five families if `--gene1only` is specified. 

2. `_realigned_tagged.bam`: This BAM file can be loaded into IGV for visualization of haplotypes (group reads by `HP` tag and color alignments by `YC` tag). All haplotypes are aligned against the gene of interest. Tutorials/Examples are provided for medically relevant genes (See below).  

3. `.json`: Output file summerizing haplotypes and variant calls for each sample. In brief, a few generally used fields are explained below.
- `final_haplotypes`: phased haplotypes
- `total_cn`: total copy number of the family (sum of gene and paralog/pseudogene)
- `two_copy_haplotypes`: haplotypes that are present in two copies based on depth. This happens when (in a small number of cases) two haplotypes are identical and we infer that there exist two of them instead of one by checking the read depth.
- `haplotype_details`: lists information about each haplotype 
  - `boundary`: the boundary of the region that is resolved on the haplotype. This is useful when a haplotype is only partially phased.
- `alleles_final`: haplotypes phased into alleles, when repeat are in tandem.
- `region_depth`: median depth of the family (include all copies of gene and paralog/pseudogene) 

Tutorials/Examples are provided for interpreting medically relevant genes: 
- [SMN1/SMN2](docs/SMN1_SMN2.md)
- [RCCX module (CYP21A2)](docs/RCCX.md)
- [PMS2](docs/PMS2.md)
- [STRC](docs/STRC.md)
- [OPN1LW/OPN1MW](docs/OPN1LW_OPN1MW.md)
- [HBA1/HBA2](docs/HBA1_HBA2.md)
- [IKBKG](docs/IKBKG.md)
- [F8](docs/F8.md)
- [NEB](docs/NEB.md)
- [NCF1](docs/NCF1.md)
