# Paraphase demo

Perform a quick Paraphase run with the following dataset and commands.

```bash
# Download human GRCh38 if you don't have one
wget https://downloads.pacbcloud.com/public/reference-genomes/human_GRCh38_no_alt_analysis_set.tar.2023-12-04.gz
tar -xpvf human_GRCh38_no_alt_analysis_set.tar.2023-12-04.gz
# Get demo input BAM by cloning the Paraphase repo
git clone https://github.com/PacificBiosciences/paraphase
# Run Paraphase for the SMN1/SMN2 region
paraphase -b ./tests/data/bams/HG01175_smn1_extracted.bam -r ./human_GRCh38_no_alt_analysis_set/human_GRCh38_no_alt_analysis_set.fasta -o ./output/ -p HG01175 -g smn1
```

Note that a warning message `Genome-wide coverage is too low or too variable; skipping depth-based correction.` is expected as this test dataset is not a WGS BAM, but a bamlet extracted from the original WGS BAM. If you would like to run the full WGS BAM for this sample, it can be downloaded from:
```bash
# BAM
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG01175/analysis/aligned_reads/hifi/GRCh38/HG01175_aligned_GRCh38_winnowmap.sorted.bam
# index
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG01175/analysis/aligned_reads/hifi/GRCh38/HG01175_aligned_GRCh38_winnowmap.sorted.bam.bai
```

## Check Paraphase outputs
Paraphase output files can be found in `./output/`:
- `HG01175.paraphase.bam` and `HG01175.paraphase.bam.bai`
- `HG01175.paraphase.json`
- `HG01175_paraphase_vcfs` folder

`HG01175.paraphase.bam` can be visualized as described in the [tutorial on SMN1/SMN2](SMN1_SMN2.md):
![HG01175 IGV](figures/HG01175_smn1.png)

`HG01175.paraphase.json` contains haplotype and copy number calls. Some fields under `region_specific_info` are listed below, showing two copies of SMN1, one copy of SMN2 and one copy of SMN with Exons7-8 deleted:
```json
"smn1": {
        "region_specific_info": {
            "smn1_cn": 2,
            "smn1_haplotypes": {
                "1121211121221222121211111111221111111111111111111112111111111211111111111111111": "smn1_smn1hap1",
                "2112111212211111211111111111111211111211122211111112111111112111111111111111111": "smn1_smn1hap2"
            },
            "smn1_read_number": 37,
            "smn2_cn": 1,
            "smn2_haplotypes": {
                "1121122121122222121122222222212122222122211122222221222221222122222222222222222": "smn1_smn2hap1"
            },
            "smn2_read_number": 18,
            "smn_del78_cn": 1,
            "smn_del78_haplotypes": {
                "2211111212211111112121213333333333333333333322222221122222222122222222222222222": "smn1_smndel78hap1"
            },
            "smn_del78_read_number": 22
        },
    }
```