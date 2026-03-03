# Variant annotation downstream of Paraphase

To address highly similar genes and common copy number changes in segmental duplication regions, the rationale of [Paraphase](https://github.com/PacificBiosciences/paraphase) is to treat gene copies of the same family as haplotypes of the same gene (pseudogenes are considered nonfunctional copies of the gene). Therefore, variants are called against the main gene selected to represent each gene family. The Paraphase VCF file contains variants called on each haplotype of a family and there are often more than two haplotypes in a family, and hence this is no longer a diploid problem. As an annotation step downstream of Paraphase, we recommend annotating the Paraphase VCF files for variant consequences and then determining the number of haplotypes free of pathogenic variants in each gene family. For recessive diseases, zero copy of functional haplotypes is analogous to the traditional scenario of biallelic mutations.

The `interpret_paraphase.py` script here provides a proof of concept for annotating Paraphase results. This script processes Paraphase VCF files and JSON output, identifies pathogenic variants using Ensembl VEP (Variant Effect Predictor) and ClinVar, and reports the number of functional haplotypes, i.e. those without pathogenic variants, per paralog region. Note that currently the script will only consider a variant as pathogenic if it has [high impact](https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html) by VEP, or is reported as `pathogenic` or `likely pathogenic` with multiple non-conflicting submissions in ClinVar. Users are welcome to modify the criteria for pathogenicity or plug in their own variant annotation engine.

Please note that this is only a proof of concept aiming to provide an example of how to interpret Paraphase results. Please reach out by email (xchen@pacificbiosciences.com) if you have suggestions or need assistance.

## Requirements

- Python 3.6+
- Required Python packages:
  - `requests`
  - `pyyaml`
  - `urllib3`

## Installation

1. This script is available as you clone or download [Paraphase](https://github.com/PacificBiosciences/paraphase).
2. Ensure the transcript file `paraphase_genes_transcript_start_end.txt` is present in the `data/` directory relative to the script `interpret_paraphase.py`. The script will automatically use this transcript file.

## Usage

### Basic Usage

```bash
python interpret_paraphase.py -i <paraphase_json_file> -o <output_directory>
```

### Required Arguments

- `-i, --input`: Path to Paraphase JSON file
  - Currently only supports Paraphase output from a GRCh38 run
  - The VCF files must be in a directory named `{sampleID}_paraphase_vcfs` in the same folder as the JSON file
  - Example: If the JSON file is called `sample1.paraphase.json`, VCFs should be in `sample1_paraphase_vcfs/`

- `-o, --out`: Output directory

### Optional Arguments

- `-r, --region`: Region(s) to annotate (comma-separated)
  - The default is to annotate predefined medical regions: `rccx`, `smn1`, `pms2`, `strc`, `cfc1`, `ikbkg`, `ncf1`, `neb`, `f8`, `hba`, `TNXB`, `OTOA`, `GBA`, `CFHclust`
  - Or specify individual regions: `rccx,smn1,pms2`
  - Or use `all` to annotate all regions in the config file (this could take a few hours)

- `-c, --config`: Path to Paraphase config file (YAML format)
  - If not provided, the script uses the default config from the Paraphase GitHub repository

- `-l, --log-level`: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
  - Default: `INFO`

- `--api-delay`: Minimum seconds between API requests for rate limiting
  - Default: `0.34` (~3 requests/second)
  - Increase if you encounter "Too Many Requests" errors

### Examples

**Annotate medical regions:**
```bash
python interpret_paraphase.py \
  -i sample1.paraphase.json \
  -o ./output \
```

**Annotate specific regions:**
```bash
python interpret_paraphase.py \
  -i sample1.paraphase.json \
  -o ./output \
  -r rccx,smn1,pms2
```

**Annotate all regions in the config:**
```bash
python interpret_paraphase.py \
  -i sample1.paraphase.json \
  -o ./output \
  -r all
```

**With custom config and debug logging:**
```bash
python interpret_paraphase.py \
  -i sample1.paraphase.json \
  -o ./output \
  -c /path/to/config.yaml \
  -l DEBUG
```

**With increased API rate limiting:**
```bash
python interpret_paraphase.py \
  -i sample1.paraphase.json \
  -o ./output \
  --api-delay 0.5
```

## Output Files

The script generates two CSV files in the output directory:

### Summary File (`{sampleID}_summary.csv`)

Contains per-gene summary statistics:
- `Region`: Region name
- `Region_median_depth`: Median depth of the region
- `Gene`: Gene symbol
- `Total_haplotypes`: Total number of haplotypes found by Paraphase
- `Num_haplotypes_with_pathogenic_variants`: Number of haplotypes containing pathogenic variants
- `Num_haplotypes_without_pathogenic_variants`: Number of haplotypes free of pathogenic variants (considered functional). **For clinically significant result, look for 0 in this field.**
- `Useful_information_from_Paraphase_JSON`: Region-specific, clinically relevant information from Paraphase JSON

### Detailed File (`{sampleID}_detailed.csv`)

Contains detailed variant annotations per haplotype:
- `Region`: Region name
- `Gene`: Gene symbol
- `Transcript`: Transcript ID
- `Transcript_start`: Transcript start coordinate
- `Transcript_end`: Transcript end coordinate
- `Haplotype`: Haplotype name
- `Haplotype_boundary`: Haplotype boundary coordinates reported by Paraphase
- `Boundary_check`: Overlap status between the haplotype and the transcript of interest. Possible values are:
  - `fully_contained`: The transcript region is fully contained in the phased haplotype, i.e. this haplotype covers the entire region needed to evaluate the gene.
  - `partial`: The haplotype only covers a portion of the transcript region. If this haplotype does not have a pathogenic variant, it is still possible that a pathogenic variant might exist in the unphased region of the gene copy.
  - `truncated`: The haplotype is truncated (softclipped) in the middle of the transcript region, indicating the end of sequence homology. A truncated haplotype is considered a nonfunctional gene copy - it most likely comes from a pseudogene.
  - `no_overlap`: The phased haplotype does not overlap with the transcript region at all, indicating that this is an irrelevant haplotype for the gene of interest. 
- `Pathogenic_variants_annotated`: Annotated pathogenic variants with ClinVar IDs, descriptions, and VEP consequences

## Rate Limiting

The script includes built-in rate limiting to prevent API rate limit errors:
- Default: ~3 requests/second (0.34 second delay)
- Automatic retry with exponential backoff on 429 errors
- Respects `Retry-After` headers from APIs
- Configurable via `--api-delay` parameter

## Notes

- The script requires internet connectivity to query Ensembl VEP and ClinVar APIs.
- Processing time depends on the number of variants and API response times.
- Large batches may take significant time due to rate limiting.
- This script currently supports GRCh38 only.
- As annotation is done per variant, this script does not consider the combined effect of multiple variants on the same haplotype. For example, a 1bp deletion and a 1bp insertion in the same codon will both be reported out as pathogenic, while their combined effect is likely not pathogenic. This is a shortcoming of this proof-of-concept script, as it would be very useful to make use of the long phasing information present in Paraphase haplotypes.
- No small variants are annotated for `f8` and `CFHclust`, as the main utility of Paraphase for them is to call SVs. The script reports those SV calls made by Paraphase in the summary.csv.  

## Useful links

- Paraphase documentation: https://github.com/PacificBiosciences/paraphase
- Ensembl VEP API: https://rest.ensembl.org/documentation/info/vep_hgvs_get
- ClinVar API: https://www.ncbi.nlm.nih.gov/clinvar/docs/api_http/

