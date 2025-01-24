# Running Paraphase on targeted data

Paraphase can work with targeted data, such as:
- Hybrid capture based enrichment data
- CRISPR-Cas9 targeted data
- Amplicon data

The config file may need to be modified based on the design of the target panel. Please reach out to Xiao Chen (xchen@pacificbiosciences.com) if you need assistance.

Paraphase provides a few options for users to better work with targeted data: 
1) Use the `--targeted` option to drop the assumption of uniform coverage across the genome.
2) Additionally there are two optional parameters one can tune for targeted data:
- `--min-read-variant`: Partially controls the number of supporting reads for a variant for identifying variants used for phasing. The cutoff for variant-supporting reads is determined by min(this number, max(5, depth\*0.11)). Default is 20. At standard WGS depth, the default value is overwritten by max(5, depth*0.11). For targeted data with high coverage, set this number relatively high to avoid picking up sequencing errors and to reduce run time. For example, if you expect your per-haplotype depth is 200, you can set `--min-read-variant` to 40 or even higher.
- `--min-read-haplotype`: Minimum number of unique supporting reads for a haplotype. Default is 4. For targeted data with high coverage, this cutoff can be increased to reduce errors and to reduce run time.
