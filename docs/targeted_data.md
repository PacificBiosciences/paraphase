# Running Paraphase on targeted data

Paraphase can work with targeted data, such as:
- Hybrid capture based enrichment data
- CRISPR-Cas9 targeted data
- Amplicon data

The config file may need to be modified based on the design of the target panel. Please reach out to Xiao Chen (xchen@pacificbiosciences.com) if you need assistance.

Paraphase provides a few options for users to better work with targeted data: 
1) Use the `--targeted` option to drop the assumption of uniform coverage across the genome.
2) Additionally there are two optional parameters designed for targeted data. The default values are expected to work well with high depth data since they are frequency-based. But users can tune them based on the depth of their data and the expected copy numbers of their regions of interest. 
- `--min-variant-frequency`:  Minimum frequency for a variant to be used for phasing. The cutoff for variant-supporting reads is determined by max(5, total_depth * min_frequency). Note that total_depth is the combined depth of all paralogs for a paralog group. Default is 0.11. 
- `--min-haplotype-frequency`: Minimum frequency of unique supporting reads for a haplotype. The cutoff for haplotype-supporting reads is determined by max(4, total_depth * min_frequency). Note that total_depth is the combined depth of all paralogs for a paralog group. Default is 0.03. This cutoff can be increased to filter out spurious and low-frequency haplotypes.
