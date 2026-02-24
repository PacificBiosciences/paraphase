# Paraphase json

Paraphase produces a single json file per sample, containing data across all analyzed regions.

This file includes summary information about the gene family, as well as details on haplotype phasing. Certain structural variants called from the presence of hybrid haplotypes, such as SVs in the `fusions_called` field described below and SVs in HBA1/2, RCCX or F8, are reported in the json file, as they are per-allele events instead of per-haplotype events that can be described in the vcf. In contrast, small variants like SNPs and InDels should be interpreted from the Paraphase [vcf file](vcf.md).

As described in the [vcf file](vcf.md) documentation, Paraphase can assign and label haplotypes to specific genes within a family in some regions where prior knowledge on paralog differentiation is provided. Certain json fields are unique to families where such prior knowledge has been used.

## json fields

The json file contains fields that are common to all Paraphase target regions and fields specific to individual region. Region-specific fields are documented separately in region-specific tutorial files. Common fields are described here:

### General fields:

- `total_cn`: total copy number of the entire gene family (sum of gene and paralog/pseudogene).
- `gene_cn`: total copy number for the gene (not including pseudogenes). Present only when Paraphase performs paralog differentiation and assigns haplotypes to genes or pseudogenes.
- `two_copy_haplotypes`: haplotypes inferred to exist in two identical copies based on depth. This happens when there is a haplotype that is present at twice the depth of other haplotypes and we infer that there exist two of them instead of one.
- `phase_region`: coordinates of the region over which phasing was performed. Can be used to quickly locate the region when visualizing Paraphase BAM with IGV.
- `genome_depth`: average sequencing depth across the whole genome.
- `region_depth`: median and 80th percentile of the depth across positions in the target region, where the depth is the combined depth of all reads, i.e. considering genes and paralogs together.
- `sample_sex`: inferred biological sex based on coverage of X and Y chromosomes.
- `genes`: gene symbols for the genes and paralogs/pseudogenes analyzed in this region.

### Information on phased haplotypes
  - `sites_for_phasing` variant sites used for phasing (currently only using SNPs for phasing).
  - `final_haplotypes` list of phased haplotypes for all gene copies in the family. Each haplotype is encoded based on the positions in `sites_for_phasing` using the following scheme:
    - `1`: reference base
    - `2`: variant base
    - `3` and `4` for positions inside a deletion, i.e. the region spanning the position is deleted on the haplotype. Currently at most two big deletions are allowed in a region, and they are labeled 3 and 4. The deletions are either called de novo by Paraphase or pre-specified in the config file. Note that these deletions are for better phasing of haplotypes, not for the purpose of calling deletions.
    - `x` for unknown. Could be due to low base quality, noisy alignment or missing coverage.
    - `0`: Some haplotypes may be shorter than the full query region because their alignment ends where sequence homology stops, i.e. they become soft-clipped. For positions beyond the softclips, variant values are marked as `0`.
  - `haplotype_details`: metadata for each haplotype, including:
    - `boundary`: Start and end positions for the resolved haplotype region (useful when only part of a haplotype is phased).
    - `is_truncated`: contains the information on whether the haplotype is truncated at the 3 or 5' end, indicating the end of homology region (soft-clipped).
    - `variants`: variants observed in the haplotype. For most accurate calling of variants, we recommend referring to the Paraphase vcf.
  - `read_details`: for each read aligning to the region, includes assignment info at each `sites_for_phasing` position.
  - `unique_supporting_reads`: reads that uniquely support one of the haplotypes in `final_haplotypes`.
  - `nonunique_supporting_reads`: reads that match to more than one haplotype and cannot uniquely be assigned. Lists all haplotypes that each read is consistent with.
  - `assembled_haplotypes`: initial set of phased haplotypes before downstream filtering.
  - `het_sites_not_used_in_phasing`: heterozygous sites not used in phasing, either because they are indels or located in low confidence regions.
  - `heterozygous_sites`: sites that differ between haplotypes.
  - `homozygous_sites`: sites that are different from the reference, but identical across reads.
  - `highest_total_cn`: internal-use field indicating the highest number of haplotypes seen across a region. This value could be higher than `total_cn` due to noisy alignments.

### Phasing haplotypes into alleles
  - `alleles_final`: haplotypes phased into alleles (this is possible when the segmental duplication is in tandem).
  - `raw_alleles`: initial (pre-filtering) alleles phased among haplotypes.
  - `haplotype links`: pairs of haplotypes supported by read-level evidence indicating linkage. When present, these links help phase haplotypes into alleles.

### Structural variants
- `fusions_called`: deletions or duplications from unequal crossing over between paralogs, called by checking the flanking sequences of phased haplotypes. Currently supported for: CYP2D6, GBA, CYP11B1, and the CFH gene cluster.