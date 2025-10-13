# Paraphase json

Paraphase produces a single json file per sample, containing data across all analyzed regions.

This file includes summary information about the gene family, as well as details on haplotype phasing. Certain 
structural variants—such as gene fusions or large gene-specific deletions (e.g., in HBA1)—are also reported in the json 
file. In contrast, small variants like SNPs and InDels should be interpreted from the accompanying `vcf` [file](vcf.md). 
Haplotypes present in two copies will be reported as separate entries in the vcf.

As described in the `vcf` [file](vcf.md) documentation, Paraphase can assign and label haplotypes to specific genes 
within a family when with prior knowledge on paralog differentiation is provided. Certain json fields are unique to 
families where such prior knowledge has been used.

## json fields

The json file contains fields that common to all genes and fields specific to individual genes. Gene specific fields are 
documented separately in each gene-specific documentation files. Common fields are described here:

### General fields:

- `total_cn`: total copy number of the entire gene family (sum of gene and paralog/pseudogene).
- `gene_cn`: total copy number for the genes (not including pseudogenes). Present only when paralog differentiation is available.
- `two_copy_haplotypes`: haplotypes inferred to exist in two identical copies based on depth, even if only one sequence is observed. 
- `genome_depth`: average sequencing depth across the whole genome.
- `region_depth`: median and 80th percentile of the total cumulative depth of all haplotypes in the region.
- `sample_sex`: inferred biological sex based on coverage of X and Y chromosomes.
- `fusions_called`: deletions or duplications from unequal crossing over between paralogs, called by checking the 
- flanking sequences of phased haplotypes. Currently supported for: CYP2D6, GBA, CYP11B1, and the CFH gene cluster.
- `highest_total_cn`: internal-use field, to be removed in future releases.

### Phasing information
  - `sites_for_phasing` variant sites used for phasing
  - `final_haplotypes` list of phased haplotypes for all gene copies in the family. Each haplotype is encoded based on 
the positions in sites_for_phasing using the following scheme:
    - `1`: ref
    - `2`: variant
    - `3` for deletion
    - `x` for unknown
  - `haplotype_details`: metadata for each haplotype, including:
    - `boundary`: Start and end positions for the resolved haplotype region (useful when only part of a haplotype is phased).
    - `is_truncated`: contains the information whether the haplotype is truncated at the 3 or 5' end
    - `variants`: the phasing variants observed in that haplotype
  - `read_details`: for each read aligning to the gene family, includes assignment info at each `sites_for_phasing` position.
  - `unique_supporting_reads`: reads that uniquely support one of the haplotype in `final_haplotypes`.
  - `nonunique_supporting_reads`: reads that match to more than one haplotype and cannot uniquely be assigned.
  - `phase_region`: coordinates of the region over which phasing was performed.
  - `alleles_final`: haplotypes phased into alleles (this is possible when the segmental duplication is in tandem).
  - `assembled_haplotypes`: initial set of phased haplotypes before downstream filtering. Some haplotypes may be shorter 
than the full query region because phasing ends where sequence homology stops. For positions outside the resolved region 
of a haplotype, variant values are marked as `0`.
  - `het_sites_not_used_in_phasing`: heterozygous sites not used in phasing, either because they are indels or located 
in low confidence regions
  - `heterozygous_sites`: sites that differ between haplotypes
  - `homozygous_sites`: sites that are different from the reference, but identical across reads
  - `raw_alleles`: initial (pre-filtering) haplotype phasing into alleles
  - `haplotype links`: pairs of haplotypes supported by read-level evidence indicating linkage. When present, these 
links help phase haplotypes into alleles.
