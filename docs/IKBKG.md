# IKBKG

In this region, Paraphase calls small variants in IKBKG. In addition, there is a known 11.7kb deletion that can occur in either IKBKG or the pseudogene. Paraphase calls this deletion and locates it to IKBKG or the pseudogene.

## Fields in the `json` file

- `deletion_haplotypes`: haplotypes carrying the 11.7kb deletion

Note that this deletion is reported in the VCF as a structural variant (SV). 

## Visualizing haplotypes

To visualize phased haplotypes, load the output bam file in IGV, group reads by the `HP` tag and color alignments by `YC` tag. Green and purple represent two alleles, i.e. all haplotypes in green are on one one allele and all haplotypes in purple are on the other allele. Reads in gray are either unassigned or consistent with more than one possible haplotype. When two haplotypes are identical over a region, there can be more than one haplotype consistent with a read, and the read is randomly assigned to a haplotype and colored in gray. 

![IKBKG example](figures/RCCX.png)