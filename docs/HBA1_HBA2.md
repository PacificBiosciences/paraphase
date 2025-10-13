# HBA1/HBA2


This [medically relevant region](https://www.ncbi.nlm.nih.gov/books/NBK1435/), is characterised by two regions of homology, represented by "A" and "B" in this simplified schematic of the region:
![HBA1_HBA2 Region](figures/HBA1-HBA2-diagram.png)

Paraphase returns the total copy number of HBA1 and HBA2 in this region. Variants are called in the VCF relative to the HBA2 reference sequence.


#### Structural Variants: 3.7 kb and 4.2 kb Deletions

Two well-known structural variants in this region are the **3.7 kb** and **4.2 kb deletions**, which occur due to **homologous recombination** between the homologous sequences:

- **3.7 kb deletion (or duplication)**  
  Results from recombination between the "B" boxes, forming a hybrid HBA1/HBA2 gene.

- **4.2 kb deletion (or duplication)**  
  Results from recombination between the "A" boxes, leading to a deletion (or duplication) of HBA2.


## Fields in the `json` file

Fields shared across all genes are defined in the general json [file](json.png). The region includes several unique fields:

- `genotype`: Reports the genotype for this family. Possible alleles include:
  - `aa`: wild-type
  - `aaa`: duplication
  - `-a`: single-gene deletion
  - `--`: double-gene deletion
- `surrounding_region_depth`: xx to do xx
- **`sv_called`**: Reports structural variants:
  - `3p7del`: 3.7 kb deletion
  - `3p7dup`: 3.7 kb duplication
  - `4p2del`: 4.2 kb deletion
  - `4p2dup`: 4.2 kb duplication  
  Along with their coordinates.
  
Paraphase labels the HBA1 and HBA2 haplotypes accordingly. Paraphase also returns a `homology_haplotype` that 
aligns to the "A" region closest to HBA2 and corresponds to the reads mapping to the left "A" region. This haplotype 
should be clipped on both sides in wild-type and is used to detect the 4.2 deletion (appears left-clipped only).

## Visualizing haplotypes

To visualize phased haplotypes, load the output bam file in IGV, group reads by the `HP` tag and color alignments by `YC` tag. Green and purple represent two alleles, i.e. all haplotypes in green are on one one allele and all haplotypes in purple are on the other allele. 

Reads in gray are either unassigned or consistent with more than one possible haplotype. When two haplotypes are identical over a region, there can be more than one haplotype consistent with a read, and the read is randomly assigned to a haplotype and colored in gray. 

![HBA example](figures/HBA.png)

- The top panel shows a sample with two copies of HBA1 and two copies of HBA2, one on each allele. 
- The bottom panel shows a sample with a `-a` allele, where there is a deletion, leaving only one copy of HBA (`hba_del_hap1`).
