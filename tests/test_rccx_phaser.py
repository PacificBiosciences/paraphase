import pytest
import os
from paraphase.genes.rccx_phaser import RccxPhaser


class TestRccxPhaser(object):

    cur_dir = os.path.dirname(__file__)
    sample_dir = os.path.join(cur_dir, "test_data")
    sample_id = "HG00733"
    phaser = RccxPhaser(sample_id, sample_dir)

    def test_annotate_var(self):
        allele_var = [[], []]
        annotated_allele = self.phaser.annotate_var(allele_var)
        assert annotated_allele == "WT"

        allele_var = [[], ["var1"]]
        annotated_allele = self.phaser.annotate_var(allele_var)
        assert annotated_allele == "WT"

        allele_var = [["var1"], ["var1", "var2"]]
        annotated_allele = self.phaser.annotate_var(allele_var)
        assert annotated_allele == "var1"

        allele_var = [["var1", "var2", "var3"], ["var1", "var2"]]
        annotated_allele = self.phaser.annotate_var(allele_var)
        assert annotated_allele == "var1,var2"

        allele_var = [[]]
        annotated_allele = self.phaser.annotate_var(allele_var)
        assert annotated_allele == "pseudogene_deletion"

        allele_var = [["var1"]]
        annotated_allele = self.phaser.annotate_var(allele_var)
        assert annotated_allele == "deletion_var1"

        allele_var = [["var1", "var2"], [], []]
        annotated_allele = self.phaser.annotate_var(allele_var)
        assert annotated_allele == "gene_duplication"

        allele_var = [
            ["var1", "var2", "var3", "var4", "var5", "var6"],
            ["var1", "var2", "var3", "var4", "var5", "var6"],
            [],
        ]
        annotated_allele = self.phaser.annotate_var(allele_var)
        assert annotated_allele == "pseudogene_duplication"

        allele_var = [["var1", "var2"], ["var1"], []]
        annotated_allele = self.phaser.annotate_var(allele_var)
        assert annotated_allele == "duplication_WT_plus_var1"

        allele_var = [
            ["var1", "var2", "var3", "var4", "var5", "var6"],
            ["var1", "var2", "var3", "var4", "var5", "var6"],
            ["var1", "var2"],
        ]
        annotated_allele = self.phaser.annotate_var(allele_var)
        assert annotated_allele == "var1,var2_pseudogene_duplication"

        allele_var = [
            ["var1", "var2", "var3", "var4", "var5", "var6"],
            ["var1", "var2", "var3", "var4"],
            ["var1", "var2"],
        ]
        annotated_allele = self.phaser.annotate_var(allele_var)
        assert annotated_allele == "duplication_var1,var2_plus_var1,var2,var3,var4"

        allele_var = [
            [],
            [],
            [],
            ["var1", "var2"],
        ]
        annotated_allele = self.phaser.annotate_var(allele_var)
        assert annotated_allele is None

    def test_annotate_alleles(self):
        ending_copies = []
        ass_haps = ["hap1", "hap2", "hap3", "hap4"]
        two_cp_haplotypes = []
        alleles = [["hap1", "hap2"], ["hap3", "hap4"]]
        hap_variants = {
            "hap1": [],
            "hap2": [],
            "hap3": [],
            "hap4": [],
        }
        annotated_alleles = self.phaser.annotate_alleles(
            True,
            alleles,
            hap_variants,
            ending_copies,
            ass_haps,
            two_cp_haplotypes,
        )
        assert annotated_alleles == ["WT", "WT"]

        annotated_alleles = self.phaser.annotate_alleles(
            False,
            alleles,
            hap_variants,
            ending_copies,
            ass_haps,
            two_cp_haplotypes,
        )
        assert annotated_alleles == []

        ending_copies = ["hap1", "hap2"]
        annotated_alleles = self.phaser.annotate_alleles(
            False,
            alleles,
            hap_variants,
            ending_copies,
            ass_haps,
            two_cp_haplotypes,
        )
        assert annotated_alleles == ["WT", "WT"]

        hap_variants = {
            "hap1": ["var1"],
            "hap2": [],
            "hap3": [],
            "hap4": [],
        }
        ending_copies = ["hap1", "hap2"]
        annotated_alleles = self.phaser.annotate_alleles(
            False,
            alleles,
            hap_variants,
            ending_copies,
            ass_haps,
            two_cp_haplotypes,
        )
        assert annotated_alleles == ["var1", "WT"]
