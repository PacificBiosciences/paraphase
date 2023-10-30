import pytest
import os
from paraphase.genes.opn1lw_phaser import Opn1lwPhaser


class TestOpn1lwPhaser(object):

    cur_dir = os.path.dirname(__file__)
    sample_dir = os.path.join(cur_dir, "test_data")
    sample_id = "HG00733"
    phaser = Opn1lwPhaser(sample_id, sample_dir)

    def test_phase_single_allele(self):
        first_copy = "h1"
        # two copies
        all_haps_on_this_allele = ["h1", "h2"]
        last_copies = ["h2"]
        middle_copies = []
        hap_links = {}
        alleles = []
        directional_links = []
        directional_links_loose = {}

        phased_allele = self.phaser.phase_single_allele(
            first_copy,
            all_haps_on_this_allele,
            last_copies,
            middle_copies,
            hap_links,
            alleles,
            directional_links,
            directional_links_loose,
        )
        assert phased_allele == ["h1", "h2"]

        # three copies, last one known
        all_haps_on_this_allele = ["h1", "h2", "h3"]
        last_copies = ["h3"]
        middle_copies = []
        hap_links = {}
        alleles = []
        directional_links = []
        directional_links_loose = {}
        phased_allele = self.phaser.phase_single_allele(
            first_copy,
            all_haps_on_this_allele,
            last_copies,
            middle_copies,
            hap_links,
            alleles,
            directional_links,
            directional_links_loose,
        )
        assert phased_allele == ["h1", "h2"]

        # three total copies. we don't know the last copy but we know the middle copy
        all_haps_on_this_allele = ["h1", "h2", "h3"]
        last_copies = []
        middle_copies = ["h3"]
        hap_links = {}
        alleles = []
        directional_links = []
        directional_links_loose = {}
        phased_allele = self.phaser.phase_single_allele(
            first_copy,
            all_haps_on_this_allele,
            last_copies,
            middle_copies,
            hap_links,
            alleles,
            directional_links,
            directional_links_loose,
        )
        assert phased_allele == ["h1", "h3"]

        # four copies and we know the one before the last copy
        all_haps_on_this_allele = ["h1", "h2", "h3", "h4"]
        last_copies = ["h4"]
        middle_copies = []
        hap_links = {"h4": ["h2"]}
        alleles = []
        directional_links = []
        directional_links_loose = {}
        phased_allele = self.phaser.phase_single_allele(
            first_copy,
            all_haps_on_this_allele,
            last_copies,
            middle_copies,
            hap_links,
            alleles,
            directional_links,
            directional_links_loose,
        )
        assert phased_allele == ["h1", "h3"]

        # four total copies. we know the directional link between the middle two copies
        all_haps_on_this_allele = ["h1", "h2", "h3", "h4"]
        last_copies = ["h4"]
        middle_copies = []
        hap_links = {}
        alleles = []
        directional_links = ["h2-h3"]
        directional_links_loose = {}
        phased_allele = self.phaser.phase_single_allele(
            first_copy,
            all_haps_on_this_allele,
            last_copies,
            middle_copies,
            hap_links,
            alleles,
            directional_links,
            directional_links_loose,
        )
        assert phased_allele == ["h1", "h2"]

        # if other haps are phased on this allele
        all_haps_on_this_allele = ["h1", "h2", "h3", "h4"]
        last_copies = ["h4"]
        middle_copies = []
        hap_links = {}
        alleles = [["h2", "h4"]]
        directional_links = []
        directional_links_loose = {}
        phased_allele = self.phaser.phase_single_allele(
            first_copy,
            all_haps_on_this_allele,
            last_copies,
            middle_copies,
            hap_links,
            alleles,
            directional_links,
            directional_links_loose,
        )
        assert phased_allele == ["h1", "h3"]

        # loose links of the first copy
        all_haps_on_this_allele = ["h1", "h2", "h3", "h4"]
        last_copies = ["h4"]
        middle_copies = []
        hap_links = {}
        alleles = []
        directional_links = []
        directional_links_loose = {"h1": {"h2": 3}}
        phased_allele = self.phaser.phase_single_allele(
            first_copy,
            all_haps_on_this_allele,
            last_copies,
            middle_copies,
            hap_links,
            alleles,
            directional_links,
            directional_links_loose,
        )
        assert phased_allele == ["h1", "h2"]
