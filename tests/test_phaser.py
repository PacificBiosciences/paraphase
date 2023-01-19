import pytest
import yaml
import os
from paraphase.phaser import Phaser


class TestPhaser(object):

    cur_dir = os.path.dirname(__file__)
    sample_dir = os.path.join(cur_dir, "test_data")
    sample_id = "HG00733"
    data_dir = os.path.join(os.path.dirname(cur_dir), "paraphase", "data")
    config_file = os.path.join(data_dir, "smn1", "config.yaml")
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)
    data_paths = config.get("data")
    for data_entry in data_paths:
        old_data_file = data_paths[data_entry]
        new_data_file = os.path.join(data_dir, "smn1", old_data_file)
        data_paths[data_entry] = new_data_file
    phaser = Phaser(sample_id, sample_dir, config)

    def test_depth_prob(self):
        prob = Phaser.depth_prob(40, 20)
        assert prob.index(max(prob)) == 1

    def test_get_pivot_site_index(self):
        self.phaser.het_sites = ["70951940_A_C", "70951946_T_G"]
        pos, found_splice = self.phaser.get_pivot_site_index()
        assert pos == 1
        assert found_splice is True

        self.phaser.het_sites = ["70951940_A_C", "70951946_T_G", "70951958_T_G"]
        pos, found_splice = self.phaser.get_pivot_site_index()
        assert pos == 1
        assert found_splice is True

        self.phaser.het_sites = ["70951947_A_C", "70951949_T_G", "70951958_T_G"]
        pos, found_splice = self.phaser.get_pivot_site_index()
        assert pos == -1
        assert found_splice is False

    def test_process_indel(self):
        ref_seq, var_seq, indel_size = self.phaser.process_indel(
            70940935, "A", "A+3CCC"
        )
        assert ref_seq == "A"
        assert var_seq == "ACCC"
        assert indel_size == 3

        ref_seq, var_seq, indel_size = self.phaser.process_indel(70940935, "A", "A-2NN")
        assert ref_seq == "ACT"
        assert var_seq == "A"
        assert indel_size == 2

    def test_simplify_read_haps(self):
        haps = {"r1": ["1", "2"], "r2": ["2", "1"]}
        haplotypes_to_reads, reads_to_haplotypes = Phaser.simplify_read_haps(haps)
        assert haplotypes_to_reads == {"12": ["r1"], "21": ["r2"]}
        assert reads_to_haplotypes == {"r1": "12", "r2": "21"}

    def test_get_start_end(self):
        hap = "xx1211xxx"
        nstart, nend = Phaser.get_start_end(hap)
        assert nstart == 2
        assert nend == 5

    def test_get_hap_variant_ranges(self):
        self.phaser.het_sites = [
            "70917101_A_C",
            "70917111_A_C",
            "70917150_A_C",
            "70917200_A_C",
            "70917300_A_C",
            "70917400_A_C",
        ]
        hap = "x121xx"
        n1, n2 = self.phaser.get_hap_variant_ranges(hap)
        assert n1 == 70917102
        assert n2 == 70917299

        hap = "121211"
        n1, n2 = self.phaser.get_hap_variant_ranges(hap)
        assert n1 == 70917100
        assert n2 == 70961220

    def test_update_reads_for_deletions(self):
        self.phaser.het_sites = [
            "70917101_A_C",
            "70917111_A_C",
            "70917150_A_C",
            "70917200_A_C",
            "70917300_A_C",
            "70917400_A_C",
        ]
        del_reads_partial = {"r1"}
        raw_read_haps = {
            "r1": "1" * 6,
        }
        raw_read_haps, het_sites = self.phaser.update_reads_for_deletions(
            raw_read_haps,
            self.phaser.het_sites,
            70917199,
            70917201,
            del_reads_partial,
            "3",
            "del",
        )
        assert raw_read_haps["r1"] == "111311"
        assert het_sites == [
            "70917101_A_C",
            "70917111_A_C",
            "70917150_A_C",
            "70917200_A_C",
            "70917300_A_C",
            "70917400_A_C",
        ]

        # del at beginning, insert del at beginning
        self.phaser.het_sites = [
            "70917101_A_C",
            "70917111_A_C",
            "70917150_A_C",
            "70917200_A_C",
            "70917300_A_C",
            "70917400_A_C",
        ]
        del_reads_partial = {"r1"}
        raw_read_haps = {
            "r1": "1" * 6,
        }
        raw_read_haps, het_sites = self.phaser.update_reads_for_deletions(
            raw_read_haps,
            self.phaser.het_sites,
            70917100,
            70917100,
            del_reads_partial,
            "3",
            "del",
        )
        assert raw_read_haps["r1"] == "3111111"
        assert het_sites == [
            "del",
            "70917101_A_C",
            "70917111_A_C",
            "70917150_A_C",
            "70917200_A_C",
            "70917300_A_C",
            "70917400_A_C",
        ]

        # del at the end, do nothing?
        self.phaser.het_sites = [
            "70917101_A_C",
            "70917111_A_C",
            "70917150_A_C",
            "70917200_A_C",
            "70917300_A_C",
            "70917400_A_C",
        ]
        del_reads_partial = {"r1"}
        raw_read_haps = {
            "r1": "1" * 6,
        }
        raw_read_haps, het_sites = self.phaser.update_reads_for_deletions(
            raw_read_haps,
            self.phaser.het_sites,
            70917500,
            70917600,
            del_reads_partial,
            "3",
            "del",
        )
        assert raw_read_haps["r1"] == "111111"
        assert het_sites == [
            "70917101_A_C",
            "70917111_A_C",
            "70917150_A_C",
            "70917200_A_C",
            "70917300_A_C",
            "70917400_A_C",
        ]
