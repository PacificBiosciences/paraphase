import pytest
import yaml
import os
from paraphase.phaser import Phaser


class TestPhaser(object):

    cur_dir = os.path.dirname(__file__)
    sample_dir = os.path.join(cur_dir, "test_data")
    sample_id = "HG00733"
    data_dir = os.path.join(os.path.dirname(cur_dir), "paraphase", "data")
    config_file = os.path.join(data_dir, "config.yaml")
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)
    config = config["smn1"]
    config.setdefault("gene", "smn1")
    realign_region = config["realign_region"]
    nchr = realign_region.split(":")[0]
    nchr_old = realign_region.replace(":", "_").replace("-", "_")
    config.setdefault("nchr", nchr)
    config.setdefault("nchr_old", nchr_old)
    data_paths = config.get("data")
    for data_entry in data_paths:
        old_data_file = data_paths[data_entry]
        new_data_file = os.path.join(data_dir, "smn1", old_data_file)
        data_paths[data_entry] = new_data_file
    data_paths.setdefault("reference", os.path.join(sample_dir, "ref.fa"))
    phaser = Phaser(sample_id, sample_dir)
    phaser.set_parameter(config)

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

    def test_check_linking_read(self):
        aln1 = "xxx12"
        aln2 = "11xxxx"
        assert self.phaser.check_linking_read(aln1, aln2) == "1-2"
        assert self.phaser.check_linking_read(aln1, aln2, reverse=True) is None

        aln1 = "11xxxx"
        aln2 = "xxx12"
        assert self.phaser.check_linking_read(aln1, aln2) == "2-1"
        assert self.phaser.check_linking_read(aln1, aln2, reverse=True) is None

        aln1 = "xxx12"
        aln2 = "xxxxxx11"
        assert self.phaser.check_linking_read(aln1, aln2) is None
        assert self.phaser.check_linking_read(aln1, aln2, reverse=True) == "0-0"

        aln1 = "xxx1x"
        aln2 = "1xxxxx11"
        assert self.phaser.check_linking_read(aln1, aln2) is None
        assert self.phaser.check_linking_read(aln1, aln2, reverse=True) is None

    def test_get_alleles_from_links(self):
        ass_haps = ["1", "2", "3", "4", "5", "6", "7"]
        read_links = {"1": ["2"], "2": ["1", "3"], "3": ["2"]}
        alleles = self.phaser.get_alleles_from_links(read_links, ass_haps)
        assert alleles == [["1", "2", "3"]]

        # test merge
        read_links = {
            "1": ["2"],
            "2": ["1", "3"],
            "4": ["5"],
            "5": ["4", "3"],
            "3": ["5"],
            "6": ["7"],
            "7": ["6"],
        }
        alleles = self.phaser.get_alleles_from_links(read_links, ass_haps)
        assert len(alleles) == 2 and (
            sorted(alleles[0]) == ["1", "2", "3", "4", "5"]
            or sorted(alleles[1]) == ["1", "2", "3", "4", "5"]
        )

    def test_get_directed_links(self):
        raw_read_haps = {
            "r1": "xx1",
            "r1_sup": "2xx",
            "r2": "xx1",
            "r2_sup": "2xx",
            "r3": "x1x",  # does not cover either end
            "r3_sup": "2xx",
        }
        new_reads = {
            "r1": [{"r1": ["111"]}, {"r1_sup": ["211"]}],
            "r2": [{"r2": ["111"]}, {"r2_sup": ["211"]}],
            "r3": [{"r3": ["111"]}, {"r3_sup": ["211"]}],
        }
        ass_haps = {"111": "hap1", "211": "hap2"}

        (
            nondirected_links,
            directed_links,
            directed_links_loose,
        ) = self.phaser.get_directed_links(new_reads, raw_read_haps, ass_haps, False)
        assert directed_links == {"hap1-hap2": [1, 1]}
        assert nondirected_links == {"hap1-hap2": [1, 1, 1]}
