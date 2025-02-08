import pytest
import os
import json
from paraphase.prepare_bam_and_vcf import VcfGenerater
from .test_phaser import update_config


class TestVcfGenerater(object):
    cur_dir = os.path.dirname(__file__)
    sample_dir = os.path.join(cur_dir, "test_data", "vcf")

    def test_get_var(self):
        all_bases = ["A"] * 10
        var_seq, dp, ad, gt, qual, counter = VcfGenerater.get_var(all_bases, "A")
        assert gt == "0"
        assert var_seq == "A"
        assert dp == 10
        assert ad == (10, 10)

        all_bases = []
        var_seq, dp, ad, gt, qual, counter = VcfGenerater.get_var(all_bases, "A")
        assert gt == "."
        assert var_seq == "A"
        assert dp == 0
        assert ad == (0, 0)

        all_bases = ["T"] * 2
        var_seq, dp, ad, gt, qual, counter = VcfGenerater.get_var(all_bases, "A")
        assert gt == "1"
        assert var_seq == "T"
        assert dp == 2
        assert ad == (0, 2)

        all_bases = ["*"] * 3
        var_seq, dp, ad, gt, qual, counter = VcfGenerater.get_var(all_bases, "A")
        assert gt == "1"
        assert var_seq == "*"
        assert dp == 3
        assert ad == (0, 3)

        all_bases = ["T", "T", "A"]
        var_seq, dp, ad, gt, qual, counter = VcfGenerater.get_var(all_bases, "A")
        assert gt == "1"
        assert var_seq == "T"
        assert dp == 3
        assert ad == (1, 2)

        all_bases = ["A+2T", "A+2T", "A"]
        var_seq, dp, ad, gt, qual, counter = VcfGenerater.get_var(all_bases, "A")
        assert gt == "1"
        assert var_seq == "A+2T"
        assert dp == 3
        assert ad == (1, 2)

        all_bases = ["A+2T", "A+2T", "A+2T", "A+2T", "A"]
        var_seq, dp, ad, gt, qual, counter = VcfGenerater.get_var(all_bases, "A")
        assert gt == "1"
        assert var_seq == "A+2T"
        assert dp == 5
        assert ad == (1, 4)

    def test_convert_alt_record(self):
        assert VcfGenerater.convert_alt_record("T", "A") == "A"
        assert VcfGenerater.convert_alt_record("T", "TACG") == "T+3ACG"
        assert VcfGenerater.convert_alt_record("TGC", "T") == "T-2NN"

    def test_modify_hapbound(self):
        assert VcfGenerater.modify_hapbound(1, 2, None) == "1-2"
        assert VcfGenerater.modify_hapbound(1, 2, ["5p"]) == "1truncated-2"
        assert VcfGenerater.modify_hapbound(1, 2, ["3p"]) == "1-2truncated"
        assert (
            VcfGenerater.modify_hapbound(1, 2, ["5p", "3p"]) == "1truncated-2truncated"
        )

    def test_run_without_realign(self):
        sample_id = "HG004"
        with open(os.path.join(self.sample_dir, "HG004.paraphase.json")) as f:
            phase_calls = json.load(f)

        # homozygous case
        config = update_config("ARL17A")
        vcf_generater = VcfGenerater(
            sample_id,
            self.sample_dir,
            phase_calls["ARL17A"],
        )
        vcf_generater.set_parameter(config, tmpdir=self.sample_dir, prog_cmd="test")
        variants_info, hap_info = vcf_generater.run_without_realign()
        assert len(hap_info) == 2
        assert hap_info[0][0] == "ARL17A_homozygous_hap1"
        assert hap_info[1][0] == "ARL17A_homozygous_hap1_cp2"

        # two-copy haplotypes
        # truncated haplotypes
        config = update_config("AMY2A")
        vcf_generater = VcfGenerater(
            sample_id,
            self.sample_dir,
            phase_calls["AMY2A"],
        )
        vcf_generater.set_parameter(config, tmpdir=self.sample_dir, prog_cmd="test")
        variants_info, hap_info = vcf_generater.run_without_realign()
        assert hap_info == [
            ["AMY2A_hap1", 103616000, 103631602, []],
            ["AMY2A_hap1_cp2", 103616000, 103631602, []],
            ["AMY2A_hap2", 103619306, 103631602, ["5p"]],
            ["AMY2A_hap3", 103619306, 103631602, ["5p"]],
        ]

        # ikbkg deletion, big SV
        config = update_config("ikbkg")
        vcf_generater = VcfGenerater(
            sample_id,
            self.sample_dir,
            phase_calls["ikbkg"],
        )
        vcf_generater.set_parameter(config, tmpdir=self.sample_dir, prog_cmd="test")
        variants_info, hap_info = vcf_generater.run_without_realign()
        assert 154558014 in variants_info
        assert ["154558014_DEL_154569698", ".", ".", [], "1", None] in variants_info[
            154558014
        ]
