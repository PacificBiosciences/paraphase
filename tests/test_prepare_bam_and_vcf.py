import pytest
import os
from paraphase.prepare_bam_and_vcf import VcfGenerater


class TestVcfGenerater(object):
    cur_dir = os.path.dirname(__file__)
    sample_dir = os.path.join(cur_dir, "test_data")
    sample_id = "HG00733"
    vcf_generater = VcfGenerater(
        sample_id,
        sample_dir,
        None,
    )

    def test_get_var(self):
        all_bases = ["A"] * 10
        var_seq, dp, ad, gt = VcfGenerater.get_var(all_bases, "A")
        assert gt == "0"
        assert var_seq == "A"
        assert dp == 10
        assert ad == 10

        all_bases = []
        var_seq, dp, ad, gt = VcfGenerater.get_var(all_bases, "A")
        assert gt == "."
        assert var_seq == "A"
        assert dp == 0
        assert ad == 0

        all_bases = ["T"] * 2
        var_seq, dp, ad, gt = VcfGenerater.get_var(all_bases, "A")
        assert gt == "."
        assert var_seq == "T"
        assert dp == 2
        assert ad == 2

        all_bases = ["T"] * 3
        var_seq, dp, ad, gt = VcfGenerater.get_var(all_bases, "A")
        assert gt == "1"
        assert var_seq == "T"
        assert dp == 3
        assert ad == 3

        all_bases = ["T", "T", "A"]
        var_seq, dp, ad, gt = VcfGenerater.get_var(all_bases, "A")
        assert gt == "1"
        assert var_seq == "T"
        assert dp == 3
        assert ad == 2

        all_bases = ["A+2T", "A+2T", "A"]
        var_seq, dp, ad, gt = VcfGenerater.get_var(all_bases, "A")
        assert gt == "."
        assert var_seq == "A+2T"
        assert dp == 3
        assert ad == 2

        all_bases = ["A+2T", "A+2T", "A+2T", "A+2T", "A"]
        var_seq, dp, ad, gt = VcfGenerater.get_var(all_bases, "A")
        assert gt == "1"
        assert var_seq == "A+2T"
        assert dp == 5
        assert ad == 4
