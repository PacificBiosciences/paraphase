import pytest
import yaml
import os
from paraphase.genes.smn1_phaser import Smn1Phaser


class TestSmn1Phaser(object):

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
    phaser = Smn1Phaser(sample_id, sample_dir)
    phaser.set_parameter(config)

    def test_check_smn1_smn2_presence(self):
        self.phaser.check_smn1_smn2_presence()
        assert self.phaser.has_smn1 is True
        assert self.phaser.has_smn2 is True

    def test_get_long_del_reads(self):
        smn2_del_reads, smn2_del_reads_partial = self.phaser.get_long_del_reads(
            self.phaser.del2_3p_pos1,
            self.phaser.del2_3p_pos2,
            self.phaser.del2_5p_pos1,
            self.phaser.del2_5p_pos2,
            self.phaser.deletion2_size,
        )
        smn1_del_reads, smn1_del_reads_partial = self.phaser.get_long_del_reads(
            self.phaser.del1_3p_pos1,
            self.phaser.del1_3p_pos2,
            self.phaser.del1_5p_pos1,
            self.phaser.del1_5p_pos2,
            self.phaser.deletion1_size,
        )
        assert smn2_del_reads != set()
        assert smn2_del_reads_partial != set()
        assert smn1_del_reads == set()
        assert smn1_del_reads_partial == set()

    def test_allow_del_bases(self):
        (
            self.phaser.smn2_del_reads,
            self.phaser.smn2_del_reads_partial,
        ) = self.phaser.get_long_del_reads(
            self.phaser.del2_3p_pos1,
            self.phaser.del2_3p_pos2,
            self.phaser.del2_5p_pos1,
            self.phaser.del2_5p_pos2,
            self.phaser.deletion2_size,
        )
        assert self.phaser.allow_del_bases(70948287) is True
        assert self.phaser.allow_del_bases(70948285) is False

    def test_assign_haps_to_gene(self):
        self.phaser.het_sites = ["70951940_A_C", "70951946_T_G", "70951958_T_G"]
        haps = ["111", "121", "131"]
        smn1_haps, smn2_haps, smn2_del_haps = self.phaser.assign_haps_to_gene(haps)
        assert smn1_haps == ["111"]
        assert smn2_haps == ["121"]
        assert smn2_del_haps == ["131"]

        # no splice site variant, has smn1 but no smn2
        self.phaser.het_sites = ["70951947_A_C", "70951949_T_G", "70951958_T_G"]
        self.phaser.has_smn1 = True
        self.phaser.smn2_reads_splice = 0
        haps = ["111", "121", "131"]
        smn1_haps, smn2_haps, smn2_del_haps = self.phaser.assign_haps_to_gene(haps)
        assert smn1_haps == ["111", "121"]
        assert smn2_haps == []
        assert smn2_del_haps == ["131"]

        # no splice site variant, has smn2 but no smn1
        self.phaser.het_sites = ["70951947_A_C", "70951949_T_G", "70951958_T_G"]
        self.phaser.has_smn1 = False
        self.phaser.smn2_reads_splice = 25
        haps = ["111", "121", "131"]
        smn1_haps, smn2_haps, smn2_del_haps = self.phaser.assign_haps_to_gene(haps)
        assert smn1_haps == []
        assert smn2_haps == ["111", "121"]
        assert smn2_del_haps == ["131"]

        # no splice site variant, has both smn1 and smn2, returns empty
        self.phaser.het_sites = ["70951947_A_C", "70951949_T_G", "70951958_T_G"]
        self.phaser.has_smn1 = True
        self.phaser.smn2_reads_splice = 25
        haps = ["111", "121", "131"]
        smn1_haps, smn2_haps, smn2_del_haps = self.phaser.assign_haps_to_gene(haps)
        assert smn1_haps == []
        assert smn2_haps == []
        assert smn2_del_haps == ["131"]

    def test_adjust_smn1_cn(self):
        ass_haps = ["11", "22"]
        smn1_haps = {"11": "smn1hap1"}
        smn1_cn, _ = self.phaser.adjust_smn1_cn(2, 1, 2, ass_haps, None, smn1_haps)
        assert smn1_cn == 2

        self.phaser.mdepth = None
        smn1_cn, _ = self.phaser.adjust_smn1_cn(1, 1, 3, ass_haps, None, smn1_haps)
        assert smn1_cn is None

        self.phaser.has_smn2 = False
        smn1_cn, _ = self.phaser.adjust_smn1_cn(1, 1, 3, ass_haps, None, smn1_haps)
        assert smn1_cn == 2

        # genome depth too low
        self.phaser.mdepth = 15
        smn1_cn, _ = self.phaser.adjust_smn1_cn(1, 1, 2, ass_haps, None, smn1_haps)
        assert smn1_cn is None

        self.phaser.mdepth = 30
        self.phaser.smn1_reads_splice = 15
        smn1_cn, _ = self.phaser.adjust_smn1_cn(1, 1, 2, ass_haps, None, smn1_haps)
        assert smn1_cn == 1

        self.phaser.mdepth = 30
        self.phaser.smn1_reads_splice = 30
        smn1_cn, _ = self.phaser.adjust_smn1_cn(1, 1, 2, ass_haps, None, smn1_haps)
        assert smn1_cn == 2

        self.phaser.mdepth = 38
        self.phaser.smn1_reads_splice = 30
        smn1_cn, _ = self.phaser.adjust_smn1_cn(1, 1, 2, ass_haps, None, smn1_haps)
        assert smn1_cn is None

        # compare smn1 and smn2 depth
        self.phaser.mdepth = None
        self.phaser.smn1_reads_splice = 30
        self.phaser.smn2_reads_splice = 32
        smn1_cn, two_cp_haps = self.phaser.adjust_smn1_cn(
            1, 2, 2, ass_haps, None, smn1_haps
        )
        assert smn1_cn == 2
        assert two_cp_haps == ["smn1hap1"]

        # compare smn1 and smn2 depth, smn2_cn = 1, only updates to None
        self.phaser.mdepth = None
        self.phaser.smn1_reads_splice = 30
        self.phaser.smn2_reads_splice = 15
        smn1_cn, _ = self.phaser.adjust_smn1_cn(1, 1, 2, ass_haps, None, smn1_haps)
        assert smn1_cn is None

        # first check genome depth, get None and then compare smn1 and smn2 depth
        self.phaser.mdepth = 38
        self.phaser.smn1_reads_splice = 30
        self.phaser.smn2_reads_splice = 32
        smn1_cn, _ = self.phaser.adjust_smn1_cn(1, 2, 2, ass_haps, None, smn1_haps)
        assert smn1_cn == 2

    def test_adjust_smn2_cn(self):
        smn2_haps = {"22": "smn2hap1"}
        # genome depth too low
        self.phaser.mdepth = 15
        smn2_cn, _ = self.phaser.adjust_smn2_cn(1, 1, smn2_haps)
        assert smn2_cn is None

        # smn2, consider smn2_del
        self.phaser.mdepth = 30
        self.phaser.smn2_reads_splice = 30
        self.phaser.smn2_del_reads_partial = set()
        smn2_cn, _ = self.phaser.adjust_smn2_cn(1, 1, smn2_haps)
        assert smn2_cn == 2

        # smn2, when there is smn2_del, only one read
        self.phaser.mdepth = 30
        self.phaser.smn2_reads_splice = 30
        self.phaser.smn2_del_reads_partial = {"r1"}
        smn2_cn, _ = self.phaser.adjust_smn2_cn(1, 1, smn2_haps)
        assert smn2_cn == 2

        # smn2, when there is smn2_del, more than one read
        self.phaser.mdepth = 30
        self.phaser.smn2_reads_splice = 30
        self.phaser.smn2_del_reads_partial = {"r1", "r2"}
        smn2_cn, _ = self.phaser.adjust_smn2_cn(1, 1, smn2_haps)
        assert smn2_cn == 1

        # compare smn1 and smn2 depth, with smn2 del78
        self.phaser.mdepth = None
        self.phaser.smn1_reads_splice = 30
        self.phaser.smn2_reads_splice = 32
        smn2_cn, _ = self.phaser.adjust_smn2_cn(2, 1, smn2_haps)
        assert smn2_cn == 1

        # compare smn1 and smn2 depth, without smn2 del78
        self.phaser.smn2_del_reads_partial = set()
        self.phaser.mdepth = None
        self.phaser.smn1_reads_splice = 30
        self.phaser.smn2_reads_splice = 32
        smn2_cn, _ = self.phaser.adjust_smn2_cn(2, 1, smn2_haps)
        assert smn2_cn == 2

        # compare smn1 and smn2 depth, without smn2 del78
        # do not consider when smn1_cn is 1
        self.phaser.smn2_del_reads_partial = set()
        self.phaser.mdepth = None
        self.phaser.smn1_reads_splice = 15
        self.phaser.smn2_reads_splice = 32
        smn2_cn, _ = self.phaser.adjust_smn2_cn(1, 1, smn2_haps)
        assert smn2_cn == 1

        self.phaser.close_handle()
