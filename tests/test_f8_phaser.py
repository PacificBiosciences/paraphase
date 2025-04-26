import pytest
import yaml
import os
from paraphase.genes.f8_phaser import F8Phaser
from .test_phaser import update_config


class TestF8Phaser(object):
    cur_dir = os.path.dirname(__file__)
    sample_dir = os.path.join(cur_dir, "test_data", "f8")

    def test_inversion(self):
        sample_id = "inv"
        genome_bam = os.path.join(self.sample_dir, f"f8_{sample_id}_genome.bam")
        phaser = F8Phaser(
            sample_id, self.sample_dir, genome_bam=genome_bam, sample_sex="male"
        )
        config = update_config("f8")
        phaser.set_parameter(config)
        f8_call = phaser.call()
        assert f8_call.sv_called == {"f8_int22invhap1": "inversion"}

    def test_deletion(self):
        sample_id = "del"
        genome_bam = os.path.join(self.sample_dir, f"f8_{sample_id}_genome.bam")
        phaser = F8Phaser(
            sample_id, self.sample_dir, genome_bam=genome_bam, sample_sex="male"
        )
        config = update_config("f8")
        phaser.set_parameter(config)
        f8_call = phaser.call()
        assert f8_call.sv_called == {"f8_int22delhap1": "deletion"}
