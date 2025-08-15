# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


from ..phaser import Phaser


class CfhClust(Phaser):
    def __init__(
        self,
        sample_id,
        outdir,
        cfh,
        cfhr3,
    ):
        Phaser.__init__(self, sample_id, outdir)
        self.cfh = cfh
        self.cfhr3 = cfhr3

    def call(self):
        haps = {}
        two_cp_haps = []
        fusions = {}
        total_cn = None
        if (
            self.cfh["final_haplotypes"] is not None
            and self.cfhr3["final_haplotypes"] is not None
        ):
            haps.update(self.cfh["final_haplotypes"])
            haps.update(self.cfhr3["final_haplotypes"])
        if (
            self.cfh["fusions_called"] is not None
            and self.cfhr3["fusions_called"] is not None
        ):
            fusions.update(self.cfh["fusions_called"])
            fusions.update(self.cfhr3["fusions_called"])
        if (
            self.cfh["two_copy_haplotypes"] is not None
            and self.cfhr3["two_copy_haplotypes"] is not None
        ):
            two_cp_haps += self.cfh["two_copy_haplotypes"]
            two_cp_haps += self.cfhr3["two_copy_haplotypes"]
        if self.cfh["total_cn"] is not None and self.cfhr3["total_cn"] is not None:
            total_cn = min(self.cfh["total_cn"], self.cfhr3["total_cn"])
            if (
                fusions != {}
                and len(self.cfh["final_haplotypes"]) >= 2
                and len(self.cfhr3["final_haplotypes"]) >= 2
            ):
                total_cn = min(
                    total_cn,
                    len(self.cfh["final_haplotypes"]),
                    len(self.cfhr3["final_haplotypes"]),
                )

        return self.GeneCall(
            total_cn,
            None,
            haps,
            two_cp_haps,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            f"{self.cfh["phase_region"]},{self.cfhr3["phase_region"]}",
            None,
            fusions,
        )
