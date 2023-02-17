# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


import numpy as np
import pysam


class GenomeDepth:
    def __init__(self, bam, config):
        self.bam = bam
        self.config = config
        self._bamh = pysam.AlignmentFile(bam, "rb")
        self.mdepth = None
        self.mad = None

    def get_genome_depth(self):
        depth = []
        with open(self.config["data"]["depth_region"]) as f:
            for line in f:
                at = line.split()
                nchr = at[0]
                pos1 = int(at[1])
                site_depth = 0
                for pos in [pos1, pos1 + 1600]:
                    for pileupcolumn in self._bamh.pileup(
                        nchr, pos - 1, pos, truncate=True
                    ):
                        site_depth = pileupcolumn.get_num_aligned()
                depth.append(site_depth)
        self.mdepth = np.median(depth)
        self.mad = np.median([abs(a - self.mdepth) for a in depth]) / self.mdepth

    def get_regional_depth(self):
        if (
            "depth_region" in self.config["coordinates"]["hg38"]
            and "nchr" in self.config["coordinates"]["hg38"]
        ):
            depth_region = self.config["coordinates"]["hg38"]["depth_region"]
            nchr = self.config["coordinates"]["hg38"]["nchr"]
            region_depth = []
            for region in depth_region:
                depth = []
                nstep = max(1, int((region[1] - region[0]) / 100))
                for pos in range(region[0], region[1], nstep):
                    for pileupcolumn in self._bamh.pileup(
                        nchr, pos - 1, pos, truncate=True
                    ):

                        site_depth = pileupcolumn.get_num_aligned()
                        depth.append(site_depth)
                region_depth.append(np.median(depth))
            return region_depth
        return None

    def call(self):
        self.get_genome_depth()
        region_depth = self.get_regional_depth()
        self._bamh.close()
        return self.mdepth, self.mad, region_depth
