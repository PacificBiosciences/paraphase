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

    def call(self):
        self.get_genome_depth()
        self._bamh.close()
        return self.mdepth, self.mad
