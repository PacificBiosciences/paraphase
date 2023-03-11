# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


import numpy as np
import pysam


class GenomeDepth:
    def __init__(self, bam, genome_depth_region_file):
        self.bam = bam
        self.genome_depth_region_file = genome_depth_region_file
        self._bamh = pysam.AlignmentFile(bam, "rb")
        self.mdepth = None
        self.mad = None

    def get_genome_depth(self):
        depth = []
        with open(self.genome_depth_region_file) as f:
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
