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
        self.x = []
        self.y = []

    def get_genome_depth(self):
        depth = []
        with open(self.genome_depth_region_file) as f:
            for line in f:
                at = line.split()
                nchr = at[0]
                pos1 = int(at[1])
                for pos in [pos1, pos1 + 1600]:
                    site_depth = self._bamh.count(
                        nchr, pos - 1, pos, read_callback="all"
                    )
                    depth.append(site_depth)
        self.mdepth = np.median(depth)
        self.mad = np.median([abs(a - self.mdepth) for a in depth]) / self.mdepth

    def call(self):
        self.get_genome_depth()
        self._bamh.close()
        return self.mdepth, self.mad

    def check_sex(self):
        """Determine sample sex based on coverage"""
        with open(self.genome_depth_region_file) as f:
            for line in f:
                at = line.split()
                nchr = at[0]
                pos1 = int(at[1])
                site_depth = self._bamh.count(nchr, pos1 - 1, pos1, read_callback="all")
                if "X" in nchr:
                    self.x.append(site_depth)
                elif "Y" in nchr:
                    self.y.append(site_depth)
        cov_x, cov_y = np.median(self.x), np.median(self.y)
        self._bamh.close()
        if cov_y / cov_x < 0.05:
            return "female"
        elif cov_y / cov_x > 0.1:
            return "male"
        return None
