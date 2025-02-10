# paraphase
# Author: Xiao Chen <xchen@pacificbiosciences.com>


import numpy as np
import pysam
import logging
from paraphase.prepare_bam_and_vcf import pysam_handle

logging.basicConfig(level=logging.INFO)


class GenomeDepth:
    def __init__(self, bam, genome_depth_region_file, reference_fasta):
        self.bam = bam
        self.genome_depth_region_file = genome_depth_region_file
        self._bamh = pysam_handle(bam, reference_fasta)
        self.mdepth = None
        self.mad = None
        self.x = []
        self.y = []
        self.chr_in_name = GenomeDepth.check_chr_name(self._bamh)

    @staticmethod
    def check_chr_name(bam_handle):
        """
        Check bam header to see if chromosome names contain "chr"
        """
        header = bam_handle.header
        header = header.to_dict()
        sq = header.get("SQ")
        if sq is not None:
            chr_names = [a.get("SN") for a in sq]
            if True not in ["chr" in a for a in chr_names]:
                return False
        return True

    def get_genome_depth(self):
        depth = []
        with open(self.genome_depth_region_file) as f:
            for line in f:
                if line[0] != "#":
                    at = line.split()
                    nchr = at[0]
                    if self.chr_in_name is False:
                        nchr = nchr.replace("chr", "")
                    pos1 = int(at[1])
                    for pos in [pos1, pos1 + 1600]:
                        try:
                            site_depth = self._bamh.count(
                                nchr, pos - 1, pos, read_callback="all"
                            )
                            depth.append(site_depth)
                        except Exception:
                            logging.info(
                                f"This BAM does not have data at {nchr}:{pos}, which is used for determining sample coverage."
                            )
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
                if line[0] != "#":
                    at = line.split()
                    nchr = at[0]
                    if self.chr_in_name is False:
                        nchr = nchr.replace("chr", "")
                    pos1 = int(at[1])
                    try:
                        site_depth = self._bamh.count(
                            nchr, pos1 - 1, pos1, read_callback="all"
                        )
                        if "X" in nchr:
                            self.x.append(site_depth)
                        elif "Y" in nchr:
                            self.y.append(site_depth)
                    except Exception:
                        logging.info(
                            f"This BAM does not have data at {nchr}:{pos1}, which is used for determining sample sex."
                        )
        cov_x, cov_y = np.median(self.x), np.median(self.y)
        self._bamh.close()

        if np.isnan(cov_x) or np.isnan(cov_y) or cov_x == 0:
            return None
        if cov_y / cov_x < 0.05:
            return "female"
        elif cov_y / cov_x > 0.1:
            if cov_x > 1.95 * cov_y:
                return "female"
            return "male"
        return None
