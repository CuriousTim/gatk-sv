"""
Set BCF genotypes to null where samples have low coverage.

Usage: python3 null_low_coverage_gt.py <inbcf> <outbcf> <bincov> <targets> <mincov>

<inbcf>    VCF/BCF to check.
<outbcf>   Where to write the modified VCF/BCF.
<bincov>   Binned coverage matrix. Must be indexed.
<targets>  TSV with sample to modify in first column, sample whose coverage
           should be checked in the second column.
<mincov>   Minimum median coverage of region for a genotype to be kept.
"""

import sys

import pandas as pd
from pysam import VariantFile
from pysam import TabixFile


class RdMatReader:
    def __init__(self, path):
        self.handle = TabixFile(path)
        # skip the first three fields '#Chr\tStart\tEnd'
        self.header = self.handle.header[0].lstrip("#").split("\t")[3:]

    def fetch(self, contig, start, end):
        for rec in self.handle.fetch(contig, start, end):
            yield RdMatReader._parse_rec(rec)

    def fetch_as_df(self, contig, start, end):
        return pd.DataFrame.from_records(
            self.fetch(contig, start, end), columns=self.header
        )

    def _parse_rec(rec):
        return tuple(int(x) for x in rec.split("\t")[3:])


def read_sample_map(path):
    tmp = dict()
    with open(path, mode="r", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip().split("\t")
            tmp[parts[0]] = parts[1]

    return tmp


def filter_gt_by_coverage(inbcf, outbcf, rdmat, sample_map, min_cov):
    for rec in inbcf.fetch():
        if rec.info["SVTYPE"] == "DEL":
            outbcf.write(rec)
            continue
        cov = rdmat.fetch_as_df(
            rec.contig, rec.pos - 1, rec.pos + rec.info["SVLEN"] - 1
        )
        median_covs = cov.median(axis=0)
        for sid in rec.samples.keys():
            target = sample_map[sid]
            if target in median_covs and (
                pd.isna(median_covs[target]) or median_covs[target] < min_cov
            ):
                rec.samples[sid]["GT"] = (None, None)
        outbcf.write(rec)


def main():
    inbcf = VariantFile(sys.argv[1], mode="r")
    outbcf = VariantFile(sys.argv[2], mode="w", header=inbcf.header)
    rdmat = RdMatReader(sys.argv[3])
    sample_map = read_sample_map(sys.argv[4])
    min_cov = int(sys.argv[5])
    filter_gt_by_coverage(inbcf, outbcf, rdmat, sample_map, min_cov)


if __name__ == "__main__":
    main()
