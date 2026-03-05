import argparse
from pathlib import Path
import sys

import pandas as pd
from pysam import VariantFile
from pysam import TabixFile


class RdMatReader:
    def __init__(self, path, targets):
        self.handle = TabixFile(str(path))
        self.targets = targets
        # skip the first three fields '#Chr\tStart\tEnd'
        samples = self.handle.header[0].lstrip("#").split("\t")[3:]
        self.indexes = [i for i, v in enumerate(samples) if v in targets]
        self.header = [samples[i] for i in self.indexes]

    def fetch(self, contig, start, end):
        for rec in self.handle.fetch(contig, start, end):
            yield self._parse_rec(rec)

    def fetch_as_df(self, contig, start, end):
        return pd.DataFrame.from_records(
            self.fetch(contig, start, end), columns=self.header
        )

    def _parse_rec(self, rec):
        fields = [int(x) for x in rec.split("\t")[3:]]
        return tuple(fields[i] for i in self.indexes)


def read_sample_map(path):
    tmp = dict()
    with open(path, mode="r", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip().split("\t")
            tmp[parts[0]] = parts[1]

    return tmp


def filter_site_gts(record, sample_map, medians, mads, min_cov, max_mad):
    for sid in record.samples.keys():
        target = sample_map[sid]
        if target not in medians:
            continue
        if (record.info["SVTYPE"] != "DEL" and (pd.isna(medians[target]) or medians[target] < min_cov)):
            record.samples[sid]["GT"] = (None, None)
        elif mads[target] > max_mad:
            record.samples[sid]["GT"] = (None, None)

    return record

def mad(x):
    return (x - x.median()).abs().median()


def do_filtering(inbcf, outbcf, rdmat, sample_map, min_cov, max_mad):
    for rec in inbcf.fetch():
        cov = rdmat.fetch_as_df(rec.contig, rec.pos - 1, rec.stop)
        site_median_covs = cov.median(axis=0)
        mads = cov.apply(mad)
        filtered_rec = filter_site_gts(rec, sample_map, site_median_covs, mads, min_cov, max_mad)
        outbcf.write(rec)

def parse_args():
    parser = argparse.ArgumentParser(description = "Update genotypes in a VCF/BCF based on site coverage")
    parser.add_argument("in_bcf", type = Path, help = "VCF/BCF to update")
    parser.add_argument("out_bcf", type = Path, help = "Where to write the updated VCF/BCF")
    parser.add_argument("rd_mat", type = Path, help = "Binned read-depth matrix")
    parser.add_argument("targets", type = Path, help = "Two-column TSV with ID of sample in VCF/BCF in first column and ID of sample to check for read-depth")
    parser.add_argument("--min-cov", type = int, default = 10, help = "Minimum median coverage a sample must have to keep a genotype")
    parser.add_argument("--max-mad", type = float, default = 10, help = "Maximum MAD of coverage a sample can have to keep a genotype")

    return parser.parse_args()


def main():
    args = parse_args()
    inbcf = VariantFile(args.in_bcf, mode="r")
    outbcf = VariantFile(args.out_bcf, mode="w", header=inbcf.header)
    sample_map = read_sample_map(args.targets)
    targets = set(sample_map[x] for x in inbcf.header.samples if x in sample_map)
    rdmat = RdMatReader(args.rd_mat, targets)
    do_filtering(inbcf, outbcf, rdmat, sample_map, args.min_cov, args.max_mad)


if __name__ == "__main__":
    main()
