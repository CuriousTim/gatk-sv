import argparse
from pathlib import Path

from pysam import VariantFile


def is_valid_ped_sample(x):
    return len(x) > 0 and x != "0"


def read_pedigree(path):
    fams = dict()
    with open(path, mode="r") as fp:
        for line in fp:
            fields = line.rstrip().split("\t")
            if (
                is_valid_ped_sample(fields[1])
                and is_valid_ped_sample(fields[2])
                and is_valid_ped_sample(fields[3])
            ):
                fams[fields[1]] = (fields[2], fields[3])

    return fams


def null_inherited_gt(header, rec, ped):
    for sid in rec.samples.keys():
        if sid not in ped:
            continue
        if 1 not in rec.samples[sid]["GT"]:
            continue
        father, mother = ped[sid]
        if (
            father not in header.samples
            or mother not in header.samples
            or 1 in rec.samples[father]["GT"]
            or 1 in rec.samples[mother]["GT"]
        ):
            rec.samples[sid]["GT"] = (None, None)

        # low GQ CNVs are generally unreliable
        if (rec.info["SVTYPE"] == "DEL" or rec.info["SVTYPE"] == "DUP") and (
            rec.samples[sid]["GQ"] == 0
            or rec.samples[sid]["GQ"] == 0
            or rec.samples[sid]["GQ"] == 0
        ):
            rec.samples[sid]["GT"] = (None, None)

    return rec


def process_bcf(bcf, ped, filtered_bcf):
    for rec in bcf.fetch():
        null_inherited_gt(bcf.header, rec, ped)
        site_has_offspring_alt = False
        for sid in rec.samples.keys():
            if sid in ped and 1 in rec.samples[sid]["GT"]:
                site_has_offspring_alt = True
                break
        if site_has_offspring_alt:
            filtered_bcf.write(rec)


def parse_args():
    parser = argparse.ArgumentParser("Remove inherited variants from a GATK-SV BCF")
    parser.add_argument("bcf", type=Path, help="A GATK-SV cohort BCF")
    parser.add_argument("pedigree", type=Path, help="PED file")
    parser.add_argument("filtered_bcf", type=Path, help="Filtered BCF")

    return parser.parse_args()


def main():
    args = parse_args()
    bcf = VariantFile(args.bcf, mode="r")
    ped = read_pedigree(args.pedigree)
    filtered_bcf = VariantFile(args.filtered_bcf, mode="wb", header=bcf.header)

    process_bcf(bcf, ped, filtered_bcf)


if __name__ == "__main__":
    main()
