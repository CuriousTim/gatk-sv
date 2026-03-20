import argparse
import gzip
from pathlib import Path

from pysam import VariantFile

OUTPUT_HEADER = "\t".join(
    [
        "chr",
        "start",
        "end",
        "svlen",
        "alt",
        "vid",
        "svtype",
        "af",
        "ac",
        "algorithm",
        "evidence",
        "sample",
        "o_gt",
        "o_ev",
        "o_gq",
        "f_gt",
        "f_ev",
        "f_gq",
        "m_gt",
        "m_ev",
        "m_gq",
    ]
)


def is_valid_ped_sample(x):
    return len(x) > 0 and x != "0"


def tuple_to_str(x, sep=","):
    return sep.join(x)


def gt_to_str(x):
    return "/".join("." if allele is None else str(allele) for allele in x)


def tsv_reader(fp):
    for line in fp:
        yield line.rstrip("\n").split("\t")


def read_pedigree(path):
    fams = dict()
    with open(path, mode="r") as fp:
        for fields in tsv_reader(fp):
            if (
                is_valid_ped_sample(fields[1])
                and is_valid_ped_sample(fields[2])
                and is_valid_ped_sample(fields[3])
            ):
                fams[fields[1]] = (fields[2], fields[3])

    return fams


def read_denovos(path):
    denovos = dict()
    with gzip.open(path, mode="rt") as fp:
        for vid, sid in tsv_reader(fp):
            carriers = denovos.get(vid, set())
            carriers.add(sid)
            denovos[vid] = carriers

    return denovos


def write_site_annotations(fp, bcf_header, rec, carriers, ped):
    for sid in carriers:
        father, mother = ped[sid]
        fp.write(f"{rec.chrom}\t{rec.start}\t{rec.stop}\t{rec.info['SVLEN']}\t")
        fp.write(f"{tuple_to_str(rec.alts)}\t{rec.id}\t{rec.info["SVTYPE"]}\t")
        fp.write(f"{rec.info["AF"][0]:.5f}\t{rec.info["AC"][0]}\t")
        fp.write(f"{tuple_to_str(rec.info['ALGORITHMS'])}\t")
        fp.write(f"{tuple_to_str(rec.info['EVIDENCE'])}\t")
        fp.write(f"{sid}\t")
        fp.write(f"{gt_to_str(rec.samples[sid]['GT'])}\t")
        fp.write(f"{tuple_to_str(rec.samples[sid]['EV'])}\t")
        fp.write(f"{rec.samples[sid]['GQ']}\t")
        fp.write(f"{gt_to_str(rec.samples[father]['GT'])}\t")
        fp.write(f"{tuple_to_str(rec.samples[father]['EV'])}\t")
        fp.write(f"{rec.samples[father]['GQ']}\t")
        fp.write(f"{gt_to_str(rec.samples[mother]['GT'])}\t")
        fp.write(f"{tuple_to_str(rec.samples[mother]['EV'])}\t")
        fp.write(f"{rec.samples[mother]['GQ']}\n")


def annotate(denovos, bcf, ped, output_path):
    with gzip.open(output_path, "wt") as fp:
        fp.write(OUTPUT_HEADER)
        fp.write("\n")
        for rec in bcf.fetch():
            if rec.id in denovos:
                write_site_annotations(fp, bcf.header, rec, denovos[rec.id], ped)


def parse_args():
    parser = argparse.ArgumentParser("Annotate de novo SVs with a BCF")
    parser.add_argument(
        "denovos", type=Path, help="TSV of de novos. Variant ID and sample ID."
    )
    parser.add_argument("pedigree", type=Path, help="PED file")
    parser.add_argument("bcf", type=Path, help="A BCF to use for annotations.")
    parser.add_argument(
        "output", type=Path, help="Where to write the annotated de novos."
    )

    return parser.parse_args()


def main():
    args = parse_args()
    denovos = read_denovos(args.denovos)
    ped = read_pedigree(args.pedigree)
    bcf = VariantFile(args.bcf, mode="r")
    output_path = args.output

    annotate(denovos, bcf, ped, output_path)


if __name__ == "__main__":
    main()
