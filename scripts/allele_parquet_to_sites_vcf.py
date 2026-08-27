#!/usr/bin/env python3
"""Write a biallelic, sites-only VCF from a Zarr allele-pointer Parquet."""

import argparse
from pathlib import Path

import pyarrow.parquet as pq


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--output-contig", required=True)
    parser.add_argument("--contig-length", required=True, type=int)
    args = parser.parse_args()

    table = pq.read_table(
        args.input,
        columns=["variant_index", "alt_index", "pos", "ref", "alt"],
    ).sort_by([("pos", "ascending"), ("variant_index", "ascending"), ("alt_index", "ascending")])
    data = table.to_pydict()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as out:
        out.write("##fileformat=VCFv4.3\n")
        out.write("##source=rare-variant-pipeline-zarr-allele-extractor\n")
        out.write(f"##contig=<ID={args.output_contig},length={args.contig_length}>\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for variant_index, alt_index, pos, ref, alt in zip(
            data["variant_index"],
            data["alt_index"],
            data["pos"],
            data["ref"],
            data["alt"],
        ):
            pointer = f"z{variant_index}_a{alt_index}"
            out.write(
                f"{args.output_contig}\t{pos}\t{pointer}\t{ref}\t{alt}"
                "\t.\tPASS\t.\n"
            )

    print(f"rows={table.num_rows:,} output={args.output}")


if __name__ == "__main__":
    main()
