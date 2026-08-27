#!/usr/bin/env python3
"""Stream every non-empty ALT allele in a VCZ store to Parquet."""

import argparse
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import zarr


SCHEMA = pa.schema([
    ("variant_index", pa.int64()),
    ("alt_index", pa.int16()),
    ("chrom", pa.string()),
    ("pos", pa.int32()),
    ("ref", pa.string()),
    ("alt", pa.string()),
    ("record_id", pa.string()),
])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--zarr", required=True, type=Path)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    root = zarr.open_group(str(args.zarr), mode="r")
    alleles = root["variant_allele"]
    positions = root["variant_position"]
    contigs = root["variant_contig"]
    contig_ids = [str(value) for value in root["contig_id"][:]]
    try:
        contig_index = contig_ids.index(args.chrom)
    except ValueError as exc:
        raise ValueError("contig not found: {}".format(args.chrom)) from exc

    args.output.parent.mkdir(parents=True, exist_ok=True)
    writer = pq.ParquetWriter(args.output, SCHEMA, compression="zstd")
    total = 0
    chunk_size = alleles.chunks[0]
    try:
        for start in range(0, alleles.shape[0], chunk_size):
            end = min(start + chunk_size, alleles.shape[0])
            allele_chunk = alleles[start:end, :]
            position_chunk = positions[start:end]
            contig_chunk = contigs[start:end]
            rows = {name: [] for name in SCHEMA.names}
            for local_index in range(end - start):
                if int(contig_chunk[local_index]) != contig_index:
                    continue
                ref = str(allele_chunk[local_index, 0])
                variant_index = start + local_index
                for alt_index in range(1, allele_chunk.shape[1]):
                    alt = str(allele_chunk[local_index, alt_index])
                    if not alt:
                        continue
                    rows["variant_index"].append(variant_index)
                    rows["alt_index"].append(alt_index)
                    rows["chrom"].append(args.chrom)
                    rows["pos"].append(int(position_chunk[local_index]))
                    rows["ref"].append(ref)
                    rows["alt"].append(alt)
                    rows["record_id"].append("z{}_a{}".format(variant_index, alt_index))
            if rows["variant_index"]:
                table = pa.Table.from_pydict(rows, schema=SCHEMA)
                writer.write_table(table)
                total += table.num_rows
    finally:
        writer.close()
    print("alleles={:,} output={}".format(total, args.output))


if __name__ == "__main__":
    main()
