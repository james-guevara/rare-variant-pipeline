#!/usr/bin/env python3
"""Add cohort AC, AN, and AF to an allele sidecar using VCZ genotypes."""

import argparse
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import zarr


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--zarr", required=True, type=Path)
    parser.add_argument("--alleles", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    table = pq.read_table(args.alleles)
    variant_indexes = table["variant_index"].to_numpy()
    alt_indexes = table["alt_index"].to_numpy()
    order = np.argsort(variant_indexes, kind="stable")

    root = zarr.open_group(str(args.zarr), mode="r")
    genotypes = root["call_genotype"]
    masks = root["call_genotype_mask"]
    chunk_size = genotypes.chunks[0]

    ac = np.zeros(len(table), dtype=np.int32)
    an = np.zeros(len(table), dtype=np.int32)
    sorted_variants = variant_indexes[order]
    for chunk_id in np.unique(sorted_variants // chunk_size):
        chunk_start = int(chunk_id * chunk_size)
        chunk_end = min(chunk_start + chunk_size, genotypes.shape[0])
        rows = order[(sorted_variants // chunk_size) == chunk_id]
        local_variants = variant_indexes[rows] - chunk_start

        gt_chunk = genotypes[chunk_start:chunk_end, :, :]
        mask_chunk = masks[chunk_start:chunk_end, :, :]
        for row, local_variant in zip(rows, local_variants):
            gt = gt_chunk[local_variant]
            called = (~mask_chunk[local_variant]) & (gt >= 0)
            an[row] = int(called.sum())
            ac[row] = int(((gt == alt_indexes[row]) & called).sum())

    af = np.divide(
        ac,
        an,
        out=np.full(len(table), np.nan, dtype=np.float64),
        where=an > 0,
    )
    result = table.append_column("cohort_ac", pa.array(ac))
    result = result.append_column("cohort_an", pa.array(an))
    result = result.append_column("cohort_af", pa.array(af))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(result, args.output, compression="zstd")

    print(f"alleles={len(table):,} output={args.output}")
    for threshold in (0.01, 0.001, 0.0001):
        kept = int(np.count_nonzero(af < threshold))
        print(
            f"AF < {threshold:g}: {kept:,}/{len(table):,} "
            f"({kept / len(table):.1%})"
        )
    for maximum_ac in (0, 1, 2, 5, 10):
        kept = int(np.count_nonzero(ac <= maximum_ac))
        print(
            f"AC <= {maximum_ac}: {kept:,}/{len(table):,} "
            f"({kept / len(table):.1%})"
        )


if __name__ == "__main__":
    main()
