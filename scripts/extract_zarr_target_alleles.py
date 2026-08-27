#!/usr/bin/env python3
"""Extract one row per ALT allele from a VCZ store within target BED intervals.

The output retains ``variant_index`` and ``alt_index`` so normalized and
annotated alleles can be mapped back to the unchanged genotype arrays.
"""

import argparse
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import zarr


def load_intervals(path: Path, chrom: str) -> list[tuple[int, int]]:
    intervals = []
    with path.open() as bed:
        for line in bed:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip().split("\t")
            if fields[0] == chrom:
                intervals.append((int(fields[1]), int(fields[2])))
    intervals.sort()
    merged: list[list[int]] = []
    for start, end in intervals:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def overlapping_indexes(
    positions: np.ndarray,
    lengths: np.ndarray,
    intervals: list[tuple[int, int]],
) -> np.ndarray:
    selected = []
    max_length = int(lengths.max(initial=1))
    for bed_start, bed_end in intervals:
        # VCF POS is 1-based; BED is zero-based, half-open. Include records that
        # begin before an interval but whose REF span overlaps it.
        lo = np.searchsorted(positions, bed_start - max_length + 2, side="left")
        hi = np.searchsorted(positions, bed_end, side="right")
        indexes = np.arange(lo, hi, dtype=np.int64)
        record_start = positions[lo:hi].astype(np.int64) - 1
        record_end = record_start + lengths[lo:hi].astype(np.int64)
        selected.extend(indexes[(record_start < bed_end) & (record_end > bed_start)])
    if not selected:
        return np.array([], dtype=np.int64)
    return np.unique(np.asarray(selected, dtype=np.int64))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--zarr", required=True, type=Path)
    parser.add_argument("--bed", required=True, type=Path)
    parser.add_argument("--chrom", default="chr22")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    root = zarr.open_group(str(args.zarr), mode="r")
    contig_ids = root["contig_id"][:].tolist()
    try:
        contig_index = contig_ids.index(args.chrom)
    except ValueError as exc:
        raise SystemExit(f"contig {args.chrom!r} is absent from the VCZ store") from exc

    positions = root["variant_position"][:]
    lengths = root["variant_length"][:]
    contigs = root["variant_contig"][:]
    intervals = load_intervals(args.bed, args.chrom)
    indexes = overlapping_indexes(positions, lengths, intervals)
    indexes = indexes[contigs[indexes] == contig_index]

    allele_array = root["variant_allele"]
    id_array = root["variant_id"]
    filter_array = root["variant_filter"]
    filter_ids = root["filter_id"][:].tolist()
    variant_chunk = allele_array.chunks[0]

    rows = {name: [] for name in (
        "variant_index", "alt_index", "chrom", "pos", "record_end",
        "ref", "alt", "variant_id", "filter",
    )}
    for chunk_id in np.unique(indexes // variant_chunk):
        chunk_start = int(chunk_id * variant_chunk)
        chunk_end = min(chunk_start + variant_chunk, len(positions))
        in_chunk = indexes[(indexes >= chunk_start) & (indexes < chunk_end)]
        local = in_chunk - chunk_start
        alleles = allele_array[chunk_start:chunk_end, :]
        ids = id_array[chunk_start:chunk_end]
        filters = filter_array[chunk_start:chunk_end, :]
        for global_index, local_index in zip(in_chunk.tolist(), local.tolist()):
            ref = str(alleles[local_index, 0])
            filt = ";".join(
                filter_ids[i] for i, present in enumerate(filters[local_index]) if present
            ) or "."
            for alt_index, alt in enumerate(alleles[local_index, 1:], start=1):
                alt = str(alt)
                if not alt:
                    continue
                rows["variant_index"].append(global_index)
                rows["alt_index"].append(alt_index)
                rows["chrom"].append(args.chrom)
                rows["pos"].append(int(positions[global_index]))
                rows["record_end"].append(
                    int(positions[global_index] - 1 + lengths[global_index])
                )
                rows["ref"].append(ref)
                rows["alt"].append(alt)
                rows["variant_id"].append(str(ids[local_index]))
                rows["filter"].append(filt)

    schema = pa.schema([
        ("variant_index", pa.int64()),
        ("alt_index", pa.int16()),
        ("chrom", pa.string()),
        ("pos", pa.int32()),
        ("record_end", pa.int32()),
        ("ref", pa.string()),
        ("alt", pa.string()),
        ("variant_id", pa.string()),
        ("filter", pa.string()),
    ])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table = pa.Table.from_pydict(rows, schema=schema)
    pq.write_table(table, args.output, compression="zstd")
    print(
        f"intervals={len(intervals):,} records={len(indexes):,} "
        f"alleles={table.num_rows:,} output={args.output}"
    )


if __name__ == "__main__":
    main()
