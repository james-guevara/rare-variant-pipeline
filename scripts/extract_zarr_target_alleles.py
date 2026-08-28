#!/usr/bin/env python3
"""Extract one row per ALT allele from a VCZ store within target BED intervals.

The output retains ``variant_index`` and ``alt_index`` so normalized and
annotated alleles can be mapped back to the unchanged genotype arrays.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
from functools import partial
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
    if not intervals or positions.size == 0:
        return np.array([], dtype=np.int64)
    positions64 = positions.astype(np.int64, copy=False)
    lengths64 = lengths.astype(np.int64, copy=False)
    bed_starts = np.fromiter((start for start, _ in intervals), dtype=np.int64)
    bed_ends = np.fromiter((end for _, end in intervals), dtype=np.int64)
    max_length = int(lengths.max(initial=1))
    # Vectorizing these bounds avoids two full NumPy dispatches for every BED
    # interval. VCF POS is 1-based; BED is zero-based, half-open. Include records
    # that begin before an interval but whose REF span overlaps it.
    lows = np.searchsorted(
        positions64, bed_starts - max_length + 2, side="left"
    )
    highs = np.searchsorted(positions64, bed_ends, side="right")
    selected = []
    for bed_start, bed_end, lo, hi in zip(
        bed_starts, bed_ends, lows, highs, strict=True
    ):
        indexes = np.arange(lo, hi, dtype=np.int64)
        record_start = positions64[lo:hi] - 1
        record_end = record_start + lengths64[lo:hi]
        selected.append(
            indexes[(record_start < bed_end) & (record_end > bed_start)]
        )
    if not selected:
        return np.array([], dtype=np.int64)
    return np.unique(np.concatenate(selected))


ROW_NAMES = (
    "variant_index", "alt_index", "chrom", "pos", "record_end",
    "ref", "alt", "variant_id", "filter",
)


def selected_chunk_rows(
    chunk_id: int,
    indexes: np.ndarray,
    variant_chunk: int,
    allele_array,
    id_array,
    filter_array,
    filter_ids: list[str],
    positions: np.ndarray,
    lengths: np.ndarray,
    chrom: str,
) -> dict[str, list]:
    """Read only selected rows from one outer Zarr shard."""
    chunk_start = int(chunk_id * variant_chunk)
    chunk_end = min(chunk_start + variant_chunk, len(positions))
    in_chunk = indexes[(indexes >= chunk_start) & (indexes < chunk_end)]
    alleles = allele_array.get_orthogonal_selection((in_chunk, slice(None)))
    ids = id_array.get_orthogonal_selection((in_chunk,))
    filters = filter_array.get_orthogonal_selection((in_chunk, slice(None)))
    rows = {name: [] for name in ROW_NAMES}
    for row_index, global_index in enumerate(in_chunk.tolist()):
        ref = str(alleles[row_index, 0])
        filt = ";".join(
            filter_ids[i]
            for i, present in enumerate(filters[row_index])
            if present
        ) or "."
        for alt_index, alt in enumerate(alleles[row_index, 1:], start=1):
            alt = str(alt)
            if not alt:
                continue
            rows["variant_index"].append(global_index)
            rows["alt_index"].append(alt_index)
            rows["chrom"].append(chrom)
            rows["pos"].append(int(positions[global_index]))
            rows["record_end"].append(
                int(positions[global_index] - 1 + lengths[global_index])
            )
            rows["ref"].append(ref)
            rows["alt"].append(alt)
            rows["variant_id"].append(str(ids[row_index]))
            rows["filter"].append(filt)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--zarr", required=True, type=Path)
    parser.add_argument("--bed", required=True, type=Path)
    parser.add_argument("--chrom", default="chr22")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--workers", type=int, default=1,
        help="number of outer Zarr shards to read concurrently (default: 1)",
    )
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be at least 1")

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

    rows = {name: [] for name in ROW_NAMES}
    chunk_ids = np.unique(indexes // variant_chunk).tolist()
    read_chunk = partial(
        selected_chunk_rows,
        indexes=indexes,
        variant_chunk=variant_chunk,
        allele_array=allele_array,
        id_array=id_array,
        filter_array=filter_array,
        filter_ids=filter_ids,
        positions=positions,
        lengths=lengths,
        chrom=args.chrom,
    )
    if args.workers == 1:
        chunk_results = map(read_chunk, chunk_ids)
        for chunk_rows in chunk_results:
            for name in ROW_NAMES:
                rows[name].extend(chunk_rows[name])
    else:
        # Executor.map preserves chunk order, so worker count cannot alter row order.
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            for chunk_rows in executor.map(read_chunk, chunk_ids):
                for name in ROW_NAMES:
                    rows[name].extend(chunk_rows[name])

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
