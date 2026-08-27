#!/usr/bin/env python3
"""Build a BED of dbNSFP sites eligible for at least one missense tier.

The four rank-score thresholds are the same values used by
``scripts/postprocess/tier_variants.py``. Regions may be padded to tolerate
representation changes before normalization. An optional existing BED can be
unioned into the result, which is how LoF-qualified exons are retained.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import duckdb

from build_target_bed import chromosome_key, merged


THRESHOLDS = {
    "ClinPred_rankscore": 0.4298,
    "AlphaMissense_rankscore": 0.9603,
    "popEVE_converted_rankscore": 0.9209,
    "MPC_rankscore": 0.8947,
}


def read_bed(path: Path) -> dict[str, list[tuple[int, int]]]:
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected at least 3 BED columns")
            intervals[fields[0]].append((int(fields[1]), int(fields[2])))
    return intervals


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dbnsfp", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--union-bed", type=Path)
    parser.add_argument("--padding", type=int, default=8)
    parser.add_argument("--add-chr-prefix", action="store_true")
    args = parser.parse_args()
    if args.padding < 0:
        parser.error("--padding must be non-negative")

    intervals = read_bed(args.union_bed) if args.union_bed else defaultdict(list)
    flags = " + ".join(
        f'CAST(COALESCE(TRY_CAST("{column}" AS DOUBLE) >= {threshold}, FALSE) AS INTEGER)'
        for column, threshold in THRESHOLDS.items()
    )
    query = f"""
        SELECT DISTINCT "#chr", TRY_CAST("pos(1-based)" AS BIGINT), ref
        FROM read_parquet(?)
        WHERE ({flags}) >= 1
          AND TRY_CAST("pos(1-based)" AS BIGINT) IS NOT NULL
        ORDER BY 1, 2
    """

    connection = duckdb.connect()
    candidates = connection.execute(query, [str(args.dbnsfp)]).fetchall()
    for source_chrom, position, ref in candidates:
        chrom = str(source_chrom)
        if args.add_chr_prefix and not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        reference_length = max(1, len(ref or ""))
        start = max(0, position - 1 - args.padding)
        end = position - 1 + reference_length + args.padding
        intervals[chrom].append((start, end))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    interval_count = 0
    covered_bases = 0
    with args.output.open("w") as handle:
        for chrom in sorted(intervals, key=chromosome_key):
            for start, end in merged(intervals[chrom]):
                handle.write(f"{chrom}\t{start}\t{end}\n")
                interval_count += 1
                covered_bases += end - start

    print(
        f"candidate_sites={len(candidates)} merged_intervals={interval_count} "
        f"covered_bases={covered_bases} padding={args.padding}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
