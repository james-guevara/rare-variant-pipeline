#!/usr/bin/env python3
"""Build a merged BED containing all protein-coding exons in a GTF."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

from build_target_bed import ATTRIBUTE_RE, chromosome_key, merged, open_text


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--chroms", default="")
    parser.add_argument("--padding", type=int, default=8)
    parser.add_argument("--add-chr-prefix", action="store_true")
    args = parser.parse_args()
    if args.padding < 0:
        parser.error("--padding must be non-negative")
    chromosomes = {value.strip() for value in args.chroms.split(",") if value.strip()} or None

    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    genes: set[str] = set()
    with open_text(args.gtf) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            if chromosomes is not None and fields[0] not in chromosomes:
                continue
            attributes = dict(ATTRIBUTE_RE.findall(fields[8]))
            biotype = attributes.get("gene_biotype", attributes.get("gene_type", ""))
            if biotype != "protein_coding":
                continue
            chrom = fields[0]
            if args.add_chr_prefix and not chrom.startswith("chr"):
                chrom = f"chr{chrom}"
            start = max(0, int(fields[3]) - 1 - args.padding)
            end = int(fields[4]) + args.padding
            intervals[chrom].append((start, end))
            genes.add(attributes.get("gene_id", ""))

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
        f"protein_coding_genes={len(genes)} merged_intervals={interval_count} "
        f"covered_bases={covered_bases} padding={args.padding}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
