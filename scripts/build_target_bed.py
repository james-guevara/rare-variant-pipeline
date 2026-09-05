#!/usr/bin/env python3
"""Build a merged BED for early targeted rare-variant extraction.

GeneBayes uses Ensembl gene IDs (``ensg``), which are also present as ``gene_id``
in Ensembl/GENCODE GTF files. Selecting on that stable key avoids gene-symbol
aliases. GTF coordinates are one-based inclusive; BED coordinates are zero-based,
half-open.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable, TextIO


ATTRIBUTE_RE = re.compile(r'(\S+)\s+"([^"]*)"')


def open_text(path: Path) -> TextIO:
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def unversioned(identifier: str) -> str:
    return identifier.split(".", 1)[0]


def selected_genes(path: Path, threshold: float) -> set[str]:
    with open_text(path) as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        required = {"ensg", "post_mean"}
        missing = required.difference(rows.fieldnames or ())
        if missing:
            raise ValueError(f"GeneBayes table lacks columns: {', '.join(sorted(missing))}")
        selected = set()
        for row in rows:
            try:
                score = float(row["post_mean"])
            except (TypeError, ValueError):
                continue
            if score >= threshold:
                selected.add(unversioned(row["ensg"]))
        return selected


def gtf_intervals(
    path: Path,
    genes: set[str],
    features: set[str],
    padding: int,
    chromosomes: set[str] | None,
    add_chr_prefix: bool,
) -> tuple[dict[str, list[tuple[int, int]]], set[str]]:
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    observed: set[str] = set()
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in features:
                continue
            source_chrom = fields[0]
            if chromosomes is not None and source_chrom not in chromosomes:
                continue
            chrom = (
                f"chr{source_chrom}"
                if add_chr_prefix and not source_chrom.startswith("chr")
                else source_chrom
            )
            attrs = dict(ATTRIBUTE_RE.findall(fields[8]))
            gene_id = unversioned(attrs.get("gene_id", ""))
            if gene_id not in genes:
                continue
            start_1 = int(fields[3])
            end_1 = int(fields[4])
            start_0 = max(0, start_1 - 1 - padding)
            end_0 = end_1 + padding
            intervals[chrom].append((start_0, end_0))
            observed.add(gene_id)
    return intervals, observed


def merged(intervals: Iterable[tuple[int, int]]) -> Iterable[tuple[int, int]]:
    ordered = sorted(intervals)
    if not ordered:
        return
    current_start, current_end = ordered[0]
    for start, end in ordered[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            yield current_start, current_end
            current_start, current_end = start, end
    yield current_start, current_end


def chromosome_key(chrom: str) -> tuple[int, int | str]:
    bare = chrom.removeprefix("chr")
    if bare.isdigit():
        return (0, int(bare))
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (0, special[bare]) if bare in special else (1, chrom)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genebayes", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-post-mean", type=float, default=0.03)
    parser.add_argument(
        "--gene-ids",
        type=Path,
        help="optional newline-delimited Ensembl gene allowlist",
    )
    parser.add_argument("--features", default="exon",
                        help="comma-separated GTF features (default: exon)")
    parser.add_argument("--padding", type=int, default=8,
                        help="bases added to both sides of every interval")
    parser.add_argument("--chroms", default="",
                        help="optional comma-separated chromosome allowlist")
    parser.add_argument(
        "--add-chr-prefix",
        action="store_true",
        help="prefix GTF sequence names with 'chr' in the output BED",
    )
    args = parser.parse_args()

    if args.padding < 0:
        parser.error("--padding must be non-negative")
    features = {value.strip() for value in args.features.split(",") if value.strip()}
    if not features:
        parser.error("--features must contain at least one feature")
    chromosomes = {c.strip() for c in args.chroms.split(",") if c.strip()} or None

    genes = selected_genes(args.genebayes, args.min_post_mean)
    if args.gene_ids:
        allowed = {
            unversioned(line.strip())
            for line in args.gene_ids.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        }
        genes.intersection_update(allowed)
    if not genes:
        raise SystemExit("No genes met the GeneBayes threshold")
    intervals, observed = gtf_intervals(
        args.gtf, genes, features, args.padding, chromosomes, args.add_chr_prefix
    )
    if not intervals:
        raise SystemExit("No selected genes matched the requested GTF features/chromosomes")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    n_intervals = 0
    covered_bases = 0
    with args.output.open("w") as out:
        for chrom in sorted(intervals, key=chromosome_key):
            for start, end in merged(intervals[chrom]):
                out.write(f"{chrom}\t{start}\t{end}\n")
                n_intervals += 1
                covered_bases += end - start

    missing = genes.difference(observed)
    print(
        f"selected_genes={len(genes)} matched_genes={len(observed)} "
        f"missing_genes={len(missing)} merged_intervals={n_intervals} "
        f"covered_bases={covered_bases}",
        file=sys.stderr,
    )
    if missing:
        preview = ",".join(sorted(missing)[:10])
        print(f"missing_gene_preview={preview}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
