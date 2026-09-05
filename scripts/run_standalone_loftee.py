#!/usr/bin/env python3
"""Run standalone LOFTEE on base-annotation rows without a Perl oracle."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from standalone_loftee import classify
from standalone_loftee.context import build_context
from standalone_loftee.resources import LofteeResources
from standalone_loftee.transcripts import TranscriptStore


LOF_CONSEQUENCES = frozenset({
    "stop_gained",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
})


def transcript_is_absent(value: str) -> bool:
    return value.strip() in {"", ".", "-"}


def is_lof_candidate(row: dict[str, str]) -> bool:
    # An absent BIOTYPE must fall through to the transcript model, matching
    # build_context's existing fallback. A present non-coding type can skip.
    biotype = row.get("BIOTYPE", "")
    if biotype and biotype != "protein_coding":
        return False
    return bool(LOF_CONSEQUENCES.intersection(row.get("Consequence", "").split("&")))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--transcripts", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--ancestor", required=True, type=Path)
    parser.add_argument("--gerp", required=True, type=Path)
    parser.add_argument("--conservation", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    with (
        args.input.open(newline="") as source,
        args.output.open("w", newline="") as target,
        TranscriptStore(args.transcripts) as transcripts,
        LofteeResources(args.reference, args.ancestor, args.gerp, args.conservation) as resources,
    ):
        reader = csv.DictReader(source, delimiter="\t")
        fields = list(reader.fieldnames or ()) + ["LoF", "LoF_filter", "LoF_flags", "LoF_info"]
        writer = csv.DictWriter(target, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        scanned = candidates = classified = 0
        for row in reader:
            scanned += 1
            if not is_lof_candidate(row):
                continue
            candidates += 1
            if transcript_is_absent(row.get("Feature", "")):
                continue
            transcript = transcripts.get(row["Feature"])
            if transcript is None:
                raise ValueError(f"missing transcript: {row['Feature']}")
            result = classify(build_context(row, transcript, resources))
            if not result.lof:
                continue
            row.update({
                "LoF": result.lof,
                "LoF_filter": result.lof_filter or ".",
                "LoF_flags": result.lof_flags or ".",
                "LoF_info": result.lof_info or ".",
            })
            writer.writerow(row)
            if result.lof:
                classified += 1
        print(
            f"rows_scanned={scanned:,} lof_candidates={candidates:,} "
            f"lof_classified={classified:,}"
        )


if __name__ == "__main__":
    main()
