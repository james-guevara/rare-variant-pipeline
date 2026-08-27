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
        for row in reader:
            transcript = transcripts.get(row["Feature"])
            if transcript is None:
                raise ValueError(f"missing transcript: {row['Feature']}")
            result = classify(build_context(row, transcript, resources))
            row.update({
                "LoF": result.lof,
                "LoF_filter": result.lof_filter or ".",
                "LoF_flags": result.lof_flags or ".",
                "LoF_info": result.lof_info or ".",
            })
            writer.writerow(row)


if __name__ == "__main__":
    main()
