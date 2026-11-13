#!/usr/bin/env python3
"""Compare standalone LOFTEE against a transcript-level Perl oracle TSV."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

from standalone_loftee import classify
from standalone_loftee.context import build_context
from standalone_loftee.resources import LofteeResources
from standalone_loftee.transcripts import TranscriptStore


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--oracle", type=Path, required=True)
    parser.add_argument("--transcripts", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--ancestor", type=Path, required=True)
    parser.add_argument("--gerp", type=Path, required=True)
    parser.add_argument("--conservation", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary-json", type=Path)
    args = parser.parse_args()

    mismatch_counts: Counter[str] = Counter()
    total = 0
    with (
        args.oracle.open(newline="") as source,
        args.output.open("w", newline="") as target,
        TranscriptStore(args.transcripts) as transcripts,
        LofteeResources(args.reference, args.ancestor, args.gerp, args.conservation) as resources,
    ):
        reader = csv.DictReader(source, delimiter="\t")
        output_fields = list(reader.fieldnames or ()) + [
            "python_LoF", "python_LoF_filter", "python_LoF_flags", "python_LoF_info"
        ]
        writer = csv.DictWriter(target, fieldnames=output_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            total += 1
            transcript = transcripts.get(row["Feature"])
            if transcript is None:
                raise ValueError(f"missing transcript: {row['Feature']}")
            result = classify(build_context(row, transcript, resources))
            predicted = {
                "python_LoF": result.lof,
                "python_LoF_filter": result.lof_filter or ".",
                "python_LoF_flags": result.lof_flags or ".",
                "python_LoF_info": result.lof_info or ".",
            }
            for oracle_field, python_field in (
                ("LoF", "python_LoF"),
                ("LoF_filter", "python_LoF_filter"),
                ("LoF_flags", "python_LoF_flags"),
                ("LoF_info", "python_LoF_info"),
            ):
                # VEP escapes commas inside CSQ subfields as ampersands. The
                # standalone tabular representation retains LOFTEE's commas.
                oracle_value = row[oracle_field].replace("&", ",")
                if oracle_value != predicted[python_field]:
                    mismatch_counts[oracle_field] += 1
            writer.writerow(row | predicted)

    print(f"rows\t{total}")
    for field in ("LoF", "LoF_filter", "LoF_flags", "LoF_info"):
        print(f"{field}_mismatches\t{mismatch_counts[field]}")
    summary_path = args.summary_json or args.output.with_suffix(".summary.json")
    summary_path.write_text(
        json.dumps(
            {
                "rows": total,
                "mismatches": {
                    field: mismatch_counts[field]
                    for field in ("LoF", "LoF_filter", "LoF_flags", "LoF_info")
                },
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    return 1 if mismatch_counts else 0


if __name__ == "__main__":
    raise SystemExit(main())
