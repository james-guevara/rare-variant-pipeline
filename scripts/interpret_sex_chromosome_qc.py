#!/usr/bin/env python3
"""Add non-diagnostic copy-number pattern flags to a ploidy audit TSV."""

import argparse
import csv
import json
from pathlib import Path


def between(value: float, lower: float, upper: float) -> bool:
    return lower <= value <= upper


def interpret(row: dict[str, str], thresholds: dict[str, float]) -> str:
    """Return a conservative evidence pattern, never a clinical karyotype."""
    x_ratio = float(row["x_autosome_dp_ratio"])
    y_ratio = float(row["y_autosome_dp_ratio"])
    y_calls = float(row["y_call_rate"])
    x_one = between(
        x_ratio,
        thresholds["one_x_dp_ratio_min"],
        thresholds["one_x_dp_ratio_max"],
    )
    x_two = between(
        x_ratio,
        thresholds["two_x_dp_ratio_min"],
        thresholds["two_x_dp_ratio_max"],
    )
    y_present = y_calls >= thresholds["y_present_call_rate_min"]
    y_absent = y_calls < thresholds["y_absent_call_rate_max"]

    if x_two and y_present:
        return "two-X-plus-Y-compatible"
    if x_one and y_present and y_ratio >= thresholds["excess_y_dp_ratio_min"]:
        return "one-X-with-excess-Y-signal"
    if x_one and y_present:
        return "one-X-plus-Y-with-GT-ploidy-discordance"
    if x_two and y_absent:
        return "two-X-no-Y-compatible"
    if x_two:
        return "two-X-with-uncertain-Y-signal"
    return "unresolved"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    config = json.loads(args.config.read_text())
    thresholds = config["aneuploidy_triage_thresholds"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.audit.open() as source, args.output.open("w", newline="") as target:
        reader = csv.DictReader(source, delimiter="\t")
        required = {
            "sample_id", "inferred_karyotype", "y_call_rate",
            "x_autosome_dp_ratio", "y_autosome_dp_ratio",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"audit is missing required columns: {sorted(missing)}")
        fieldnames = list(reader.fieldnames or []) + [
            "copy_number_evidence_pattern", "sex_chromosome_analysis_eligible",
        ]
        writer = csv.DictWriter(target, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in reader:
            row["copy_number_evidence_pattern"] = interpret(row, thresholds)
            row["sex_chromosome_analysis_eligible"] = (
                "1" if row["inferred_karyotype"] in {"XX-like", "XY-like"} else "0"
            )
            writer.writerow(row)


if __name__ == "__main__":
    main()
