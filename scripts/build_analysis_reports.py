#!/usr/bin/env python3
"""Create deterministic missingness and provenance reports for an analysis table."""

import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--dictionary", required=True, type=Path)
    parser.add_argument("--qc", required=True, type=Path)
    parser.add_argument("--cohort-id", required=True)
    parser.add_argument("--pgs-mode", required=True, choices=("computed", "precomputed", "disabled"))
    parser.add_argument("--rare-mode", required=True, choices=("computed", "precomputed", "disabled"))
    parser.add_argument("--cnv-mode", required=True, choices=("computed", "precomputed", "disabled"))
    parser.add_argument("--missingness", required=True, type=Path)
    parser.add_argument("--run-manifest", required=True, type=Path)
    args = parser.parse_args()

    with args.dataset.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        rows = list(reader)
    with args.dictionary.open(newline="") as handle:
        definitions = list(csv.DictReader(handle, delimiter="\t"))
    if [row["variable"] for row in definitions] != fields:
        raise ValueError("analysis dictionary variables do not exactly match dataset columns")

    with args.missingness.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("variable", "missing_participants", "observed_participants"))
        for field in fields:
            missing = sum(row[field] == "" for row in rows)
            writer.writerow((field, missing, len(rows) - missing))

    analysis_qc = json.loads(args.qc.read_text())
    args.run_manifest.write_text(json.dumps({
        "schema_version": 1,
        "cohort_id": args.cohort_id,
        "participants": len(rows),
        "variables": len(fields),
        "component_modes": {"pgs": args.pgs_mode, "rare": args.rare_mode, "cnv": args.cnv_mode},
        "missing_policies": {
            name: component["missing_policy"]
            for name, component in analysis_qc.get("components", {}).items()
        },
        "artifacts": {
            "analysis_dataset.tsv": sha256(args.dataset),
            "analysis_dataset_dictionary.tsv": sha256(args.dictionary),
            "analysis_qc.json": sha256(args.qc),
        },
    }, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
