#!/usr/bin/env python3
"""Package chromosome-level targeted outputs as portable Nextflow artifacts."""

import argparse
import hashlib
import json
import shutil
from pathlib import Path


SCIENTIFIC_OUTPUTS = {
    "plof_counts": "12.plof-per-sample-counts.tsv",
    "missense_counts": "12.missense-per-sample-counts.tsv",
    "plof_primary_burden": "12.plof-primary-sample-burden.tsv",
    "plof_sensitivity_burden": "12.plof-sensitivity-sample-burden.tsv",
}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--science-manifest", required=True, type=Path)
    parser.add_argument("--bindings", required=True, type=Path)
    parser.add_argument("--run-root", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()

    science = json.loads(args.science_manifest.read_text())
    bindings = json.loads(args.bindings.read_text())
    run_root = args.run_root or Path(bindings["run_root"])
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    packaged = {}
    if not args.preflight_only:
        success = run_root / "_SUCCESS"
        if not success.is_file():
            raise ValueError(f"completed run is missing _SUCCESS: {run_root}")
        for logical_name, filename in SCIENTIFIC_OUTPUTS.items():
            source = run_root / filename
            if not source.is_file():
                continue
            destination = output_dir / filename
            shutil.copyfile(source, destination)
            packaged[logical_name] = {
                "file": filename,
                "sha256": file_sha256(destination),
                "size_bytes": destination.stat().st_size,
            }

        if "plof_counts" not in packaged:
            raise ValueError(f"completed run is missing required LoF counts: {run_root}")

    manifest = {
        "schema_version": 1,
        "cohort": science["cohort"],
        "chromosome": science["chromosome"],
        "status": "PREFLIGHT_ONLY" if args.preflight_only else "SUCCEEDED",
        "files": packaged,
    }
    (output_dir / "targeted-output-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
