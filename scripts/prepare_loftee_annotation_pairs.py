#!/usr/bin/env python3
"""Create paired base-annotation TSVs for standalone LOFTEE comparison."""

import argparse
import csv
from pathlib import Path

from compare_vep_fastvep_csq import bare_transcript, read_csq


FIELDS = (
    "record_id", "CHROM", "POS", "REF", "ALT", "Feature", "Gene",
    "BIOTYPE", "Consequence", "EXON", "INTRON", "CDS_position",
    "Protein_position", "HGVS_OFFSET",
)


def write_rows(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in FIELDS} for row in rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vep", required=True, type=Path)
    parser.add_argument("--fastvep", required=True, type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    args = parser.parse_args()

    vep = read_csq(args.vep)[1]
    fast = read_csq(args.fastvep)[1]
    vep_rows = []
    fast_rows = []
    hybrid_rows = []
    for record_id, selected_rows in vep.items():
        if len(selected_rows) != 1:
            raise ValueError(f"{record_id}: expected one VEP-picked row")
        selected = selected_rows[0]
        transcript = bare_transcript(selected["Feature"])
        matches = [
            row for row in fast[record_id]
            if bare_transcript(row["Feature"]) == transcript
        ]
        if len(matches) != 1:
            raise ValueError(f"{record_id}: expected one fastVEP row for {transcript}")
        fast_row = matches[0]
        # Use the versioned Feature expected by the transcript database and make
        # base identity identical so only annotation fields can affect LOFTEE.
        fast_row = fast_row.copy()
        fast_row.update({
            field: selected[field]
            for field in ("record_id", "CHROM", "POS", "REF", "ALT", "Feature")
        })
        hybrid_row = fast_row.copy()
        hybrid_row["HGVS_OFFSET"] = selected.get("HGVS_OFFSET", "")
        vep_rows.append(selected)
        fast_rows.append(fast_row)
        hybrid_rows.append(hybrid_row)

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    write_rows(Path(str(args.output_prefix) + ".vep.tsv"), vep_rows)
    write_rows(Path(str(args.output_prefix) + ".fastvep.tsv"), fast_rows)
    write_rows(Path(str(args.output_prefix) + ".fastvep-vep-offset.tsv"), hybrid_rows)
    print(f"rows={len(vep_rows):,}")


if __name__ == "__main__":
    main()
