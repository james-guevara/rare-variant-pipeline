#!/usr/bin/env python3
"""Attach sex-chromosome QC and analysis-policy fields to carrier rows."""

import argparse
import csv
import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


def in_par(chrom: str, position: int, regions: dict) -> bool:
    return any(
        region["start"] <= position <= region["end"]
        for region in regions["pseudoautosomal_regions"].get(chrom, [])
    )


def policy(chrom: str, position: int, karyotype: str, regions: dict) -> dict:
    par = in_par(chrom, position, regions)
    canonical = regions["par_canonical_representation"]
    duplicate_par = par and chrom != canonical
    if duplicate_par:
        primary = False
        reason = "noncanonical-PAR-representation"
    elif par:
        primary = True
        reason = "PAR-autosomal-policy"
    elif chrom == "chrY":
        primary = karyotype == "XY-like"
        reason = "Y-nonPAR-karyotype-policy"
    else:
        primary = karyotype in {"XX-like", "XY-like"}
        reason = "X-nonPAR-karyotype-policy"
    return {
        "sex_chromosome_region": "PAR" if par else "nonPAR",
        "par_duplicate_excluded": duplicate_par,
        # Noncanonical chrY PAR rows are retained for provenance, but are not
        # burden counts because their chrX homolog is the canonical record.
        "burden_count_available": not duplicate_par,
        "primary_analysis_eligible": primary,
        "frequency_denominator_eligible": primary,
        "sex_chromosome_policy_reason": reason,
        "sensitivity_analysis_group": (
            "ambiguous_karyotype" if karyotype == "ambiguous" else "primary"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--sample-qc", required=True, type=Path)
    parser.add_argument("--regions", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    region_config = json.loads(args.regions.read_text())
    if region_config.get("schema_version") != 1:
        raise ValueError("sex-chromosome region config must use schema_version 1")
    with args.sample_qc.open() as handle:
        qc_rows = list(csv.DictReader(handle, delimiter="\t"))
    qc = {row["sample_id"]: row for row in qc_rows}
    if len(qc) != len(qc_rows):
        raise ValueError("sample QC contains duplicate sample_id values")

    rows = pq.read_table(args.input).to_pylist()
    output = []
    for row in rows:
        sample_id = str(row["sample_id"])
        if sample_id not in qc:
            raise ValueError(f"carrier sample lacks sex-chromosome QC: {sample_id}")
        chrom = str(row.get("#CHROM") or row.get("CHROM"))
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        if chrom not in {"chrX", "chrY"}:
            raise ValueError(f"expected chrX or chrY carrier row, found {chrom}")
        sample_qc = qc[sample_id]
        annotated = dict(row)
        annotated.update({
            "inferred_karyotype": sample_qc["inferred_karyotype"],
            "copy_number_evidence_pattern": sample_qc.get(
                "copy_number_evidence_pattern", ""
            ),
        })
        annotated.update(policy(
            chrom, int(row["POS"]), sample_qc["inferred_karyotype"], region_config
        ))
        output.append(annotated)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    if output:
        table = pa.Table.from_pylist(output)
    else:
        table = pq.read_schema(args.input).empty_table()
    pq.write_table(table, args.output, compression="zstd")
    print(f"carrier_rows={len(output):,}")
    print(f"primary_eligible={sum(r['primary_analysis_eligible'] for r in output):,}")
    print(f"sensitivity_only={sum(not r['primary_analysis_eligible'] and r['burden_count_available'] for r in output):,}")


if __name__ == "__main__":
    main()
