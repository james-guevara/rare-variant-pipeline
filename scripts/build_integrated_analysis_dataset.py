#!/usr/bin/env python3
"""Join a complete rare-burden table to a PGS analysis dataset by IID."""

import argparse
import csv
import json
from pathlib import Path


TIERS = ("lof_t1", "lof_t2", "miss_t1", "miss_t2", "miss_t3", "miss_t4")


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"missing TSV header: {path}")
        return reader.fieldnames, list(reader)


def index_unique(path: Path, key: str):
    fields, rows = read_tsv(path)
    if key not in fields:
        raise ValueError(f"{path} is missing {key}")
    indexed = {}
    for row in rows:
        value = row[key]
        if not value:
            raise ValueError(f"{path} contains an empty {key}")
        if value in indexed:
            raise ValueError(f"duplicate {key}={value} in {path}")
        indexed[value] = row
    return fields, rows, indexed


def compatible(pgs, rare, field, iid):
    left, right = pgs.get(field, ""), rare.get(field, "")
    if left and right and left != right:
        raise ValueError(f"{field} differs for IID={iid}: PGS={left!r} rare={right!r}")
    return right or left


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pgs-dataset", required=True, type=Path)
    parser.add_argument("--pgs-dictionary", required=True, type=Path)
    parser.add_argument("--rare-burdens", required=True, type=Path)
    parser.add_argument("--variable-template", required=True, type=Path)
    parser.add_argument("--cohort-id", required=True)
    parser.add_argument(
        "--missing-rare-policy", choices=("error", "exclude"), default="error",
        help="error on PGS participants absent from rare burdens, or exclude and report them",
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--dictionary", required=True, type=Path)
    parser.add_argument("--qc", required=True, type=Path)
    parser.add_argument("--exclusions", required=True, type=Path)
    args = parser.parse_args()

    pgs_fields, pgs_rows, pgs = index_unique(args.pgs_dataset, "IID")
    rare_fields, rare_rows, rare = index_unique(args.rare_burdens, "IID")
    missing_tiers = sorted(set(TIERS) - set(rare_fields))
    if missing_tiers:
        raise ValueError(f"rare burden table is missing tiers: {missing_tiers}")
    missing_rare = sorted(set(pgs) - set(rare))
    if missing_rare and args.missing_rare_policy == "error":
        raise ValueError(
            "PGS participants are absent from the completed rare-burden table: "
            f"count={len(missing_rare)} examples={missing_rare[:5]}"
        )

    pgs_data_fields = [field for field in pgs_fields if field not in {"FID", "IID", "SEX"}]
    output_fields = ["FID", "IID", "SEX", *pgs_data_fields, *TIERS]
    output_rows = []
    for pgs_row in pgs_rows:
        iid = pgs_row["IID"]
        if iid not in rare:
            continue
        rare_row = rare[iid]
        result = {
            "FID": compatible(pgs_row, rare_row, "FID", iid),
            "IID": iid,
            "SEX": compatible(pgs_row, rare_row, "SEX", iid),
        }
        result.update({field: pgs_row[field] for field in pgs_data_fields})
        for tier in TIERS:
            value = rare_row[tier]
            if value == "":
                raise ValueError(
                    f"incomplete rare burden {tier} for PGS participant IID={iid}"
                )
            try:
                int(value)
            except ValueError as error:
                raise ValueError(f"non-integer rare burden {tier}={value!r} for IID={iid}") from error
            result[tier] = value
        output_rows.append(result)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, output_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    template_fields, template_rows = read_tsv(args.variable_template)
    if "variable" not in template_fields:
        raise ValueError("variable template is missing variable")
    template = {row["variable"]: row for row in template_rows}
    pgs_dictionary_fields, pgs_dictionary_rows = read_tsv(args.pgs_dictionary)
    if "variable" not in pgs_dictionary_fields:
        raise ValueError("PGS dictionary is missing variable")
    pgs_dictionary = {row["variable"]: row for row in pgs_dictionary_rows}
    dictionary_rows = []
    for variable in output_fields:
        if variable in template:
            dictionary_rows.append(template[variable])
            continue
        source = pgs_dictionary.get(variable)
        if source is None:
            raise ValueError(f"no dictionary definition for output variable {variable}")
        dictionary_rows.append({
            "variable": variable,
            "role": "predictor" if variable.startswith("PGS_") else "covariate",
            "data_type": source.get("data_type", "string"),
            "nullable": source.get("nullable", "true"),
            "default_analysis_use": "polygenic predictor" if variable.startswith("PGS_") else "covariate",
            "source": source.get("source", "PGS analysis dataset"),
            "description": source.get("description", ""),
        })
    dictionary_fields = [
        "variable", "role", "data_type", "nullable",
        "default_analysis_use", "source", "description",
    ]
    with args.dictionary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, dictionary_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(dictionary_rows)

    with args.exclusions.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, ["IID", "reason"], delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(
            {"IID": iid, "reason": "missing_from_rare_burdens"}
            for iid in missing_rare
        )

    args.qc.write_text(json.dumps({
        "schema_version": 1,
        "cohort_id": args.cohort_id,
        "analysis_universe": "pgs_analysis_dataset",
        "missing_rare_policy": args.missing_rare_policy,
        "pgs_participants": len(pgs_rows),
        "rare_burden_participants": len(rare_rows),
        "integrated_participants": len(output_rows),
        "rare_only_participants_excluded": len(set(rare) - set(pgs)),
        "pgs_participants_missing_rare_burdens": len(missing_rare),
        "burden_variables": list(TIERS),
    }, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
