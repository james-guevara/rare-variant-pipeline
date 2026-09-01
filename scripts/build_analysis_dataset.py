#!/usr/bin/env python3
"""Assemble optional PGS, rare-variant, and CNV tables on a declared cohort universe."""

import argparse
import csv
import json
from pathlib import Path


IDENTIFIERS = ("FID", "IID", "SEX")
RARE_TIERS = ("lof_t1", "lof_t2", "miss_t1", "miss_t2", "miss_t3", "miss_t4")
DICTIONARY_FIELDS = (
    "variable", "role", "data_type", "nullable",
    "default_analysis_use", "source", "description",
)


def normalize_sex(value):
    normalized = value.strip().upper()
    if normalized in {"1", "M", "MALE"}:
        return "M"
    if normalized in {"2", "F", "FEMALE"}:
        return "F"
    if normalized in {"", "0", "NA", "N/A", ".", "UNKNOWN"}:
        return ""
    return value


def read_tsv(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"missing TSV header: {path}")
        fields = [field.removeprefix("#") for field in reader.fieldnames]
        rows = []
        for raw in reader:
            rows.append({key.removeprefix("#"): value for key, value in raw.items()})
        return fields, rows


def index_unique(path):
    fields, rows = read_tsv(path)
    if "IID" not in fields:
        raise ValueError(f"{path} is missing IID")
    indexed = {}
    for row in rows:
        iid = row["IID"]
        if not iid:
            raise ValueError(f"{path} contains an empty IID")
        if iid in indexed:
            raise ValueError(f"duplicate IID={iid} in {path}")
        indexed[iid] = row
    return fields, rows, indexed


def read_dictionary(path):
    fields, rows = read_tsv(path)
    if "variable" not in fields:
        raise ValueError(f"{path} dictionary is missing variable")
    result = {}
    for row in rows:
        variable = row["variable"]
        if not variable or variable in result:
            raise ValueError(f"invalid or duplicate dictionary variable={variable!r} in {path}")
        result[variable] = row
    return result


def identity_value(manifest_row, component_row, field, component):
    left = manifest_row.get(field, "")
    right = component_row.get(field, "")
    if field == "SEX":
        left = normalize_sex(left)
        right = normalize_sex(right)
    if left and right and left != right:
        raise ValueError(
            f"{field} differs for IID={manifest_row['IID']}: manifest={left!r} "
            f"{component}={right!r}"
        )
    return left or right


def normalize_dictionary(row):
    return {field: row.get(field, "") for field in DICTIONARY_FIELDS}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--participant-manifest", required=True, type=Path)
    parser.add_argument("--pgs-dataset", type=Path)
    parser.add_argument("--pgs-dictionary", type=Path)
    parser.add_argument("--rare-burdens", type=Path)
    parser.add_argument("--cnv-dataset", type=Path)
    parser.add_argument("--cnv-dictionary", type=Path)
    parser.add_argument("--variable-template", required=True, type=Path)
    parser.add_argument("--cohort-id", required=True)
    for component in ("pgs", "rare", "cnv"):
        parser.add_argument(
            f"--missing-{component}-policy",
            choices=("error", "allow", "exclude"), default="allow",
        )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--dictionary", required=True, type=Path)
    parser.add_argument("--qc", required=True, type=Path)
    parser.add_argument("--exclusions", required=True, type=Path)
    args = parser.parse_args()

    if bool(args.pgs_dataset) != bool(args.pgs_dictionary):
        parser.error("--pgs-dataset and --pgs-dictionary must be supplied together")
    if bool(args.cnv_dataset) != bool(args.cnv_dictionary):
        parser.error("--cnv-dataset and --cnv-dictionary must be supplied together")

    manifest_fields, manifest_rows, manifest = index_unique(args.participant_manifest)
    if not set(IDENTIFIERS).issubset(manifest_fields):
        raise ValueError("participant manifest must contain FID, IID, and SEX")
    for row in manifest_rows:
        row["SEX"] = normalize_sex(row.get("SEX", ""))

    template = read_dictionary(args.variable_template)
    components = []
    if args.pgs_dataset:
        fields, rows, indexed = index_unique(args.pgs_dataset)
        components.append(("pgs", fields, rows, indexed, read_dictionary(args.pgs_dictionary)))
    if args.rare_burdens:
        fields, rows, indexed = index_unique(args.rare_burdens)
        missing = sorted(set(RARE_TIERS) - set(fields))
        if missing:
            raise ValueError(f"rare burden table is missing tiers: {missing}")
        components.append(("rare", fields, rows, indexed, template))
    if args.cnv_dataset:
        fields, rows, indexed = index_unique(args.cnv_dataset)
        components.append(("cnv", fields, rows, indexed, read_dictionary(args.cnv_dictionary)))
    if not components:
        raise ValueError("at least one PGS, rare, or CNV dataset is required")

    used_fields = set(IDENTIFIERS)
    component_fields = {}
    dictionary_sources = {}
    component_qc = {}
    excluded = {}
    for name, fields, rows, indexed, definitions in components:
        data_fields = [field for field in fields if field not in IDENTIFIERS]
        collisions = sorted(used_fields.intersection(data_fields))
        if collisions:
            raise ValueError(f"duplicate variables across components: {collisions}")
        undefined = sorted(set(data_fields) - set(definitions))
        if undefined:
            raise ValueError(f"{name} dictionary lacks variables: {undefined}")
        used_fields.update(data_fields)
        component_fields[name] = data_fields
        dictionary_sources.update({field: definitions[field] for field in data_fields})
        missing_ids = [row["IID"] for row in manifest_rows if row["IID"] not in indexed]
        policy = getattr(args, f"missing_{name}_policy")
        if missing_ids and policy == "error":
            raise ValueError(
                f"participant manifest IIDs absent from {name}: "
                f"count={len(missing_ids)} examples={missing_ids[:5]}"
            )
        if policy == "exclude":
            for iid in missing_ids:
                excluded.setdefault(iid, []).append(f"missing_from_{name}")
        component_qc[name] = {
            "participants": len(rows),
            "participants_in_manifest": len(set(indexed).intersection(manifest)),
            "manifest_participants_missing": len(missing_ids),
            "participants_outside_manifest": len(set(indexed) - set(manifest)),
            "missing_policy": policy,
            "variables": data_fields,
        }

    output_fields = [*IDENTIFIERS]
    for name, *_ in components:
        output_fields.extend(component_fields[name])
    output_rows = []
    for manifest_row in manifest_rows:
        iid = manifest_row["IID"]
        if iid in excluded:
            continue
        result = {field: manifest_row.get(field, "") for field in IDENTIFIERS}
        for name, _, _, indexed, _ in components:
            row = indexed.get(iid)
            if row is None:
                result.update({field: "" for field in component_fields[name]})
                continue
            result["FID"] = identity_value(result, row, "FID", name)
            result["SEX"] = identity_value(result, row, "SEX", name)
            for field in component_fields[name]:
                value = row.get(field, "")
                if name == "rare" and field in RARE_TIERS:
                    if value == "":
                        raise ValueError(f"incomplete rare burden {field} for IID={iid}")
                    try:
                        int(value)
                    except ValueError as error:
                        raise ValueError(
                            f"non-integer rare burden {field}={value!r} for IID={iid}"
                        ) from error
                result[field] = value
        output_rows.append(result)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, output_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(output_rows)

    dictionary_rows = []
    for variable in output_fields:
        source = template.get(variable) or dictionary_sources.get(variable)
        if source is None:
            raise ValueError(f"no dictionary definition for output variable {variable}")
        dictionary_rows.append(normalize_dictionary(source))
    with args.dictionary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, DICTIONARY_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(dictionary_rows)

    with args.exclusions.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, ("IID", "reason"), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(
            {"IID": iid, "reason": reason}
            for iid in (row["IID"] for row in manifest_rows)
            for reason in excluded.get(iid, [])
        )

    args.qc.write_text(json.dumps({
        "schema_version": 1,
        "cohort_id": args.cohort_id,
        "analysis_universe": "participant_manifest",
        "manifest_participants": len(manifest_rows),
        "integrated_participants": len(output_rows),
        "excluded_participants": len(excluded),
        "components": component_qc,
    }, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
