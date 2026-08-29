#!/usr/bin/env python3
"""Initialize a portable rare-variant cohort from VCF input and a PLINK PSAM."""

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


GRCH38_LENGTHS = {
    **{f"chr{i}": length for i, length in enumerate((
        248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
        159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
        114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
        58617616, 64444167, 46709983, 50818468,
    ), start=1)},
    "chrX": 156040895,
    "chrY": 57227415,
}


def parse_chromosomes(value):
    chromosomes = [item.strip() for item in value.split(",") if item.strip()]
    normalized = [item if item.startswith("chr") else f"chr{item}" for item in chromosomes]
    unknown = sorted(set(normalized) - set(GRCH38_LENGTHS))
    if unknown:
        raise ValueError(f"unsupported GRCh38 chromosomes: {unknown}")
    if len(normalized) != len(set(normalized)):
        raise ValueError("chromosome list contains duplicates")
    return normalized


def normalize_sex(value):
    text = value.strip().upper()
    if text in {"1", "M", "MALE"}:
        return "M"
    if text in {"2", "F", "FEMALE"}:
        return "F"
    if text in {"", "0", "NA", "N/A", ".", "UNKNOWN", "-9"}:
        return ""
    return value.strip()


def read_psam(path):
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"empty PSAM: {path}")
    header = lines[0].split()
    header[0] = header[0].removeprefix("#")
    if "IID" not in header:
        raise ValueError("PSAM must contain IID")
    rows, seen = [], set()
    for number, line in enumerate(lines[1:], start=2):
        values = line.split()
        if len(values) != len(header):
            raise ValueError(f"PSAM row {number} has {len(values)} values; expected {len(header)}")
        row = dict(zip(header, values))
        iid = row["IID"]
        if not iid or iid in seen:
            raise ValueError(f"empty or duplicate PSAM IID at row {number}: {iid!r}")
        seen.add(iid)
        rows.append({
            "FID": "" if row.get("FID", "") in {"0", "."} else row.get("FID", ""),
            "IID": iid,
            "SEX": normalize_sex(row.get("SEX", "")),
        })
    if not rows:
        raise ValueError("PSAM contains no samples")
    return rows


def write_tsv(path, fields, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", required=True)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--joint-vcf", help="one joint VCF containing all requested chromosomes")
    inputs.add_argument("--vcf-template", help="per-chromosome VCF path/URI containing {chrom}")
    parser.add_argument("--psam", required=True, type=Path)
    parser.add_argument("--reference-build", default="GRCh38", choices=("GRCh38",))
    parser.add_argument("--chromosomes", default=",".join(GRCH38_LENGTHS))
    parser.add_argument("--shared-resources-root", required=True)
    parser.add_argument("--cohort-root", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--pass-policy", choices=("auto", "required", "ignore"), default="auto",
        help="whether the eventual genotype pipeline requires VCF FILTER=PASS",
    )
    args = parser.parse_args()
    if args.vcf_template and "{chrom}" not in args.vcf_template:
        raise ValueError("--vcf-template must contain {chrom}")

    chromosomes = parse_chromosomes(args.chromosomes)
    samples = read_psam(args.psam)
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    sample_manifest = output / "sample_manifest.tsv"
    write_tsv(sample_manifest, ["FID", "IID", "SEX"], samples)

    shared = args.shared_resources_root.rstrip("/")
    cohort_root = args.cohort_root.rstrip("/")
    rows = []
    for chrom in chromosomes:
        input_vcf = args.joint_vcf or args.vcf_template.format(chrom=chrom)
        rows.append({
            "chromosome": chrom,
            "contig": "auto",
            "contig_length": GRCH38_LENGTHS[chrom],
            "input_vcf": input_vcf,
            "zarr_store": f"{cohort_root}/zarr/{chrom}.sharded-v3.zarr",
            "target_bed": f"{shared}/targeted-annotation/inputs/lof-plus-missense-candidates.{chrom}.bed",
            "missense_candidates": f"{cohort_root}/candidates/{chrom}.observed-missense-candidates.parquet",
            "annotation_root": f"{shared}/targeted-annotation/ensembl-115",
            "loftee_root": shared,
            "genebayes": f"{shared}/targeted-annotation/GeneBayes.Supplementary_Table_1.tsv",
            "postprocess_config": f"{cohort_root}/postprocess/{chrom}/config.json",
            "run_root": f"{cohort_root}/workflows/{chrom}",
            "manifest_output": f"{output}/manifests/{args.cohort}-{chrom}.json",
            "binding_output": f"{output}/bindings/{args.cohort}-{chrom}.json",
            "preparation_state": "PENDING_DERIVED_RESOURCES",
        })
    plan = output / "chromosome_preparation.tsv"
    write_tsv(plan, list(rows[0]), rows)

    declaration = {
        "schema_version": 1,
        "cohort": args.cohort,
        "reference_build": args.reference_build,
        "input": {
            "kind": "joint_vcf" if args.joint_vcf else "per_chromosome_vcf_template",
            "location": args.joint_vcf or args.vcf_template,
            "psam": str(args.psam),
        },
        "chromosomes": chromosomes,
        "shared_resources_root": args.shared_resources_root,
        "cohort_root": args.cohort_root,
        "pass_policy": args.pass_policy,
        "generated": {
            "sample_manifest": str(sample_manifest),
            "chromosome_preparation": str(plan),
        },
    }
    (output / "cohort.json").write_text(json.dumps(declaration, indent=2, sort_keys=True) + "\n")
    sex_counts = Counter(row["SEX"] or "unknown" for row in samples)
    (output / "initialization_qc.json").write_text(json.dumps({
        "schema_version": 1,
        "cohort": args.cohort,
        "samples": len(samples),
        "family_ids_nonempty": sum(bool(row["FID"]) for row in samples),
        "sex_counts": dict(sorted(sex_counts.items())),
        "chromosomes_requested": len(chromosomes),
        "derived_resources_ready": False,
        "next_stage": "inspect_vcf_and_prepare_zarr",
    }, indent=2, sort_keys=True) + "\n")
    print(f"cohort={args.cohort}")
    print(f"samples={len(samples)}")
    print(f"chromosomes={len(chromosomes)}")
    print("derived_resources_ready=false")
    print(f"output_dir={output}")


if __name__ == "__main__":
    main()
