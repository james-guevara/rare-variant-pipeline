#!/usr/bin/env python3
"""Inspect cohort VCF headers and sample concordance before Zarr conversion."""

import argparse
import csv
import json
import re
import subprocess
from pathlib import Path


CONTIG_RE = re.compile(r"^##contig=<ID=([^,>]+)(?:,length=([0-9]+))?")
FORMAT_RE = re.compile(r"^##FORMAT=<ID=([^,>]+)")
REQUIRED_FORMATS = {"GT", "GQ", "DP"}


def read_tsv(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"missing TSV header: {path}")
        return reader.fieldnames, list(reader)


def header(location):
    result = subprocess.run(
        ["bcftools", "view", "--header-only", location],
        capture_output=True, text=True,
    )
    if result.returncode:
        raise ValueError(f"cannot read VCF header {location}: {result.stderr.strip()}")
    contigs, formats, samples = {}, set(), None
    for line in result.stdout.splitlines():
        match = CONTIG_RE.match(line)
        if match:
            contigs[match.group(1)] = int(match.group(2)) if match.group(2) else None
        match = FORMAT_RE.match(line)
        if match:
            formats.add(match.group(1))
        if line.startswith("#CHROM"):
            fields = line.split("\t")
            samples = fields[9:]
    if samples is None:
        raise ValueError(f"VCF header has no #CHROM line: {location}")
    return contigs, formats, samples


def resolve_contig(chromosome, contigs):
    bare = chromosome.removeprefix("chr")
    candidates = [chromosome, bare]
    present = [candidate for candidate in candidates if candidate in contigs]
    if len(present) != 1:
        raise ValueError(
            f"cannot uniquely resolve {chromosome}; candidates present={present}"
        )
    return present[0]


def write_tsv(path, fields, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preparation-plan", required=True, type=Path)
    parser.add_argument("--sample-manifest", required=True, type=Path)
    parser.add_argument("--output-plan", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()

    plan_fields, plan = read_tsv(args.preparation_plan)
    required_plan = {"chromosome", "contig_length", "input_vcf", "preparation_state"}
    if not required_plan.issubset(plan_fields):
        raise ValueError(f"preparation plan is missing {sorted(required_plan-set(plan_fields))}")
    _, sample_rows = read_tsv(args.sample_manifest)
    expected_samples = [row["IID"] for row in sample_rows]
    if not expected_samples or len(expected_samples) != len(set(expected_samples)):
        raise ValueError("sample manifest must contain unique IID values")

    inspections, cache = [], {}
    for row in plan:
        location = row["input_vcf"]
        if location not in cache:
            cache[location] = header(location)
        contigs, formats, samples = cache[location]
        missing_formats = sorted(REQUIRED_FORMATS - formats)
        if missing_formats:
            raise ValueError(f"VCF {location} is missing FORMAT fields {missing_formats}")
        if "AD" in formats:
            allele_depth_format = "AD"
        elif {"LAD", "LAA"}.issubset(formats):
            allele_depth_format = "LAD+LAA"
        else:
            raise ValueError(
                f"VCF {location} requires FORMAT AD or the localized pair LAD+LAA"
            )
        missing_samples = sorted(set(expected_samples) - set(samples))
        extra_samples = sorted(set(samples) - set(expected_samples))
        if missing_samples or extra_samples:
            raise ValueError(
                f"VCF/PSAM sample mismatch for {location}: "
                f"missing={missing_samples[:5]} extra={extra_samples[:5]}"
            )
        if len(samples) != len(expected_samples):
            raise ValueError(f"VCF {location} contains duplicate sample IDs")
        contig = resolve_contig(row["chromosome"], contigs)
        expected_length = int(row["contig_length"])
        observed_length = contigs[contig]
        if observed_length is not None and observed_length != expected_length:
            raise ValueError(
                f"contig length differs for {row['chromosome']}: "
                f"expected={expected_length} observed={observed_length}"
            )
        row["contig"] = contig
        row["preparation_state"] = "READY_FOR_ZARR"
        inspections.append({
            "chromosome": row["chromosome"],
            "input_vcf": location,
            "resolved_contig": contig,
            "expected_length": expected_length,
            "header_length": observed_length,
            "samples": len(samples),
            "required_formats": sorted(REQUIRED_FORMATS),
            "allele_depth_format": allele_depth_format,
            "status": "PASS",
        })

    args.output_plan.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output_plan, plan_fields, plan)
    args.report.write_text(json.dumps({
        "schema_version": 1,
        "status": "PASS",
        "unique_vcfs_inspected": len(cache),
        "samples": len(expected_samples),
        "chromosomes": inspections,
        "next_stage": "convert_vcf_to_zarr",
    }, indent=2, sort_keys=True) + "\n")
    print("vcf_inspection=PASS")
    print(f"unique_vcfs={len(cache)}")
    print(f"chromosomes={len(inspections)}")
    print(f"samples={len(expected_samples)}")


if __name__ == "__main__":
    main()
