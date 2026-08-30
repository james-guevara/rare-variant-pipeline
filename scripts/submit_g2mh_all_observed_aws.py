#!/usr/bin/env python3
"""Submit the complete G2MH all-observed chromosome workflow to AWS Batch."""

import argparse
import json
import subprocess
from pathlib import Path

from initialize_cohort import GRCH38_LENGTHS, parse_chromosomes


DEFAULT_CHROMOSOMES = ",".join([*(f"chr{i}" for i in range(1, 23)), "chrX", "chrY"])
FSX = "/fsx/loftee-parity"
VCZ = "/fsx/rare-variant-pilot/g2mh-vcz-v3/v1"


def candidate_path(chromosome: str) -> str:
    name = (
        "g2mh.chr22.observed-missense-candidates.parquet"
        if chromosome == "chr22"
        else f"{chromosome}.observed-missense-candidates.parquet"
    )
    return f"{FSX}/resources/targeted-annotation/inputs/{name}"


def postprocess_path(chromosome: str) -> str:
    if chromosome in {"chrX", "chrY"}:
        return f"{FSX}/resources/postprocess/g2mh-{chromosome}/missense-config.json"
    return "/opt/rvp/config/postprocess/aws-g2mh-autosomes.json"


def environment(chromosome: str, run_id: str) -> list[dict[str, str]]:
    values = {
        "ALL_OBSERVED": "1",
        "ANNOTATION_ROOT": f"{FSX}/resources/targeted-annotation/ensembl-115",
        "CHROMOSOME": chromosome,
        "COHORT": "g2mh",
        "COHORT_AF_MAX": "0.01",
        "CONTIG": chromosome.removeprefix("chr"),
        "CONTIG_LENGTH": str(GRCH38_LENGTHS[chromosome]),
        "GENEBAYES": f"{FSX}/resources/targeted-annotation/GeneBayes.Supplementary_Table_1.tsv",
        "LOFTEE_ROOT": f"{FSX}/resources",
        "POPULATION_AF_MAX": "0.01",
        "POSTPROCESS_CONFIG": postprocess_path(chromosome),
        "RUN_ROOT": f"{FSX}/workflows/g2mh/{run_id}/{chromosome}",
        "ZARR_STORE": f"{VCZ}/{chromosome}.sharded-v3.zarr",
    }
    if chromosome.removeprefix("chr").isdigit() and 2 <= int(chromosome[3:]) <= 20:
        values["MISSENSE_DBNSFP"] = (
            f"{FSX}/resources/dbNSFP/5.3.1a/"
            f"parquet_expanded_mane_select/{chromosome}.parquet"
        )
    else:
        values["MISSENSE_CANDIDATES"] = candidate_path(chromosome)
    if chromosome in {"chrX", "chrY"}:
        values.update({
            "SAMPLE_SEX_QC": f"{FSX}/resources/sample-qc/g2mh.sex-chromosome-qc.tsv",
            "SEX_CHROMOSOME_REGIONS": (
                f"{FSX}/resources/sample-qc/grch38-sex-chromosome-regions.json"
            ),
        })
    return [{"name": name, "value": value} for name, value in sorted(values.items())]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--chromosomes", default=DEFAULT_CHROMOSOMES)
    parser.add_argument("--queue", default="rare-variant-vcz-fsx")
    parser.add_argument("--job-definition", default="rare-variant-targeted-portable:11")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    chromosomes = parse_chromosomes(args.chromosomes)
    rows = []
    for chromosome in chromosomes:
        overrides = {
            "command": ["bash", "/opt/rvp/scripts/run_targeted_chromosome.sh"],
            "environment": environment(chromosome, args.run_id),
        }
        request = [
            "aws", "batch", "submit-job",
            "--job-name", f"rvp-g2mh-{chromosome.removeprefix('chr')}-{args.run_id}",
            "--job-queue", args.queue,
            "--job-definition", args.job_definition,
            "--container-overrides", json.dumps(overrides, separators=(",", ":")),
            "--output", "json",
        ]
        job_id = "DRY_RUN"
        if not args.dry_run:
            result = subprocess.run(request, check=True, capture_output=True, text=True)
            job_id = json.loads(result.stdout)["jobId"]
        run_root = next(
            item["value"] for item in overrides["environment"] if item["name"] == "RUN_ROOT"
        )
        rows.append((chromosome, job_id, run_root))
        print(f"submitted\t{chromosome}\t{job_id}\t{run_root}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        "chromosome\tjob_id\trun_root\n"
        + "".join("\t".join(row) + "\n" for row in rows)
    )


if __name__ == "__main__":
    main()
