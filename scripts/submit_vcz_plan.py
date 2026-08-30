#!/usr/bin/env python3
"""Create or submit AWS Batch VCZ jobs from an inspected cohort preparation plan."""

import argparse
import csv
import json
import re
import subprocess
from pathlib import Path


ADAPTER = r'''set -euo pipefail
case "$CHROMOSOME" in chr[1-9]|chr1[0-9]|chr2[0-2]|chrX|chrY) ;; *) echo "invalid chromosome" >&2; exit 2;; esac
token=${CHROMOSOME#chr}
source_root="$VCZ_WORK_ROOT/source"
conversion_work="$VCZ_WORK_ROOT/conversion"
adapted="$source_root/g2mh_1065_chr${token}.vcf.gz"
subset_checkpoint="${adapted}.complete"
expected="$VCZ_OUTPUT_ROOT/chr${token}.sharded-v3.zarr"
test "$expected" = "$ZARR_STORE" || { echo "output contract mismatch: $expected != $ZARR_STORE" >&2; exit 2; }
mkdir -p "$source_root" "$conversion_work" "$VCZ_OUTPUT_ROOT"
record_count=$(bcftools index -n "$adapted" 2>/dev/null || printf '0')
if test -f "$subset_checkpoint" && test "$record_count" -gt 0; then
  echo "checkpoint_resume\tstage=subset"
else
  rm -f "$subset_checkpoint"
  resolved_contig=""
  for candidate in "$CONTIG" "$CHROMOSOME" "$token"; do
    if bcftools index -s "$INPUT_VCF" | awk -v contig="$candidate" '$1 == contig && $3 > 0 { found=1 } END { exit !found }'; then
      resolved_contig="$candidate"
      break
    fi
  done
  test -n "$resolved_contig" || { echo "no nonempty contig matching $CONTIG, $CHROMOSOME, or $token" >&2; exit 2; }
  echo "contig_resolution\trequested=$CONTIG\tresolved=$resolved_contig"
  tmp="${adapted}.tmp.$$"
  bcftools view --threads "$VCZ_WORKERS" --regions "$resolved_contig" --output-type z --output "$tmp" "$INPUT_VCF"
  mv "$tmp" "$adapted"
  bcftools index --force --threads "$VCZ_WORKERS" "$adapted"
  record_count=$(bcftools index -n "$adapted")
  test "$record_count" -gt 0 || { echo "subset unexpectedly contains zero records" >&2; exit 2; }
  echo "subset_complete\trecords=$record_count"
  touch "$subset_checkpoint"
fi
source_fingerprint=$(sha256sum "$adapted" | awk '{print $1}')
fingerprint_file="$conversion_work/.chr${token}.source.sha256"
previous_fingerprint=$(cat "$fingerprint_file" 2>/dev/null || true)
if test "$previous_fingerprint" != "$source_fingerprint"; then
  echo "checkpoint_invalidate\tstage=conversion\treason=source_changed"
  rm -rf "$conversion_work/chr${token}" "$ZARR_STORE" "${ZARR_STORE}.complete"
  printf '%s\n' "$source_fingerprint" > "${fingerprint_file}.tmp.$$"
  mv "${fingerprint_file}.tmp.$$" "$fingerprint_file"
else
  echo "checkpoint_resume\tstage=conversion\tsource_fingerprint=$source_fingerprint"
fi
export VCZ_SOURCE_ROOT="$source_root"
export VCZ_WORK_ROOT="$conversion_work"
exec /usr/local/bin/run_chromosome.sh "$token"
'''


def read_plan(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("preparation plan has no header")
        rows = list(reader)
    required = {
        "chromosome", "contig", "input_vcf", "zarr_store",
        "vcz_work_root", "preparation_state",
    }
    if not required.issubset(reader.fieldnames):
        raise ValueError(f"preparation plan is missing {sorted(required-set(reader.fieldnames))}")
    return rows


def safe_job_name(cohort, chromosome):
    text = re.sub(r"[^A-Za-z0-9_-]+", "-", f"rv-vcz-{cohort}-{chromosome}")
    return text[:128]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preparation-plan", required=True, type=Path)
    parser.add_argument("--cohort", required=True)
    parser.add_argument("--queue", required=True)
    parser.add_argument("--job-definition", required=True)
    parser.add_argument("--region", default="us-east-1")
    parser.add_argument("--workers", type=int, choices=(16, 32, 64), default=32)
    parser.add_argument("--memory", default="96G")
    parser.add_argument("--validation", choices=("sampled", "full"), default="sampled")
    parser.add_argument("--submit", action="store_true", help="submit jobs; default is dry-run")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    jobs = []
    for row in read_plan(args.preparation_plan):
        if row["preparation_state"] != "READY_FOR_ZARR":
            raise ValueError(
                f"{row['chromosome']} is {row['preparation_state']}, expected READY_FOR_ZARR"
            )
        environment = {
            "CHROMOSOME": row["chromosome"],
            "CONTIG": row["contig"],
            "INPUT_VCF": row["input_vcf"],
            "ZARR_STORE": row["zarr_store"],
            "VCZ_OUTPUT_ROOT": str(Path(row["zarr_store"]).parent),
            "VCZ_WORK_ROOT": row["vcz_work_root"],
            "VCZ_WORKERS": str(args.workers),
            "VCZ_MEMORY": args.memory,
            "VCZ_VALIDATION": args.validation,
        }
        command = [
            "aws", "batch", "submit-job", "--region", args.region,
            "--job-name", safe_job_name(args.cohort, row["chromosome"]),
            "--job-queue", args.queue,
            "--job-definition", args.job_definition,
            "--parameters", json.dumps({"adapter_command": ADAPTER}, separators=(",", ":")),
            "--container-overrides", json.dumps({
                "environment": [
                    {"name": key, "value": value} for key, value in environment.items()
                ],
                "resourceRequirements": [
                    {"type": "VCPU", "value": str(args.workers)},
                    {"type": "MEMORY", "value": str(args.workers * 3584)},
                ],
            }, separators=(",", ":")),
            "--output", "json",
        ]
        record = {
            "chromosome": row["chromosome"],
            "job_name": safe_job_name(args.cohort, row["chromosome"]),
            "input_vcf": row["input_vcf"],
            "zarr_store": row["zarr_store"],
            "status": "PLANNED",
        }
        if args.submit:
            result = subprocess.run(command, capture_output=True, text=True)
            if result.returncode:
                raise RuntimeError(
                    f"Batch submission failed for {row['chromosome']}: {result.stderr.strip()}"
                )
            response = json.loads(result.stdout)
            record.update(status="SUBMITTED", job_id=response["jobId"])
        jobs.append(record)

    document = {
        "schema_version": 1,
        "cohort": args.cohort,
        "mode": "SUBMITTED" if args.submit else "DRY_RUN",
        "queue": args.queue,
        "job_definition": args.job_definition,
        "workers": args.workers,
        "memory": args.memory,
        "validation": args.validation,
        "jobs": jobs,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    print(f"mode={document['mode']}")
    print(f"jobs={len(jobs)}")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
