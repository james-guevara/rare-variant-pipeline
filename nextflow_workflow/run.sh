#!/bin/bash
# Run Nextflow pipeline
# Usage: ./run.sh [CHROM]

CHROM="${1:-chrY}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

cd "${SCRIPT_DIR}"
mkdir -p reports

echo "Running Nextflow pipeline for ${CHROM}"

micromamba run -n nf_env nextflow run main.nf --chrom "${CHROM}" -resume
