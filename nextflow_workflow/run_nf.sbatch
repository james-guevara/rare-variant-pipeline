#!/bin/bash
#SBATCH --job-name=nf_rare
#SBATCH --account=ddp195
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=08:00:00
#SBATCH --output=logs/nf_orchestrator_%j.out
#SBATCH --error=logs/nf_orchestrator_%j.err

CHROM="${1:-chrY}"
SCRIPT_DIR="/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/rare_variant_pipeline/nextflow_workflow"

cd "${SCRIPT_DIR}"
mkdir -p logs reports

echo "Starting Nextflow for ${CHROM} at $(date)"

micromamba run -n nf_env nextflow run main.nf \
    --chrom "${CHROM}" \
    --outdir "${SCRIPT_DIR}/output" \
    -resume

echo "Completed at $(date)"
