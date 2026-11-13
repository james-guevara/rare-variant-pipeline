#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Generic launcher: runs ONE pipeline stage for a cohort and chromosome list.
# This is the Nextflow *driver*; it stays cheap (2 cpus / 8 GB) because the real
# work happens in the tasks it submits as separate SLURM jobs.
#
#   sbatch -J norm   launch/run_stage.sh g2mh normalize  chr21,chr22
#   sbatch -J annot  launch/run_stage.sh g2mh annotate   chr21,chr22
#   sbatch -J carr   launch/run_stage.sh g2mh carriers   chr21,chr22
#   sbatch -J export launch/run_stage.sh g2mh export     chr21,chr22
#   sbatch -J post   launch/run_stage.sh g2mh postprocess chr21,chr22
#
# Stages must run in that order; each reads the previous stage's published output.
# See docs/running-g2mh.md.
set -euo pipefail

COHORT="${1:?usage: run_stage.sh <cohort-profile> <stage> <chroms>}"
STAGE="${2:?usage: run_stage.sh <cohort-profile> <stage> <chroms>}"
CHROMS="${3:?usage: run_stage.sh <cohort-profile> <stage> <chroms>}"

# Where per-stage outputs go. Override BASE to relocate everything.
BASE="${BASE:-/expanse/projects/sebat1/s3/data/sebat}"
NORMED="$BASE/normed_vcfs/$COHORT"
SITES="$BASE/sites_vcfs/$COHORT"
VEP="$BASE/vep_vcfs/$COHORT"
SPLITVEP="$BASE/split_vep/$COHORT"
RV="$BASE/rv_output/$COHORT"          # reformat/ scatter/ carriers/ merged/ indexed/ parquet/
POST="$BASE/rv_postprocess/$COHORT"

mkdir -p logs
eval "$(micromamba shell hook --shell bash)"
micromamba activate "${NF_ENV:-rvp_env}"
export NXF_TEMP="${NXF_TEMP:-/expanse/projects/sebat1/$USER/nextflow_temp}"
mkdir -p "$NXF_TEMP"

run () { echo "### $(date) :: nextflow run main.nf $*"; nextflow run main.nf "$@"; }

case "$STAGE" in
  normalize)
    run -entry RUN_NORMALIZE -profile "$COHORT" --chroms "$CHROMS" \
        --outdir "$NORMED" --trace_prefix "norm_" ;;
  annotate)
    # Four steps, each with its own outdir and trace so a failure is attributable.
    run -entry RUN_SITES_ONLY -profile "$COHORT" --chroms "$CHROMS" \
        --normed_vcf_dir "$NORMED/norm" --outdir "$SITES" --trace_prefix "sites_"
    run -entry RUN_VEP_ONLY -profile "$COHORT" --chroms "$CHROMS" \
        --sites_vcf_dir "$SITES/sites" --outdir "$VEP" --trace_prefix "vep_"
    run -entry RUN_SPLIT_VEP_ONLY -profile "$COHORT" --chroms "$CHROMS" \
        --vep_vcf_dir "$VEP/vep" --outdir "$SPLITVEP" --trace_prefix "split_vep_"
    run -entry RUN_PREPARE_VARIANTS_ONLY -profile "$COHORT" --chroms "$CHROMS" \
        --split_vep_dir "$SPLITVEP/variants" --mode "${MODE:-coding_and_utr}" \
        --max_cohort_af "${MAX_COHORT_AF:-1.0}" \
        --outdir "$RV" --trace_prefix "prepare_variants_" ;;
  carriers)
    run -entry RUN_CARRIER_EXTRACTION -profile "$COHORT" --chroms "$CHROMS" \
        --normed_vcf_dir "$NORMED/norm" --max_cohort_af "${MAX_COHORT_AF:-1.0}" \
        --outdir "$RV" --trace_prefix "carrier_" ;;
  export)
    run -entry RUN_EXPORT -profile "$COHORT" --chroms "$CHROMS" \
        --outdir "$RV" --trace_prefix "export_" ;;
  postprocess)
    run -entry RUN_POSTPROCESS -profile "$COHORT" --chroms "$CHROMS" \
        --parquet_dir "$RV/parquet" --outdir "$POST" --trace_prefix "postprocess_" ;;
  *)
    echo "unknown stage: $STAGE (normalize|annotate|carriers|export|postprocess)" >&2
    exit 2 ;;
esac

echo "### $(date) :: $STAGE done for $CHROMS"
