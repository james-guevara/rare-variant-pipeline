#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=pull_python_container_%j.out
#SBATCH --error=pull_python_container_%j.err
#SBATCH --job-name=pull_rvp_py
#
# One-time setup: fetch the Python-stages container published by CI and write it
# next to the other pipeline containers.
#
# Run this as a JOB, not on a login node. `singularity pull` shells out to mksquashfs,
# which dies with "Out of memory (frag_thrd)" under login-node limits. It succeeds
# with 32 GB in a job.
#
# Usage:  sbatch scripts/pull_python_container.sh [IMAGE_REF]
set -euo pipefail

module load singularitypro/4.1.2

IMAGE="${1:-docker://ghcr.io/james-guevara/rare-variant-pipeline-python:main}"
CONTAINERS="${CONTAINERS_DIR:-/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/rare_variant_pipeline/containers}"
SIF="$CONTAINERS/rvp_python.sif"

export SINGULARITY_CACHEDIR="/tmp/sc_${SLURM_JOB_ID:-$$}"
export SINGULARITY_TMPDIR="/tmp/st_${SLURM_JOB_ID:-$$}"
mkdir -p "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR" "$CONTAINERS"

echo "### image: $IMAGE"
echo "### sif:   $SIF"
singularity pull --force "$SIF" "$IMAGE"

echo; echo "### verify"
ls -la "$SIF"
singularity exec "$SIF" python -c "
import duckdb, polars, polars_bio, pyarrow
print('  duckdb    ', duckdb.__version__)
print('  polars    ', polars.__version__)
print('  polars_bio', polars_bio.__version__)
print('  pyarrow   ', pyarrow.__version__)
"
echo "### expected: duckdb 1.5.2  polars 1.39.3  polars_bio 0.28.0  pyarrow 22.0.0"

rm -rf "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR"
echo "### done. Point params.python_container at: $SIF"
