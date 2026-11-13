#!/usr/bin/env bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --job-name=pull-targeted-container

set -euo pipefail

image=${1:?usage: pull_targeted_container_expanse.sh IMAGE SIF}
sif=${2:?usage: pull_targeted_container_expanse.sh IMAGE SIF}

module load singularitypro/4.1.2
command -v singularity >/dev/null
mkdir -p "$(dirname "$sif")"

cache_root=/tmp/rvp-targeted-pull-${SLURM_JOB_ID:?}
export SINGULARITY_CACHEDIR=$cache_root/cache
export SINGULARITY_TMPDIR=$cache_root/tmp
mkdir -p "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR"
trap 'rm -rf "$cache_root"' EXIT

if test -s "$sif"; then
    printf 'checkpoint\tSKIP_PULL\t%s\n' "$sif"
else
    singularity pull "$sif" "$image"
fi
singularity exec "$sif" bash -c '
    set -euo pipefail
    cd /opt/rvp
    fastvep --version
    bcftools --version | head -n1
    python -c "import duckdb, pyarrow, pyBigWig, pysam, zarr"
    python scripts/extract_zarr_target_alleles.py --help >/dev/null
    test -r resources/g2mh-chr1-targeted-regression.json
    test -r resources/g2mh-chrX-targeted-regression.json
'
sha256sum "$sif"
