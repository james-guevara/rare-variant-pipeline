#!/usr/bin/env bash
#SBATCH --partition=ind-shared
#SBATCH --account=ddp195
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --job-name=target-extract-bench

set -euo pipefail

module load singularitypro/4.1.2
command -v singularity >/dev/null

benchmark_root=${BENCHMARK_ROOT:?set BENCHMARK_ROOT}
extractor=${EXTRACTOR:?set EXTRACTOR}
zarr_store=${ZARR_STORE:?set ZARR_STORE}
target_bed=${TARGET_BED:?set TARGET_BED}
container=${CONTAINER:?set CONTAINER}
chrom=${CHROM:-chr1}

mkdir -p "$benchmark_root"

for workers in 1 4; do
    output="$benchmark_root/optimized-w${workers}.parquet"
    start=$SECONDS
    singularity exec -B /expanse "$container" \
        python "$extractor" \
        --zarr "$zarr_store" \
        --bed "$target_bed" \
        --chrom "$chrom" \
        --workers "$workers" \
        --output "$output"
    elapsed=$((SECONDS - start))
    printf 'workers=%s elapsed_seconds=%s output=%s\n' \
        "$workers" "$elapsed" "$output"
done
