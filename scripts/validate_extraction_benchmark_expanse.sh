#!/usr/bin/env bash
#SBATCH --partition=ind-shared
#SBATCH --account=ddp195
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --job-name=target-extract-validate

set -euo pipefail

module load singularitypro/4.1.2
command -v singularity >/dev/null

benchmark_root=${BENCHMARK_ROOT:?set BENCHMARK_ROOT}
container=${CONTAINER:?set CONTAINER}
comparator=$benchmark_root/compare_parquet_regression.py

for workers in 1 4; do
    singularity exec -B /expanse "$container" \
        python "$comparator" \
        --new "$benchmark_root/optimized-w${workers}.parquet" \
        --reference "$benchmark_root/current.parquet"
done
