#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --array=1-24
#SBATCH --output=logs/convert_parquet_%A_%a.out
#SBATCH --error=logs/convert_parquet_%A_%a.err

# Usage:
#   INPUT_DIR=/path/to/indexed sbatch scripts/convert_to_parquet_array.sh
#
# Converts all per-chromosome .tsv.gz files to .parquet in the same directory.
# Skips chromosomes that already have a .parquet file.

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate python3.12_env_default

CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
CHR="${CHROMS[$SLURM_ARRAY_TASK_ID - 1]}"

INPUT_DIR="${INPUT_DIR:?ERROR: Set INPUT_DIR environment variable}"
INPUT="${INPUT_DIR}/${CHR}.merged.tsv.gz"
OUTPUT="${INPUT_DIR}/${CHR}.merged.parquet"

if [[ ! -f "$INPUT" ]]; then
    echo "Skipping ${CHR}: input not found (${INPUT})"
    exit 0
fi

if [[ -f "$OUTPUT" ]]; then
    echo "Skipping ${CHR}: parquet already exists"
    exit 0
fi

echo "Converting ${CHR}..."
python scripts/tsv_to_parquet.py "$INPUT" -o "$OUTPUT"
echo "Done: ${CHR}"
