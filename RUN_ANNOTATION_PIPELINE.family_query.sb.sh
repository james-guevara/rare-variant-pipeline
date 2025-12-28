#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=ind-shared
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --array=1-24
#SBATCH -J family_query
#SBATCH -o logs_family_query/fam_%A_%a.out
#SBATCH -e logs_family_query/fam_%A_%a.err

set -euo pipefail
SECONDS=0

mkdir -p logs_family_query
mkdir -p out

# map array index to chromosome label
CHR_INDEX=${SLURM_ARRAY_TASK_ID}
if   [[ $CHR_INDEX -le 22 ]]; then
  CHR=$CHR_INDEX
elif [[ $CHR_INDEX -eq 23 ]]; then
  CHR="X"
elif [[ $CHR_INDEX -eq 24 ]]; then
  CHR="Y"
fi

VCF="vcfs/chr${CHR}_jointcall_VQSR_combined.vcf.gz"
PED="wgs.psam"
BED="out/chr${CHR}.bed"
OUT="out/chr${CHR}.family_genotypes.tsv"

eval "$(micromamba shell hook --shell=bash)"
micromamba activate python3.12_env_default

python family_query.py \
  --vcf "${VCF}" \
  --ped "${PED}" \
  --out "${OUT}" \
  --region "${BED}"

echo "Job finished in $SECONDS seconds."

