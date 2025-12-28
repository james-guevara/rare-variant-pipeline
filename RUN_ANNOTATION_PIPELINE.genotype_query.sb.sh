#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=ind-shared
#SBATCH --cpus-per-task=2
#SBATCH --time=04:00:00
#SBATCH --array=1-24
#SBATCH -J bcftools_query
#SBATCH -o logs_bcftools_query/q_%A_%a.out
#SBATCH -e logs_bcftools_query/q_%A_%a.err

set -euo pipefail
SECONDS=0

module load singularitypro

mkdir -p logs_bcftools_query
mkdir -p out

BIND_FOLDER="/expanse/projects/sebat1/"
BCFTOOLS_CONTAINER="containers/bcftools:1.22--h3a4d415_1"

# map array index to chromosome label (and build chr prefix)
CHR_INDEX=${SLURM_ARRAY_TASK_ID}
if   [[ $CHR_INDEX -le 22 ]]; then
  CHR="chr${CHR_INDEX}"
elif [[ $CHR_INDEX -eq 23 ]]; then
  CHR="chrX"
elif [[ $CHR_INDEX -eq 24 ]]; then
  CHR="chrY"
else
  echo "Bad SLURM_ARRAY_TASK_ID: ${CHR_INDEX}" >&2
  exit 1
fi

BED="out/${CHR}.bed"
VCF="vcfs/${CHR}_jointcall_VQSR_combined.vcf.gz"
OUT="out/${CHR}.genotypes.tsv.gz"

if [[ ! -s "${BED}" ]]; then
  echo "ERROR: missing BED: ${BED}" >&2
  exit 1
fi
if [[ ! -s "${VCF}" ]]; then
  echo "ERROR: missing VCF: ${VCF}" >&2
  exit 1
fi

echo "Creating ${CHR} genotypes file..."

# NOTE: remove the space after the backslash on the singularity exec line (it breaks bash)
singularity exec --bind "${BIND_FOLDER}" "${BCFTOOLS_CONTAINER}" \
  bcftools query \
    -R "${BED}" \
    -HH \
    -i 'GT="alt"' \
    -f '[%CHROM\t%POS0\t%END\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\t%GQ\t%DP\t%AD\n]' \
    "${VCF}" \
  | bgzip > "${OUT}"

echo "Wrote ${OUT}"
echo "Job finished in $SECONDS seconds."

