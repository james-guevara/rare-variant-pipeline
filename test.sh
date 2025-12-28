VCF="vcfs/chr22_jointcall_VQSR_combined.vcf.gz"
PED="wgs.psam"
OUT="tests/test.family_genotypes.txt"
BED="tests/chr22.variants.merged.bed"

python family_query.py \
  --vcf "${VCF}" \
  --ped "${PED}" \
  --out "${OUT}" \
  --region "${BED}"

