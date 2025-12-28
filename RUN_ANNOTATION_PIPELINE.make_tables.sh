#!/bin/bash
set -euo pipefail

mkdir -p out

#for c in {1..22} X Y; do
#  CHR="chr${c}"
#  IN="out/${CHR}.variants.tsv.gz"
#  OUT="out/${CHR}.variants.reformatted.tsv.gz"
#
#  echo "Reformatting ${CHR}..."
#  python3 reformat_variants.py "${IN}" /dev/stdout | bgzip > "${OUT}"
#done

#for c in {1..22} X Y; do
#  CHR="chr${c}"
#  IN="out/${CHR}.genotypes.tsv.gz"
#  OUT="out/${CHR}.genotypes.resolved.tsv.gz"
#
#  echo "Resolving ${CHR}..."
#  python3 resolve_genotypes.py "${IN}" /dev/stdout | bgzip > "${OUT}"
#done

#for c in {1..22} X Y; do
#  CHR="chr${c}"
#  IN="out/${CHR}.family_genotypes.tsv"
#  OUT="out/${CHR}.family_genotypes.resolved.tsv.gz"
#
#  echo "Resolving ${CHR}..."
#  python3 resolve_family_genotypes.py "${IN}" /dev/stdout | bgzip > "${OUT}" 
#done


#for c in {1..22} X Y; do
#  CHR="chr${c}"
#  echo "Merging ${CHR}..."
#
#  python3 merge_family_variants.py \
#    --family  "out/${CHR}.family_genotypes.resolved.tsv.gz" \
#    --variants "out/${CHR}.variants.reformatted.tsv.gz" \
#    --out /dev/stdout \
#  | bgzip > "out/${CHR}.family_merged.tsv.gz"
#done

