#!/bin/bash
set -euo pipefail

OUT_TSV="out/all.family_merged.tsv"
OUT_GZ="out/all.family_merged.tsv.gz"

rm -f "${OUT_TSV}" "${OUT_GZ}"

# header from chr1 (avoid pipefail SIGPIPE on zcat)
set +o pipefail
gzip -cd out/chr1.family_merged.tsv.gz | head -n 1 > "${OUT_TSV}"
set -o pipefail

# bodies (no headers)
for c in {1..22} X Y; do
  CHR="chr${c}"
  echo "Concat ${CHR}..."
  gzip -cd "out/${CHR}.family_merged.tsv.gz" | tail -n +2 >> "${OUT_TSV}"
done

bgzip -f "${OUT_TSV}"   # creates out/all.family_merged.tsv.gz
echo "Wrote ${OUT_GZ}"

