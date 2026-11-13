#!/usr/bin/env bash
set -euo pipefail

output=${1:-resources/gene-sets/raw/legacy-pipeline-2026-01-15}
remote=${2:-expanse:/expanse/projects/sebat1/s3/data/sebat/nf_rare_spark_wgs/resources/gene_sets}
mkdir -p "$output"

# Exact raw exports and derived lists used by the historical pipeline.
scp -q "$remote"/'*' "$output"/

# Primary publication supplements underlying the compiled Fu/Satterstrom workbook.
curl --fail --location --retry 3 --max-time 300 \
  --output "$output/satterstrom-2020-supplementary-table-2.xlsx" \
  'https://ars.els-cdn.com/content/image/1-s2.0-S0092867419313984-mmc2.xlsx'
curl --fail --location --retry 3 --max-time 300 \
  --output "$output/fu-2022-supplementary-tables.xlsx" \
  'https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-022-01104-0/MediaObjects/41588_2022_1104_MOESM3_ESM.xlsx'

(
  cd "$output"
  find . -maxdepth 1 -type f ! -name SHA256SUMS -print0 \
    | sort -z \
    | xargs -0 shasum -a 256 > SHA256SUMS
)
