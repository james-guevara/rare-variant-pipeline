#!/usr/bin/env bash
set -euo pipefail

output=${1:-resources/gene-sets/raw/2026-08-29}
mkdir -p "$output"

fetch() {
  local name=$1
  local url=$2
  if test ! -s "$output/$name"; then
    curl --fail --location --retry 3 --max-time 300 --output "$output/$name" "$url"
  fi
}

fetch sfari-gene-scoring-2026Q2.csv \
  'https://gene-development.sfari.org/wp-content/themes/sfari-gene/utilities/download-csv.php?api-endpoint=genes'

fetch panelapp-intellectual-disability-panel-285-v11.26.json \
  'https://panelapp.genomicsengland.co.uk/api/v1/panels/285/?version=11.26'
fetch panelapp-early-onset-or-syndromic-epilepsy-panel-402-v9.73.json \
  'https://panelapp.genomicsengland.co.uk/api/v1/panels/402/?version=9.73'

fetch schema-phase2-2026-04-23-gene-results.tsv.bgz \
  'https://storage.googleapis.com/exome-results-browsers-public/downloads/2026-04-23/SCHEMA2/SCHEMA2_gene_results.tsv.bgz'
fetch schema-phase2-2026-05-07-gene-results.tsv.bgz \
  'https://storage.googleapis.com/exome-results-browsers-public/downloads/2026-05-07/SCHEMA2/SCHEMA2_gene_results.tsv.bgz'
fetch schema-phase2-2026-05-26-gene-results.tsv.bgz \
  'https://storage.googleapis.com/exome-results-browsers-public/downloads/2026-05-26/SCHEMA2/SCHEMA2_gene_results.tsv.bgz'
fetch schema-phase2-2026-07-15-gene-results.tsv.bgz \
  'https://storage.googleapis.com/exome-results-browsers-public/downloads/2026-07-15/SCHEMA2/SCHEMA2_gene_results.tsv.bgz'
fetch schema-phase2-2026-08-07-gene-results.tsv.bgz \
  'https://storage.googleapis.com/exome-results-browsers-public/downloads/2026-08-07/SCHEMA/SCHEMA_gene_results.tsv.bgz'
fetch schema-phase2-2026-08-21-gene-results.tsv.bgz \
  'https://storage.googleapis.com/exome-results-browsers-public/downloads/2026-08-21/SCHEMA/SCHEMA_gene_results.tsv.bgz'
fetch epi25-2022-12-01-gene-results.tsv.bgz \
  'https://storage.googleapis.com/exome-results-browsers-public/downloads/2022-12-01/Epi25/Epi25_gene_results.tsv.bgz'

fetch kosmicki-height-2026-supplementary-tables.xlsx \
  'https://www.medrxiv.org/content/medrxiv/early/2026/06/24/2026.06.22.26355163/DC2/embed/media-2.xlsx?download=true'
fetch schema-published-2022-article.html \
  'https://pmc.ncbi.nlm.nih.gov/articles/PMC9805802/'
fetch bmi-exome-2021-article.html \
  'https://pmc.ncbi.nlm.nih.gov/articles/PMC10275396/'

(
  cd "$output"
  find . -maxdepth 1 -type f ! -name SHA256SUMS -print0 \
    | sort -z \
    | xargs -0 shasum -a 256 > SHA256SUMS
)
