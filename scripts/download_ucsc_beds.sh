#!/usr/bin/env bash
#
# Download UCSC Genome Browser repeat/region BED files for hg38.
#
# Source: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/
#
# The .txt.gz files are MySQL table dumps. Column 0 is always `bin`
# (a UCSC internal spatial index) which we skip. Column positions
# differ per table — see schema files (.sql) on the same server.
#
# Usage:
#   bash scripts/download_ucsc_beds.sh [output_dir]
#
# Default output_dir: resources/repeats/

set -euo pipefail

OUTDIR="${1:-resources/repeats}"
BASE_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "Downloading to: $(pwd)"

# --------------------------------------------------------------------------
# 1. Segmental Duplications
#    Schema: bin, chrom, chromStart, chromEnd, ...
#    Extract: chrom, chromStart, chromEnd (cols 2-4, i.e. cut -f2-4)
# --------------------------------------------------------------------------
echo "Downloading genomicSuperDups..."
curl -sS -O "$BASE_URL/genomicSuperDups.txt.gz"
gunzip -c genomicSuperDups.txt.gz \
  | cut -f2-4 \
  | sort -k1,1 -k2,2n \
  > genomicSuperDups.bed
echo "  $(wc -l < genomicSuperDups.bed) regions"

# --------------------------------------------------------------------------
# 2. Simple Repeats (Tandem Repeats Finder)
#    Schema: bin, chrom, chromStart, chromEnd, name, ...
#    Extract: chrom, chromStart, chromEnd, name (cols 2-5)
# --------------------------------------------------------------------------
echo "Downloading simpleRepeat..."
curl -sS -O "$BASE_URL/simpleRepeat.txt.gz"
gunzip -c simpleRepeat.txt.gz \
  | cut -f2-5 \
  | sort -k1,1 -k2,2n \
  > simpleRepeat.bed
echo "  $(wc -l < simpleRepeat.bed) regions"

# --------------------------------------------------------------------------
# 3. RepeatMasker
#    Schema: bin, swScore, milliDiv, milliDel, milliIns,
#            genoName, genoStart, genoEnd, genoLeft, strand,
#            repName, repClass, repFamily, ...
#    Extract: genoName, genoStart, genoEnd, repName, repClass, repFamily
#             (cols 6-8, 11-13 in 1-indexed)
# --------------------------------------------------------------------------
echo "Downloading rmsk (RepeatMasker)..."
curl -sS -O "$BASE_URL/rmsk.txt.gz"
gunzip -c rmsk.txt.gz \
  | awk -F'\t' 'BEGIN{OFS="\t"} {print $6,$7,$8,$11,$12,$13}' \
  | sort -k1,1 -k2,2n \
  > rmsk.bed
echo "  $(wc -l < rmsk.bed) regions"

# --------------------------------------------------------------------------
# Cleanup raw downloads
# --------------------------------------------------------------------------
echo "Cleaning up .txt.gz files..."
rm -f genomicSuperDups.txt.gz simpleRepeat.txt.gz rmsk.txt.gz

echo ""
echo "Done. Output files:"
ls -lh *.bed
