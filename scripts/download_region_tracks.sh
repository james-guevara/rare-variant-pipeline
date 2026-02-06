#!/usr/bin/env bash
#
# Download region tracks listed in resources/region_tracks.tsv
#
# Reads the manifest and for each track:
#   1. Downloads the file
#   2. Extracts columns per the extraction rule
#   3. Sorts by chrom + start
#   4. Cleans up the raw download
#
# Usage:
#   bash scripts/download_region_tracks.sh [manifest] [output_dir]
#
# Defaults:
#   manifest:   resources/region_tracks.tsv
#   output_dir: resources/regions/
#
# To download only specific tracks:
#   bash scripts/download_region_tracks.sh resources/region_tracks.tsv resources/regions/ genomicSuperDups simpleRepeat

set -euo pipefail

MANIFEST="${1:-resources/region_tracks.tsv}"
OUTDIR="${2:-resources/regions}"
shift 2 2>/dev/null || true
SELECTED=("$@")  # remaining args are track names to download (empty = all)

if [[ ! -f "$MANIFEST" ]]; then
    echo "ERROR: Manifest not found: $MANIFEST" >&2
    exit 1
fi

mkdir -p "$OUTDIR"

echo "Manifest:   $MANIFEST"
echo "Output dir: $OUTDIR"
if [[ ${#SELECTED[@]} -gt 0 ]]; then
    echo "Tracks:     ${SELECTED[*]}"
else
    echo "Tracks:     all"
fi
echo ""

count=0
skipped=0

while IFS=$'\t' read -r name url category extract description; do
    # Skip comments and blank lines
    [[ -z "$name" || "$name" == \#* ]] && continue

    # If specific tracks requested, skip others
    if [[ ${#SELECTED[@]} -gt 0 ]]; then
        match=0
        for sel in "${SELECTED[@]}"; do
            [[ "$name" == "$sel" ]] && match=1 && break
        done
        [[ $match -eq 0 ]] && continue
    fi

    outfile="$OUTDIR/${name}.bed"

    # Skip if already downloaded
    if [[ -f "$outfile" ]]; then
        lines=$(wc -l < "$outfile")
        echo "SKIP  $name ($lines lines already exist)"
        skipped=$((skipped + 1))
        continue
    fi

    echo "GET   $name"
    echo "      $url"

    # Download
    tmpfile=$(mktemp)
    curl -sS -L -o "$tmpfile" "$url"

    # Decompress if gzipped
    if file "$tmpfile" | grep -q gzip; then
        mv "$tmpfile" "${tmpfile}.gz"
        gunzip "${tmpfile}.gz"
    fi

    # Extract columns and sort
    if [[ "$extract" == "none" ]]; then
        # Already BED format
        sort -k1,1 -k2,2n "$tmpfile" > "$outfile"
    elif [[ "$extract" == cut:* ]]; then
        fields="${extract#cut:}"
        cut -f"$fields" "$tmpfile" | sort -k1,1 -k2,2n > "$outfile"
    elif [[ "$extract" == awk:* ]]; then
        expr="${extract#awk:}"
        awk -F'\t' "BEGIN{OFS=\"\t\"} $expr" "$tmpfile" | sort -k1,1 -k2,2n > "$outfile"
    else
        echo "      WARNING: Unknown extract rule '$extract', copying as-is" >&2
        cp "$tmpfile" "$outfile"
    fi

    # Cleanup
    rm -f "$tmpfile"

    lines=$(wc -l < "$outfile")
    echo "      -> $outfile ($lines lines)"
    count=$((count + 1))

done < "$MANIFEST"

echo ""
echo "Done. Downloaded $count, skipped $skipped."
ls -lhS "$OUTDIR"/*.bed 2>/dev/null || echo "(no files)"
