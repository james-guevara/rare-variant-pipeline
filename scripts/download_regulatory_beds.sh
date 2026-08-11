#!/usr/bin/env bash
#
# Download regulatory BED files for --mode regulatory
#
# Downloads 5 tracks to resources/regulatory/:
#   1. phastConsElements.bed    - from UCSC (copy from regions/)
#   2. encodeCcreCombined.bed   - from UCSC (copy from regions/)
#   3. chromhmm_fetal_brain.bed - Roadmap Epigenomics (E081 + E082, active states)
#   4. abc_enhancers.bed        - Nasser et al. 2021 ABC enhancers (brain cell types)
#   5. psychencode.bed          - PsychENCODE brain enhancers
#
# Prerequisites:
#   - bedtools (for merging ChromHMM intervals)
#   - Run download_region_tracks.sh first (for phastCons and cCRE source files)
#   - ~/bin/liftOver + hg19ToHg38 chain file (downloaded automatically if needed for ABC)
#
# Usage:
#   bash scripts/download_regulatory_beds.sh [output_dir]
#   Default output_dir: resources/regulatory/

set -euo pipefail

OUTDIR="${1:-resources/regulatory}"
REGIONS_DIR="${2:-resources/regions}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

mkdir -p "$OUTDIR"

echo "=== Downloading regulatory BED files ==="
echo "Output dir: $OUTDIR"
echo ""

downloaded=0
skipped=0
failed=0

# Helper: report status
report() {
    local name="$1" lines="$2"
    echo "      -> $OUTDIR/$name ($lines lines)"
    downloaded=$((downloaded + 1))
}

skip() {
    local name="$1"
    local lines
    lines=$(wc -l < "$OUTDIR/$name")
    echo "SKIP  $name ($lines lines already exist)"
    skipped=$((skipped + 1))
}

fail() {
    local name="$1" msg="$2"
    echo "FAIL  $name: $msg" >&2
    failed=$((failed + 1))
}

# ============================================================
# 1. phastConsElements.bed (copy from regions/)
# ============================================================
echo "--- phastConsElements.bed ---"
if [[ -s "$OUTDIR/phastConsElements.bed" ]]; then
    skip "phastConsElements.bed"
elif [[ -f "$REGIONS_DIR/phastConsElements100way.bed" ]]; then
    echo "COPY  phastConsElements100way.bed -> phastConsElements.bed"
    cp "$REGIONS_DIR/phastConsElements100way.bed" "$OUTDIR/phastConsElements.bed"
    report "phastConsElements.bed" "$(wc -l < "$OUTDIR/phastConsElements.bed")"
else
    fail "phastConsElements.bed" "Source not found: $REGIONS_DIR/phastConsElements100way.bed (run download_region_tracks.sh first)"
fi

# ============================================================
# 2. encodeCcreCombined.bed (copy from regions/)
# ============================================================
echo "--- encodeCcreCombined.bed ---"
if [[ -s "$OUTDIR/encodeCcreCombined.bed" ]]; then
    skip "encodeCcreCombined.bed"
elif [[ -f "$REGIONS_DIR/encodeCcreCombined.bed" ]]; then
    lines=$(wc -l < "$REGIONS_DIR/encodeCcreCombined.bed")
    if [[ "$lines" -lt 100 ]]; then
        fail "encodeCcreCombined.bed" "Source file has only $lines lines (expected ~926K). Re-download via bigBedToBed."
    else
        echo "COPY  encodeCcreCombined.bed"
        cp "$REGIONS_DIR/encodeCcreCombined.bed" "$OUTDIR/encodeCcreCombined.bed"
        report "encodeCcreCombined.bed" "$lines"
    fi
else
    fail "encodeCcreCombined.bed" "Source not found: $REGIONS_DIR/encodeCcreCombined.bed (run download_region_tracks.sh first)"
fi

# ============================================================
# 3. chromhmm_fetal_brain.bed (Roadmap Epigenomics E081 + E082)
# ============================================================
echo "--- chromhmm_fetal_brain.bed ---"
if [[ -s "$OUTDIR/chromhmm_fetal_brain.bed" ]]; then
    skip "chromhmm_fetal_brain.bed"
else
    if ! command -v bedtools &>/dev/null; then
        fail "chromhmm_fetal_brain.bed" "bedtools not found (required for merging intervals)"
    else
        echo "GET   ChromHMM 15-state: E081 (Fetal Brain Female) + E082 (Fetal Brain Male)"
        tmpdir=$(mktemp -d)

        # Download both fetal brain epigenomes
        echo "      Downloading E081..."
        curl -sS -L -o "$tmpdir/E081.bed.gz" \
            "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E081_15_coreMarks_hg38lift_mnemonics.bed.gz"

        echo "      Downloading E082..."
        curl -sS -L -o "$tmpdir/E082.bed.gz" \
            "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E082_15_coreMarks_hg38lift_mnemonics.bed.gz"

        # Verify downloads aren't HTML error pages
        for f in "$tmpdir/E081.bed.gz" "$tmpdir/E082.bed.gz"; do
            if ! file "$f" | grep -q gzip; then
                fail "chromhmm_fetal_brain.bed" "Download failed (not gzip): $f"
                rm -rf "$tmpdir"
                continue 2  # This doesn't work in bash for outer loops, handle below
            fi
        done

        # Filter to active chromatin states: 1_TssA, 6_EnhG, 7_Enh
        echo "      Filtering to active states (1_TssA, 6_EnhG, 7_Enh) and merging..."
        gunzip -c "$tmpdir/E081.bed.gz" "$tmpdir/E082.bed.gz" \
            | awk -F'\t' '$4 ~ /^(1_TssA|6_EnhG|7_Enh)$/' \
            | cut -f1-3 \
            | sort -k1,1 -k2,2n \
            | bedtools merge \
            > "$OUTDIR/chromhmm_fetal_brain.bed"

        rm -rf "$tmpdir"
        report "chromhmm_fetal_brain.bed" "$(wc -l < "$OUTDIR/chromhmm_fetal_brain.bed")"
    fi
fi

# ============================================================
# 4. abc_enhancers.bed (Nasser et al. 2021)
# ============================================================
echo "--- abc_enhancers.bed ---"
if [[ -s "$OUTDIR/abc_enhancers.bed" ]]; then
    skip "abc_enhancers.bed"
else
    echo "GET   ABC enhancer-gene predictions (Nasser et al. 2021)"
    echo "      WARNING: AllPredictions file is ~3.5GB compressed — this may take a while"
    tmpdir=$(mktemp -d)

    echo "      Downloading AllPredictions.AvgHiC.ABC0.015..."
    curl -sS -L -o "$tmpdir/AllPredictions.txt.gz" \
        "https://mitra.stanford.edu/engreitz/oak/public/Nasser2021/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz"

    if ! file "$tmpdir/AllPredictions.txt.gz" | grep -q gzip; then
        fail "abc_enhancers.bed" "Download failed (not gzip)"
        rm -rf "$tmpdir"
    else
        # Check coordinate system (hg19 vs hg38) from first data line
        echo "      Checking coordinate system..."
        first_chr=$(gunzip -c "$tmpdir/AllPredictions.txt.gz" | head -2 | tail -1 | cut -f1)
        echo "      First chromosome: $first_chr"

        # Determine if hg19 (chr1) or hg38 (chr1) — both use "chr" prefix
        # The Nasser 2021 data uses hg19. We need to liftOver to hg38.
        # Extract brain-related biosamples first, then liftOver.

        echo "      Filtering to brain biosamples..."
        # Header line has column names. Extract brain-related cell types.
        # Key columns: chr(1), start(2), end(3), TargetGene(7), ABC.Score(21), CellType(10)
        # Get header to find exact column indices
        header=$(gunzip -c "$tmpdir/AllPredictions.txt.gz" | head -1)

        # Use Python for robust column-index finding and filtering
        python3 -c "
import gzip, sys

with gzip.open('$tmpdir/AllPredictions.txt.gz', 'rt') as f:
    header = f.readline().strip().split('\t')
    # Find column indices
    cols = {name: i for i, name in enumerate(header)}
    chr_i = cols['chr']
    start_i = cols['start']
    end_i = cols['end']
    gene_i = cols['TargetGene']
    score_i = cols['ABC.Score']
    cell_i = cols['CellType']

    brain_terms = ['brain', 'neuron', 'astrocyte', 'oligodendrocyte', 'fetal_brain',
                   'cortex', 'cerebellum', 'hippocampus', 'neural']

    out = open('$tmpdir/abc_brain_hg19.bed', 'w')
    count = 0
    for line in f:
        fields = line.strip().split('\t')
        cell_type = fields[cell_i].lower()
        if any(term in cell_type for term in brain_terms):
            out.write(f'{fields[chr_i]}\t{fields[start_i]}\t{fields[end_i]}\t{fields[gene_i]}\t{fields[score_i]}\n')
            count += 1
    out.close()
    print(f'      Extracted {count} brain enhancer-gene predictions', file=sys.stderr)
" 2>&1

        brain_lines=$(wc -l < "$tmpdir/abc_brain_hg19.bed")
        echo "      Brain predictions: $brain_lines"

        if [[ "$brain_lines" -eq 0 ]]; then
            fail "abc_enhancers.bed" "No brain biosamples found in ABC data"
            rm -rf "$tmpdir"
        else
            # LiftOver from hg19 to hg38
            LIFTOVER="$HOME/bin/liftOver"
            CHAIN="$HOME/bin/hg19ToHg38.over.chain.gz"

            if [[ ! -x "$LIFTOVER" ]]; then
                echo "      Downloading liftOver..."
                curl -sS -L -o "$LIFTOVER" "https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/liftOver"
                chmod +x "$LIFTOVER"
            fi

            if [[ ! -f "$CHAIN" ]]; then
                echo "      Downloading hg19ToHg38 chain file..."
                curl -sS -L -o "$CHAIN" "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
            fi

            echo "      Running liftOver hg19 -> hg38..."
            "$LIFTOVER" "$tmpdir/abc_brain_hg19.bed" "$CHAIN" \
                "$tmpdir/abc_brain_hg38.bed" "$tmpdir/abc_unmapped.bed"

            unmapped=$(wc -l < "$tmpdir/abc_unmapped.bed")
            lifted=$(wc -l < "$tmpdir/abc_brain_hg38.bed")
            echo "      Lifted: $lifted, unmapped: $((unmapped / 2))"

            sort -k1,1 -k2,2n "$tmpdir/abc_brain_hg38.bed" > "$OUTDIR/abc_enhancers.bed"
            rm -rf "$tmpdir"
            report "abc_enhancers.bed" "$(wc -l < "$OUTDIR/abc_enhancers.bed")"
        fi
    fi
fi

# ============================================================
# 5. psychencode.bed (PsychENCODE Resource Portal)
# ============================================================
echo "--- psychencode.bed ---"
if [[ -s "$OUTDIR/psychencode.bed" ]]; then
    skip "psychencode.bed"
else
    echo "GET   PsychENCODE brain enhancers (DER-04a, hg38 liftover)"
    tmpfile=$(mktemp)

    # Try primary source
    echo "      Trying resource.psychencode.org..."
    http_code=$(curl -sS -L -o "$tmpfile" -w "%{http_code}" \
        "http://resource.psychencode.org/Datasets/Derived/DER-04a_hg38lft_PEC_enhancers.bed" 2>/dev/null || echo "000")

    if [[ "$http_code" == "200" ]] && [[ -s "$tmpfile" ]] && ! head -1 "$tmpfile" | grep -q "<!DOCTYPE"; then
        echo "      Downloaded from primary source"
        sort -k1,1 -k2,2n "$tmpfile" > "$OUTDIR/psychencode.bed"
        rm -f "$tmpfile"
        report "psychencode.bed" "$(wc -l < "$OUTDIR/psychencode.bed")"
    else
        rm -f "$tmpfile"
        # Try Synapse fallback
        echo "      Primary source failed (HTTP $http_code). Trying Synapse..."
        if command -v synapse &>/dev/null; then
            synapse get syn12162741 --downloadLocation "$OUTDIR/" 2>/dev/null
            if [[ -f "$OUTDIR/DER-04a_hg38lft_PEC_enhancers.bed" ]]; then
                sort -k1,1 -k2,2n "$OUTDIR/DER-04a_hg38lft_PEC_enhancers.bed" > "$OUTDIR/psychencode.bed"
                rm -f "$OUTDIR/DER-04a_hg38lft_PEC_enhancers.bed"
                report "psychencode.bed" "$(wc -l < "$OUTDIR/psychencode.bed")"
            else
                fail "psychencode.bed" "Synapse download failed"
            fi
        else
            fail "psychencode.bed" "Primary source down and synapse CLI not installed. Install: pip install synapseclient"
        fi
    fi
fi

# ============================================================
# Summary
# ============================================================
echo ""
echo "=== Summary ==="
echo "Downloaded: $downloaded"
echo "Skipped:    $skipped"
echo "Failed:     $failed"
echo ""

if [[ -d "$OUTDIR" ]]; then
    echo "Files in $OUTDIR:"
    ls -lhS "$OUTDIR"/*.bed 2>/dev/null || echo "(no files)"
fi

if [[ $failed -gt 0 ]]; then
    echo ""
    echo "WARNING: $failed track(s) failed. Check errors above."
    exit 1
fi
