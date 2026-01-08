#!/bin/bash
# Strip CSQ annotation from INFO column

INPUT_DIR="${1:-out_merge}"
OUTPUT_DIR="${2:-out_stripped}"

mkdir -p "$OUTPUT_DIR"

for f in "$INPUT_DIR"/*.merged.tsv; do
    chrom=$(basename "$f" .merged.tsv)
    out="$OUTPUT_DIR/${chrom}.merged.tsv"
    echo "Processing $chrom..."
    sed 's/;CSQ=[^\t]*//; s/;*\t/\t/g' "$f" > "$out"
done

echo "Done!"
