"""
Join FILTER column from per-chrom VCFs onto merged pipeline output.

For each chromosome, extracts (CHROM, POS, REF, ALT, FILTER) from the
source VCF (split to biallelic with bcftools norm -m-), then left-joins
onto the merged TSV on (CHROM, POS, REF, ALT_specific).

Usage:
    python join_filter.py \
        --vcf-dir /path/to/vcfs \
        --vcf-pattern '{chrom}_jointcall_VQSR_combined.vcf.gz' \
        --merged-dir /path/to/output/indexed \
        --out-dir /path/to/output/with_filter \
        --chroms chr1,chr2,...,chr22,chrX,chrY
"""
import polars as pl
import argparse
import subprocess
import sys
from pathlib import Path


def extract_filter_from_vcf(vcf_path: str) -> pl.DataFrame:
    """Extract (CHROM, POS, REF, ALT, FILTER) from VCF, split multi-allelics."""
    cmd = [
        "bcftools", "norm", "-m-", vcf_path,
        "-Ou",  # uncompressed BCF to pipe
        "|",
        "bcftools", "query", "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\n",
    ]
    result = subprocess.run(
        f"bcftools norm -m- {vcf_path} -Ou | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\n'",
        shell=True, capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"ERROR: bcftools failed: {result.stderr}", file=sys.stderr)
        sys.exit(1)

    lines = result.stdout.strip().split("\n")
    if not lines or lines == [""]:
        return pl.DataFrame(schema={
            "#CHROM": pl.Utf8, "POS": pl.Utf8, "REF": pl.Utf8,
            "ALT_specific": pl.Utf8, "FILTER": pl.Utf8,
        })

    df = pl.read_csv(
        result.stdout.encode(),
        separator="\t",
        has_header=False,
        new_columns=["#CHROM", "POS", "REF", "ALT_specific", "FILTER"],
        infer_schema_length=0,
    )
    # Deduplicate — same site can appear multiple times after norm
    return df.unique(subset=["#CHROM", "POS", "REF", "ALT_specific"])


def main():
    parser = argparse.ArgumentParser(description="Join FILTER from VCFs onto merged output")
    parser.add_argument("--vcf-dir", required=True, help="Directory containing per-chrom VCFs")
    parser.add_argument("--vcf-pattern", required=True,
                        help="VCF filename pattern with {chrom} placeholder")
    parser.add_argument("--merged-dir", required=True, help="Directory with merged TSV.gz files")
    parser.add_argument("--out-dir", required=True, help="Output directory for TSV.gz with FILTER")
    parser.add_argument("--chroms", required=True, help="Comma-separated chromosome list")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    chroms = args.chroms.split(",")

    for chrom in chroms:
        vcf_path = Path(args.vcf_dir) / args.vcf_pattern.replace("{chrom}", chrom)
        merged_path = Path(args.merged_dir) / f"{chrom}.merged.tsv.gz"

        if not vcf_path.exists():
            print(f"WARNING: VCF not found: {vcf_path}, skipping {chrom}", file=sys.stderr)
            continue
        if not merged_path.exists():
            print(f"WARNING: Merged TSV not found: {merged_path}, skipping {chrom}", file=sys.stderr)
            continue

        print(f"{chrom}: extracting FILTER from {vcf_path.name}...", file=sys.stderr)
        filter_df = extract_filter_from_vcf(str(vcf_path))

        print(f"{chrom}: reading {merged_path.name}...", file=sys.stderr)
        merged = pl.read_csv(str(merged_path), separator="\t", infer_schema_length=10000)

        # Ensure POS is string for join (VCF output is string)
        if merged["POS"].dtype != pl.Utf8:
            merged = merged.with_columns(pl.col("POS").cast(pl.Utf8))

        # Left join FILTER
        join_cols = ["#CHROM", "POS", "REF", "ALT_specific"]
        result = merged.join(filter_df, on=join_cols, how="left")

        # Put FILTER after ALT_specific
        cols = result.columns
        if "FILTER" in cols:
            cols.remove("FILTER")
            insert_idx = cols.index("ALT_specific") + 1
            cols.insert(insert_idx, "FILTER")
            result = result.select(cols)

        out_path = out_dir / f"{chrom}.merged.tsv"
        print(f"{chrom}: writing {out_path.name} ({result.height} rows, "
              f"FILTER joined: {result['FILTER'].drop_nulls().len()}/{result.height})...",
              file=sys.stderr)
        result.write_csv(str(out_path), separator="\t")
        print(f"{chrom}: done. bgzip + tabix if needed.", file=sys.stderr)

    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
