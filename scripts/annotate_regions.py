#!/usr/bin/env python3
"""
Add region-based BED annotations to variant output files.

Reads a manifest (TSV) of region tracks and annotates variants with
overlap counts for each track. Works on parquet or TSV input.

Usage:
    # Annotate with all tracks in manifest
    python annotate_regions.py input.parquet -o output.parquet \
        --manifest resources/region_tracks.tsv \
        --regions-dir resources/regions/

    # Annotate with specific tracks only
    python annotate_regions.py input.parquet -o output.parquet \
        --regions-dir resources/regions/ \
        --beds genomicSuperDups.bed simpleRepeat.bed rmsk.bed

    # Annotate a TSV file
    python annotate_regions.py input.tsv.gz -o output.parquet \
        --manifest resources/region_tracks.tsv \
        --regions-dir resources/regions/

    # Custom coordinate columns
    python annotate_regions.py input.parquet -o output.parquet \
        --manifest resources/region_tracks.tsv \
        --regions-dir resources/regions/ \
        --chrom-col chrom --start-col start --end-col end
"""
import polars as pl
import polars_bio as pb
import argparse
import sys
import time
from pathlib import Path


def load_input(path: str) -> pl.LazyFrame:
    """Load input file as LazyFrame (parquet or TSV)."""
    p = Path(path)
    if p.suffix == ".parquet":
        return pl.scan_parquet(path)
    elif p.suffix == ".gz" or p.suffixes[-2:] == [".tsv", ".gz"]:
        return pl.scan_csv(path, separator="\t", infer_schema_length=10000)
    elif p.suffix in (".tsv", ".csv", ".txt"):
        return pl.scan_csv(path, separator="\t", infer_schema_length=10000)
    else:
        sys.exit(f"ERROR: Unsupported file format: {p.suffix}")


def add_bed_overlap(
    lf: pl.LazyFrame,
    bed_path: str,
    col_name: str,
    chrom_col: str,
    start_col: str,
    end_col: str,
) -> pl.LazyFrame:
    """Add overlap count with a BED file as a new column."""
    coords = lf.select([chrom_col, start_col, end_col]).unique().collect()

    bed = pl.scan_csv(bed_path, separator="\t", has_header=False).select([
        pl.col("column_1").alias("chrom"),
        pl.col("column_2").cast(pl.Int64).alias("start"),
        pl.col("column_3").cast(pl.Int64).alias("end"),
    ])

    v_iv = coords.select([
        pl.col(chrom_col).alias("chrom"),
        pl.col(start_col).cast(pl.Int64).alias("start"),
        pl.col(end_col).cast(pl.Int64).alias("end"),
    ])

    counts = pb.count_overlaps(v_iv, bed, use_zero_based=True, output_type="polars.DataFrame")
    coord_counts = coords.with_columns(counts.select(pl.col("count").alias(col_name)))

    return lf.join(coord_counts.lazy(), on=[chrom_col, start_col, end_col], how="left")


def add_rmsk_overlap(
    lf: pl.LazyFrame,
    bed_path: str,
    chrom_col: str,
    start_col: str,
    end_col: str,
) -> pl.LazyFrame:
    """Add RepeatMasker overlap count plus repClass and repFamily columns.

    rmsk.bed expected format: chrom, start, end, swScore, strand, repName, repClass, repFamily
    """
    coords = lf.select([chrom_col, start_col, end_col]).unique().collect()

    bed = pl.scan_csv(bed_path, separator="\t", has_header=False).select([
        pl.col("column_1").alias("chrom"),
        pl.col("column_2").cast(pl.Int64).alias("start"),
        pl.col("column_3").cast(pl.Int64).alias("end"),
        pl.col("column_4").cast(pl.Int64).alias("swScore"),
        pl.col("column_7").alias("repClass"),
    ])

    v_iv = coords.select([
        pl.col(chrom_col).alias("chrom"),
        pl.col(start_col).cast(pl.Int64).alias("start"),
        pl.col(end_col).cast(pl.Int64).alias("end"),
    ])

    # Get overlapping pairs
    pairs = pb.overlap(v_iv, bed, use_zero_based=True, output_type="polars.DataFrame")

    if pairs.height == 0:
        return lf.with_columns(
            pl.lit(0).alias("rmsk"),
            pl.lit(None).cast(pl.Utf8).alias("rmsk_repClass"),
        )

    # Per-variant: count overlaps, pick the best hit (highest swScore) for class
    summary = (
        pairs.group_by(["chrom_1", "start_1", "end_1"])
        .agg([
            pl.len().alias("rmsk"),
            pl.col("repClass_2").sort_by("swScore_2", descending=True).first().alias("rmsk_repClass"),
        ])
        .rename({"chrom_1": chrom_col, "start_1": start_col, "end_1": end_col})
    )

    # Cast back to match original types
    summary = summary.with_columns([
        pl.col(start_col).cast(coords.schema[start_col]),
        pl.col(end_col).cast(coords.schema[end_col]),
    ])

    return lf.join(summary.lazy(), on=[chrom_col, start_col, end_col], how="left").with_columns(
        pl.col("rmsk").fill_null(0),
    )


def add_gap_overlap(
    lf: pl.LazyFrame,
    bed_path: str,
    chrom_col: str,
    start_col: str,
    end_col: str,
) -> pl.LazyFrame:
    """Add assembly gap overlap count plus gap type column.

    gap.bed expected format: chrom, start, end, type
    """
    coords = lf.select([chrom_col, start_col, end_col]).unique().collect()

    bed = pl.scan_csv(bed_path, separator="\t", has_header=False).select([
        pl.col("column_1").alias("chrom"),
        pl.col("column_2").cast(pl.Int64).alias("start"),
        pl.col("column_3").cast(pl.Int64).alias("end"),
        pl.col("column_4").alias("gap_type"),
    ])

    v_iv = coords.select([
        pl.col(chrom_col).alias("chrom"),
        pl.col(start_col).cast(pl.Int64).alias("start"),
        pl.col(end_col).cast(pl.Int64).alias("end"),
    ])

    pairs = pb.overlap(v_iv, bed, use_zero_based=True, output_type="polars.DataFrame")

    if pairs.height == 0:
        return lf.with_columns(
            pl.lit(0).alias("gap"),
            pl.lit(None).cast(pl.Utf8).alias("gap_type"),
        )

    summary = (
        pairs.group_by(["chrom_1", "start_1", "end_1"])
        .agg([
            pl.len().alias("gap"),
            pl.col("gap_type_2").first().alias("gap_type"),
        ])
        .rename({"chrom_1": chrom_col, "start_1": start_col, "end_1": end_col})
    )

    summary = summary.with_columns([
        pl.col(start_col).cast(coords.schema[start_col]),
        pl.col(end_col).cast(coords.schema[end_col]),
    ])

    return lf.join(summary.lazy(), on=[chrom_col, start_col, end_col], how="left").with_columns(
        pl.col("gap").fill_null(0),
    )


def parse_manifest(manifest_path: str) -> list[dict]:
    """Parse region_tracks.tsv manifest, return list of track dicts."""
    tracks = []
    with open(manifest_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            tracks.append({
                "name": parts[0],
                "url": parts[1],
                "category": parts[2],
                "extract": parts[3],
                "description": parts[4],
            })
    return tracks


def main():
    parser = argparse.ArgumentParser(
        description="Add region-based BED annotations to variant files."
    )
    parser.add_argument("input", help="Input file (parquet, TSV, or gzipped TSV)")
    parser.add_argument("-o", "--output", required=True, help="Output file (parquet or TSV)")
    parser.add_argument("--manifest", help="Region tracks manifest TSV (resources/region_tracks.tsv)")
    parser.add_argument("--regions-dir", required=True, help="Directory containing BED files")
    parser.add_argument("--beds", nargs="+", help="Specific BED files to use (overrides manifest)")
    parser.add_argument("--chrom-col", default="#CHROM", help="Chromosome column name (default: #CHROM)")
    parser.add_argument("--start-col", default="POS0", help="Start position column name (default: POS0)")
    parser.add_argument("--end-col", default="END", help="End position column name (default: END)")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip annotations for columns that already exist in the input")
    args = parser.parse_args()

    regions_dir = Path(args.regions_dir)
    if not regions_dir.is_dir():
        sys.exit(f"ERROR: Regions directory not found: {regions_dir}")

    # Determine which BED files to use
    if args.beds:
        bed_files = []
        for b in args.beds:
            p = regions_dir / b if not Path(b).is_absolute() else Path(b)
            if not p.exists():
                print(f"WARNING: BED file not found, skipping: {p}", file=sys.stderr)
                continue
            name = p.stem  # filename without extension
            bed_files.append((name, str(p)))
    elif args.manifest:
        tracks = parse_manifest(args.manifest)
        bed_files = []
        for t in tracks:
            p = regions_dir / f"{t['name']}.bed"
            if p.exists():
                bed_files.append((t["name"], str(p)))
            else:
                print(f"SKIP  {t['name']} (not downloaded: {p})", file=sys.stderr)
    else:
        sys.exit("ERROR: Provide either --manifest or --beds")

    if not bed_files:
        sys.exit("ERROR: No BED files found to annotate with")

    # Load input
    print(f"Loading {args.input}...", file=sys.stderr)
    lf = load_input(args.input)

    # Check existing columns
    existing_cols = set(lf.collect_schema().names())

    # Validate coordinate columns
    for col in [args.chrom_col, args.start_col, args.end_col]:
        if col not in existing_cols:
            sys.exit(f"ERROR: Column '{col}' not found. Available: {sorted(existing_cols)[:20]}...")

    # Annotate
    total_start = time.time()
    for name, bed_path in bed_files:
        if args.skip_existing and name in existing_cols:
            print(f"SKIP  {name} (column already exists)", file=sys.stderr)
            continue

        t0 = time.time()
        print(f"ADD   {name} <- {bed_path}", file=sys.stderr, end="", flush=True)

        if name == "rmsk":
            lf = add_rmsk_overlap(lf, bed_path, args.chrom_col, args.start_col, args.end_col)
        elif name == "gap":
            lf = add_gap_overlap(lf, bed_path, args.chrom_col, args.start_col, args.end_col)
        else:
            lf = add_bed_overlap(lf, bed_path, name, args.chrom_col, args.start_col, args.end_col)

        print(f"  ({time.time() - t0:.1f}s)", file=sys.stderr)

    # Write output
    out_path = Path(args.output)
    print(f"Writing {out_path}...", file=sys.stderr)

    if out_path.suffix == ".parquet":
        lf.sink_parquet(str(out_path))
    else:
        lf.sink_csv(str(out_path), separator="\t")

    elapsed = time.time() - total_start
    print(f"Done in {elapsed:.1f}s. Annotated with {len(bed_files)} tracks -> {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
