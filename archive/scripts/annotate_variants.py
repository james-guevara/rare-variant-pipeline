#!/usr/bin/env python3
"""
Add variant-level annotations (REVEL, popEVE, pext) to pipeline output parquet files.

Joins on chrom/pos/ref/alt (REVEL, popEVE) or chrom/pos (pext).
Runs as a post-pipeline step on annotated parquet files.
All joins are performed in DuckDB for memory-efficient streaming.

Usage:
    # Add all three annotations
    python annotate_variants.py input.parquet -o output.parquet \
        --revel /path/to/dbNSFP5.3a/dbNSFP5.3a_variant.chr{chrom}.gz \
        --popeve /path/to/popEVE/grch38_popEVE_ukbb_20250715.vcf.gz \
        --pext /path/to/gnomAD/pext/gnomad.pext.gtex_v10.base_level.tsv.gz

    # Add only REVEL
    python annotate_variants.py input.parquet -o output.parquet \
        --revel /path/to/dbNSFP5.3a/dbNSFP5.3a_variant.chr{chrom}.gz

    # Custom coordinate columns
    python annotate_variants.py input.parquet -o output.parquet \
        --revel /path/to/dbNSFP5.3a_variant.chr{chrom}.gz \
        --chrom-col '#CHROM' --pos-col POS --ref-col REF --alt-col ALT
"""
import duckdb
import argparse
import sys
import time
from pathlib import Path


def resolve_revel_path(template: str, chrom: str) -> str | None:
    """Resolve REVEL file path, trying bare and prefixed chromosome."""
    chrom_bare = chrom.replace("chr", "")
    path = template.replace("{chrom}", chrom_bare)
    if Path(path).exists():
        return path
    path = template.replace("{chrom}", chrom)
    if Path(path).exists():
        return path
    print(f"  WARNING: REVEL file not found: {path}", file=sys.stderr)
    return None


def add_revel(
    con: duckdb.DuckDBPyConnection,
    current: str,
    revel_path: str,
    chrom_col: str,
    pos_col: str,
    ref_col: str,
    alt_col: str,
    tx_col: str = "Feature",
) -> str:
    """Add REVEL scores via DuckDB LEFT JOINs. Returns new view name.

    Adds three columns:
      - REVEL_score_tx:  matched to the specific transcript (Feature ↔ Ensembl_transcriptid)
      - REVEL_score_max: max across all transcripts at this variant
      - REVEL_score_min: min across all transcripts at this variant
    """
    # First create a view of the raw dbNSFP REVEL data
    con.execute(f"""
        CREATE VIEW _revel_raw AS
        SELECT
            'chr' || "#chr" AS _chrom,
            "pos(1-based)" AS _pos,
            ref AS _ref,
            alt AS _alt,
            Ensembl_transcriptid AS _tx,
            TRY_CAST(REVEL_score AS DOUBLE) AS REVEL_score
        FROM read_csv_auto('{revel_path}', delim='\t', header=true,
             all_varchar=true)
        WHERE REVEL_score != '.'
    """)

    # Aggregate min/max per variant
    con.execute("""
        CREATE VIEW _revel_agg AS
        SELECT _chrom, _pos, _ref, _alt,
            MAX(REVEL_score) AS REVEL_score_max,
            MIN(REVEL_score) AS REVEL_score_min
        FROM _revel_raw
        GROUP BY _chrom, _pos, _ref, _alt
    """)

    # Join: tx-matched + aggregated
    next_view = f"{current}_revel"
    con.execute(f"""
        CREATE VIEW {next_view} AS
        SELECT c.*,
            tx.REVEL_score AS REVEL_score_tx,
            agg.REVEL_score_max,
            agg.REVEL_score_min
        FROM {current} c
        LEFT JOIN (
            SELECT DISTINCT _chrom, _pos, _ref, _alt, _tx, REVEL_score
            FROM _revel_raw
        ) tx
            ON c."{chrom_col}" = tx._chrom
            AND CAST(c."{pos_col}" AS VARCHAR) = tx._pos
            AND c."{ref_col}" = tx._ref
            AND c."{alt_col}" = tx._alt
            AND c."{tx_col}" = tx._tx
        LEFT JOIN _revel_agg agg
            ON c."{chrom_col}" = agg._chrom
            AND CAST(c."{pos_col}" AS VARCHAR) = agg._pos
            AND c."{ref_col}" = agg._ref
            AND c."{alt_col}" = agg._alt
    """)
    return next_view


def add_popeve(
    con: duckdb.DuckDBPyConnection,
    current: str,
    vcf_path: str,
    chrom: str,
    chrom_col: str,
    pos_col: str,
    ref_col: str,
    alt_col: str,
) -> str:
    """Add popEVE score via DuckDB LEFT JOIN. Returns new view name."""
    chrom_bare = chrom.replace("chr", "")
    next_view = f"{current}_popeve"
    # comment='#' skips both ## metadata and #CHROM header, so provide names
    con.execute(f"""
        CREATE VIEW {next_view} AS
        SELECT c.*, p.popEVE_score
        FROM {current} c
        LEFT JOIN (
            SELECT
                CASE WHEN CHROM LIKE 'chr%' THEN CHROM
                     ELSE 'chr' || CHROM END AS _chrom,
                POS AS _pos,
                REF AS _ref,
                ALT AS _alt,
                TRY_CAST(
                    regexp_extract(INFO, 'popEVE=([^;]+)', 1) AS DOUBLE
                ) AS popEVE_score
            FROM read_csv('{vcf_path}', delim='\t', comment='#',
                 names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'],
                 all_varchar=true)
            WHERE CHROM = '{chrom_bare}' OR CHROM = '{chrom}'
        ) p ON c."{chrom_col}" = p._chrom
            AND CAST(c."{pos_col}" AS VARCHAR) = p._pos
            AND c."{ref_col}" = p._ref
            AND c."{alt_col}" = p._alt
    """)
    return next_view


def add_pext(
    con: duckdb.DuckDBPyConnection,
    current: str,
    pext_path: str,
    chrom: str,
    chrom_col: str,
    pos_col: str,
) -> str:
    """Add pext scores via DuckDB LEFT JOIN. Returns new view name."""
    brain_tissues = [
        "Brain_Amygdala", "Brain_Anteriorcingulatecortex_BA24",
        "Brain_Caudate_basalganglia", "Brain_CerebellarHemisphere",
        "Brain_Cerebellum", "Brain_Cortex", "Brain_FrontalCortex_BA9",
        "Brain_Hippocampus", "Brain_Hypothalamus",
        "Brain_Nucleusaccumbens_basalganglia", "Brain_Putamen_basalganglia",
        "Brain_Spinalcord_cervicalc_1", "Brain_Substantianigra",
    ]
    brain_max_expr = "GREATEST(" + ", ".join(
        f'TRY_CAST("{t}" AS DOUBLE)' for t in brain_tissues
    ) + ")"

    next_view = f"{current}_pext"
    con.execute(f"""
        CREATE VIEW {next_view} AS
        SELECT c.*, pext.pext_mean, pext.pext_brain_max
        FROM {current} c
        LEFT JOIN (
            SELECT
                split_part(locus, ':', 1) AS _chrom,
                split_part(locus, ':', 2) AS _pos,
                MAX(TRY_CAST(exp_prop_mean AS DOUBLE)) AS pext_mean,
                MAX({brain_max_expr}) AS pext_brain_max
            FROM read_csv('{pext_path}', delim='\t', header=true,
                 all_varchar=true)
            WHERE locus LIKE '{chrom}:%'
            GROUP BY 1, 2
        ) pext ON c."{chrom_col}" = pext._chrom
              AND CAST(c."{pos_col}" AS VARCHAR) = pext._pos
    """)
    return next_view


def main():
    parser = argparse.ArgumentParser(
        description="Add variant-level annotations (REVEL, popEVE, pext) to parquet files."
    )
    parser.add_argument("input", help="Input file (parquet or TSV)")
    parser.add_argument("-o", "--output", required=True, help="Output file (parquet)")
    parser.add_argument("--revel", metavar="PATH",
                        help="dbNSFP per-chrom file path with {chrom} placeholder")
    parser.add_argument("--popeve", metavar="PATH", help="popEVE VCF file path")
    parser.add_argument("--pext", metavar="PATH", help="gnomAD pext base-level TSV path")
    parser.add_argument("--chrom-col", default="#CHROM",
                        help="Chromosome column (default: #CHROM)")
    parser.add_argument("--pos-col", default="POS",
                        help="Position column, 1-based (default: POS)")
    parser.add_argument("--ref-col", default="REF",
                        help="Reference allele column (default: REF)")
    parser.add_argument("--alt-col", default="ALT",
                        help="Alternate allele column (default: ALT)")
    parser.add_argument("--tx-col", default="Feature",
                        help="Transcript ID column for tx-matched scores (default: Feature)")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip annotations for columns that already exist")
    args = parser.parse_args()

    if not any([args.revel, args.popeve, args.pext]):
        sys.exit("ERROR: Provide at least one of --revel, --popeve, --pext")

    print(f"Loading {args.input}...", file=sys.stderr)

    con = duckdb.connect()

    # Register input
    p = Path(args.input)
    if p.suffix == ".parquet":
        con.execute(f"CREATE VIEW input AS SELECT * FROM read_parquet('{args.input}')")
    else:
        con.execute(
            f"CREATE VIEW input AS SELECT * FROM read_csv_auto('{args.input}',"
            f" delim='\\t', header=true)"
        )

    # Check existing columns
    existing_cols = {row[0] for row in con.execute("DESCRIBE input").fetchall()}

    # Detect chromosome
    chrom = con.execute(f'SELECT "{args.chrom_col}" FROM input LIMIT 1').fetchone()[0]

    current = "input"
    total_start = time.time()

    if args.revel:
        if args.skip_existing and "REVEL_score_tx" in existing_cols:
            print("SKIP  REVEL (columns already exist)", file=sys.stderr)
        else:
            revel_path = resolve_revel_path(args.revel, chrom)
            if revel_path:
                print("ADD   REVEL (tx, max, min)", file=sys.stderr)
                current = add_revel(
                    con, current, revel_path,
                    args.chrom_col, args.pos_col, args.ref_col, args.alt_col,
                    args.tx_col,
                )

    if args.popeve:
        if args.skip_existing and "popEVE_score" in existing_cols:
            print("SKIP  popEVE (column already exists)", file=sys.stderr)
        else:
            if not Path(args.popeve).exists():
                print(f"  WARNING: popEVE file not found: {args.popeve}",
                      file=sys.stderr)
            else:
                print("ADD   popEVE", file=sys.stderr)
                current = add_popeve(
                    con, current, args.popeve, chrom,
                    args.chrom_col, args.pos_col, args.ref_col, args.alt_col,
                )

    if args.pext:
        if args.skip_existing and "pext_mean" in existing_cols:
            print("SKIP  pext (column already exists)", file=sys.stderr)
        else:
            if not Path(args.pext).exists():
                print(f"  WARNING: pext file not found: {args.pext}",
                      file=sys.stderr)
            else:
                print("ADD   pext", file=sys.stderr)
                current = add_pext(
                    con, current, args.pext, chrom,
                    args.chrom_col, args.pos_col,
                )

    # Write output — all actual work (reads, joins) happens here
    out_path = Path(args.output)
    print(f"Writing {out_path}...", file=sys.stderr)

    if out_path.suffix == ".parquet":
        con.execute(f"COPY {current} TO '{out_path}' (FORMAT PARQUET)")
    else:
        con.execute(f"COPY {current} TO '{out_path}' (HEADER, DELIMITER '\\t')")

    elapsed = time.time() - total_start
    print(f"Done in {elapsed:.1f}s -> {out_path}", file=sys.stderr)
    con.close()


if __name__ == "__main__":
    main()
