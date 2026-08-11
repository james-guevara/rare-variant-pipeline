#!/usr/bin/env python3
"""
Score rare regulatory variants with AlphaGenome API.

Reads a regulatory parquet file, filters to a targeted subset of rare
noncoding variants, queries AlphaGenome for variant effect predictions
(focused on brain tissues), and outputs an annotated parquet/TSV.

Filtering strategy (applied in order):
  1. gnomAD AF < threshold (default 0.001) or missing
  2. In at least one brain-relevant regulatory track (abc, psychencode, chromhmm)
  3. Optionally: in a specific gene list (via --genes or --gene-file)
     Matches against abc_target_gene column.
  4. Optionally: CADD phred >= threshold

Usage:
    # Score all rare brain-regulatory variants
    python alphagenome_score.py regulatory.parquet -o scored.parquet \
        --api-key YOUR_KEY

    # Score only variants near NDD genes
    python alphagenome_score.py regulatory.parquet -o scored.parquet \
        --api-key YOUR_KEY --gene-file ndd_genes.txt

    # Dry run: just filter and report count (no API calls)
    python alphagenome_score.py regulatory.parquet -o scored.parquet --dry-run

    # Custom AF threshold and CADD filter
    python alphagenome_score.py regulatory.parquet -o scored.parquet \
        --api-key YOUR_KEY --max-af 0.0001 --min-cadd 20
"""
import argparse
import sys
import time
from pathlib import Path

import pandas as pd

# Brain-related ontology terms for AlphaGenome queries (UBERON + CL)
BRAIN_ONTOLOGY_TERMS = [
    "UBERON:0000955",  # brain
    "UBERON:0001890",  # forebrain
    "UBERON:0001950",  # neocortex
    "UBERON:0002298",  # cerebral cortex
    "UBERON:0001870",  # frontal cortex
    "UBERON:0002421",  # hippocampus
    "UBERON:0001876",  # amygdala
    "UBERON:0001873",  # caudate nucleus
    "UBERON:0001874",  # putamen
    "UBERON:0002037",  # cerebellum
    "CL:0000540",      # neuron
    "CL:0000127",      # astrocyte
    "CL:0000128",      # oligodendrocyte
]

# Sequence context length for AlphaGenome predictions
# Longer = more context = better predictions, but slower
SEQUENCE_LENGTH = 100_000  # 100kb, good balance of speed and accuracy


def filter_variants(input_path: str, args) -> pd.DataFrame:
    """Filter parquet to targeted rare regulatory variants."""
    import duckdb

    con = duckdb.connect()

    # Build WHERE clauses
    conditions = []

    # AF filter: rare or missing
    af_col = args.af_col
    conditions.append(
        f'("{af_col}" IS NULL OR "{af_col}" = \'.\' OR '
        f'TRY_CAST("{af_col}" AS DOUBLE) < {args.max_af})'
    )

    # Brain regulatory track filter
    if not args.no_brain_filter:
        brain_tracks = []
        for col in ["in_abc", "in_psychencode", "in_chromhmm"]:
            brain_tracks.append(f'("{col}" = \'1\' OR "{col}" = \'true\')')
        conditions.append(f"({' OR '.join(brain_tracks)})")

    # CADD filter
    if args.min_cadd:
        conditions.append(
            f'(TRY_CAST("CADD_phred" AS DOUBLE) >= {args.min_cadd})'
        )

    where_clause = " AND ".join(conditions)

    # Deduplicate to unique variant sites (chrom/pos/ref/alt)
    query = f"""
        SELECT DISTINCT
            "#CHROM" AS chrom,
            CAST("POS" AS INTEGER) AS pos,
            "REF" AS ref,
            "ALT_specific" AS alt,
            "abc_target_gene",
            "abc_score",
            "CADD_phred",
            "{af_col}" AS gnomad_af,
            "in_abc", "in_psychencode", "in_chromhmm",
            "in_phastcons", "in_ccre",
            "Consequence", "SYMBOL", "IMPACT"
        FROM read_parquet('{input_path}')
        WHERE {where_clause}
    """

    df = con.execute(query).fetchdf()
    con.close()

    # Gene list filter (applied in pandas for flexibility)
    if args.gene_file or args.genes:
        genes = set()
        if args.gene_file:
            with open(args.gene_file) as f:
                for line in f:
                    gene = line.strip()
                    if gene and not gene.startswith("#"):
                        genes.add(gene)
        if args.genes:
            genes.update(args.genes)

        print(f"  Gene list: {len(genes)} genes", file=sys.stderr)
        df = df[df["abc_target_gene"].isin(genes) | df["SYMBOL"].isin(genes)]

    return df


def score_variants(df: pd.DataFrame, api_key: str, args) -> pd.DataFrame:
    """Query AlphaGenome API for each variant and extract brain scores."""
    from alphagenome.data import genome
    from alphagenome.models import dna_client

    model = dna_client.create(api_key)

    results = []
    n = len(df)
    t0 = time.time()

    for i, row in df.iterrows():
        idx = len(results) + 1
        chrom = row["chrom"]
        pos = int(row["pos"])
        ref = row["ref"]
        alt = row["alt"]

        print(
            f"  [{idx}/{n}] {chrom}:{pos} {ref}>{alt}",
            end="",
            flush=True,
            file=sys.stderr,
        )

        half_len = SEQUENCE_LENGTH // 2
        interval = genome.Interval(
            chromosome=chrom,
            start=max(0, pos - half_len),
            end=pos + half_len,
        )
        variant = genome.Variant(
            chromosome=chrom,
            position=pos,
            reference_bases=ref,
            alternate_bases=alt,
        )

        try:
            from alphagenome.models import variant_scorers

            scorers = variant_scorers.get_recommended_scorers(
                organism=genome.Organism.HUMAN
            )
            scores = model.score_variant(
                interval=interval,
                variant=variant,
                variant_scorers=scorers,
            )
            tidy = variant_scorers.tidy_scores(scores)

            # Filter to brain-related tracks
            brain_mask = tidy["ontology_term"].isin(BRAIN_ONTOLOGY_TERMS) | tidy[
                "description"
            ].str.contains("brain|neuron|cortex|cerebel", case=False, na=False)

            brain_scores = tidy[brain_mask]

            # Summarize: max absolute score across brain tracks
            summary = {
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "ag_n_brain_tracks": len(brain_scores),
                "ag_max_abs_score": brain_scores["score"].abs().max()
                if len(brain_scores) > 0
                else None,
                "ag_max_quantile": brain_scores["quantile_score"].max()
                if "quantile_score" in brain_scores.columns
                and len(brain_scores) > 0
                else None,
                "ag_top_track": brain_scores.loc[
                    brain_scores["score"].abs().idxmax(), "description"
                ]
                if len(brain_scores) > 0
                else None,
                "ag_top_output_type": brain_scores.loc[
                    brain_scores["score"].abs().idxmax(), "output_type"
                ]
                if len(brain_scores) > 0
                else None,
            }
            results.append(summary)
            elapsed = time.time() - t0
            rate = idx / elapsed if elapsed > 0 else 0
            print(f" -> {summary['ag_max_abs_score']:.3f} ({rate:.1f}/s)", file=sys.stderr)

        except Exception as e:
            err_msg = str(e)
            if "Rate limit" in err_msg or "quota" in err_msg.lower():
                print(f" -> RATE LIMITED, waiting 60s...", file=sys.stderr)
                time.sleep(60)
                # Retry once
                try:
                    scores = model.score_variant(
                        interval=interval,
                        variant=variant,
                        variant_scorers=scorers,
                    )
                    tidy = variant_scorers.tidy_scores(scores)
                    brain_mask = tidy["ontology_term"].isin(
                        BRAIN_ONTOLOGY_TERMS
                    ) | tidy["description"].str.contains(
                        "brain|neuron|cortex|cerebel", case=False, na=False
                    )
                    brain_scores = tidy[brain_mask]
                    summary = {
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "ag_n_brain_tracks": len(brain_scores),
                        "ag_max_abs_score": brain_scores["score"].abs().max()
                        if len(brain_scores) > 0
                        else None,
                        "ag_max_quantile": brain_scores["quantile_score"].max()
                        if "quantile_score" in brain_scores.columns
                        and len(brain_scores) > 0
                        else None,
                        "ag_top_track": None,
                        "ag_top_output_type": None,
                    }
                    results.append(summary)
                    print(f" -> {summary['ag_max_abs_score']:.3f} (retry ok)", file=sys.stderr)
                except Exception as e2:
                    print(f" -> FAILED after retry: {e2}", file=sys.stderr)
                    results.append({
                        "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                        "ag_n_brain_tracks": None, "ag_max_abs_score": None,
                        "ag_max_quantile": None, "ag_top_track": None,
                        "ag_top_output_type": None,
                    })
            else:
                print(f" -> ERROR: {e}", file=sys.stderr)
                results.append({
                    "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                    "ag_n_brain_tracks": None, "ag_max_abs_score": None,
                    "ag_max_quantile": None, "ag_top_track": None,
                    "ag_top_output_type": None,
                })

    elapsed = time.time() - t0
    print(f"\nScored {len(results)} variants in {elapsed:.0f}s", file=sys.stderr)

    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(
        description="Score rare regulatory variants with AlphaGenome API."
    )
    parser.add_argument("input", help="Input regulatory parquet file")
    parser.add_argument("-o", "--output", required=True, help="Output file (parquet or tsv)")
    parser.add_argument("--api-key", help="AlphaGenome API key (or set ALPHAGENOME_API_KEY env var)")
    parser.add_argument("--max-af", type=float, default=0.001,
                        help="Max gnomAD AF (default: 0.001)")
    parser.add_argument("--af-col", default="gnomAD4.1_joint_AF",
                        help="AF column name (default: gnomAD4.1_joint_AF)")
    parser.add_argument("--min-cadd", type=float, default=None,
                        help="Min CADD phred score (default: no filter)")
    parser.add_argument("--genes", nargs="+", help="Gene symbols to filter on")
    parser.add_argument("--gene-file", help="File with one gene symbol per line")
    parser.add_argument("--no-brain-filter", action="store_true",
                        help="Don't require brain regulatory track overlap")
    parser.add_argument("--dry-run", action="store_true",
                        help="Filter only, report count, no API calls")
    parser.add_argument("--save-filtered", metavar="PATH",
                        help="Save filtered variant list before scoring")
    args = parser.parse_args()

    import os
    api_key = args.api_key or os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key and not args.dry_run:
        sys.exit("ERROR: Provide --api-key or set ALPHAGENOME_API_KEY env var")

    # Step 1: Filter
    print("Filtering variants...", file=sys.stderr)
    df = filter_variants(args.input, args)
    print(f"  {len(df)} variants after filtering", file=sys.stderr)

    if len(df) == 0:
        print("No variants passed filters. Nothing to score.", file=sys.stderr)
        return

    if args.save_filtered:
        p = Path(args.save_filtered)
        if p.suffix == ".parquet":
            df.to_parquet(p, index=False)
        else:
            df.to_csv(p, sep="\t", index=False)
        print(f"  Saved filtered variants to {p}", file=sys.stderr)

    if args.dry_run:
        print(f"\nDry run: {len(df)} variants would be scored.", file=sys.stderr)
        print("\nTop genes (by abc_target_gene):", file=sys.stderr)
        if "abc_target_gene" in df.columns:
            counts = df["abc_target_gene"].value_counts().head(20)
            for gene, n in counts.items():
                print(f"  {gene}: {n}", file=sys.stderr)
        return

    # Step 2: Score
    print(f"\nScoring {len(df)} variants with AlphaGenome...", file=sys.stderr)
    scores_df = score_variants(df, api_key, args)

    # Step 3: Merge scores back with filtered variant info
    merged = df.merge(scores_df, on=["chrom", "pos", "ref", "alt"], how="left")

    # Step 4: Write output
    out_path = Path(args.output)
    if out_path.suffix == ".parquet":
        merged.to_parquet(out_path, index=False)
    else:
        merged.to_csv(out_path, sep="\t", index=False)

    print(f"Done -> {out_path} ({len(merged)} rows)", file=sys.stderr)

    # Summary
    scored = merged["ag_max_abs_score"].notna().sum()
    print(f"  Successfully scored: {scored}/{len(merged)}", file=sys.stderr)
    if scored > 0:
        print(f"  Max absolute score: {merged['ag_max_abs_score'].max():.4f}", file=sys.stderr)
        print(f"  Mean absolute score: {merged['ag_max_abs_score'].mean():.4f}", file=sys.stderr)


if __name__ == "__main__":
    main()
