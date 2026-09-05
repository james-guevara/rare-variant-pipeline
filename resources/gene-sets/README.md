# Candidate gene-set catalog

This directory deliberately preserves multiple defensible gene-set definitions rather than
selecting one prematurely. Raw, version-pinned sources live under `raw/<retrieval-date>/`;
derived flat tables live under `processed/<retrieval-date>/`.
Exact historical pipeline lists are retained under `legacy/2026-01-15/`.

Rebuild with:

```bash
bash scripts/download_gene_set_sources.sh
uv run --script scripts/build_gene_set_catalog.py
```

The ignored raw historical bundle can be recovered from Expanse, together with the primary
Fu and Satterstrom publication supplements, using `scripts/recover_legacy_gene_set_sources.sh`.

The processed release contains:

- `gene_set_membership.tsv`: one row per gene/set membership, including source and evidence.
- `gene_set_summary.tsv`: definition and size of every nonempty set.
- `gene_set_overlap.tsv`: pairwise intersection, union, and Jaccard similarity.
- `gene_set_sources.tsv`: source release, URL, retrieval date, and usage notes.
- `schema_release_comparison.tsv`: gene-level shared/gained/lost calls for published-versus-current
  and consecutive Phase II Bonferroni sets.

For sharing, `share/2026-08-29/all_gene_sets.csv` is a denormalized CSV containing every gene-set
membership together with its set size, definition, release, source URL, and usage note.
`share/2026-08-29/gene_sets_wide.csv` is the compact companion table: one column per gene set,
with its actual gene symbols beneath it. `share/2026-08-29/gene_set_summary.csv` has one row per
gene set with its phenotype, size, definition, source, release, URL, and usage note. Regenerate all
three with
`python3 scripts/export_gene_set_catalog_csv.py`.

Current alternatives include current SFARI score strata, broad and green PanelApp ID/epilepsy
panels, the published SCHEMA-10, six Phase II SCHEMA snapshots at both Bonferroni and 5% FDR,
Epi25 consequence/phenotype/significance strata, the 207-gene Kosmicki height set and signed
singleton-pLoF subsets, and signed exome-wide BMI sets.

Important: SCHEMA Phase II carries time-sensitive consortium usage terms. Review the current
terms at <https://schema.broadinstitute.org/terms> before publishing analyses based on its global,
exome-wide, or large gene-set results. The catalog records candidates; it does not constitute a
scientific or publication-policy decision.
