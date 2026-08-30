# Candidate gene-set catalog

This directory deliberately preserves multiple defensible gene-set definitions rather than
selecting one prematurely. Raw, version-pinned sources live under `raw/<retrieval-date>/`;
derived flat tables live under `processed/<retrieval-date>/`.

Rebuild with:

```bash
bash scripts/download_gene_set_sources.sh
uv run --script scripts/build_gene_set_catalog.py
```

The processed release contains:

- `gene_set_membership.tsv`: one row per gene/set membership, including source and evidence.
- `gene_set_summary.tsv`: definition and size of every nonempty set.
- `gene_set_overlap.tsv`: pairwise intersection, union, and Jaccard similarity.
- `gene_set_sources.tsv`: source release, URL, retrieval date, and usage notes.
- `schema_release_comparison.tsv`: gene-level shared/gained/lost calls for published-versus-current
  and consecutive Phase II Bonferroni sets.

Current alternatives include current SFARI score strata, broad and green PanelApp ID/epilepsy
panels, the published SCHEMA-10, six Phase II SCHEMA snapshots at both Bonferroni and 5% FDR,
Epi25 consequence/phenotype/significance strata, the 207-gene Kosmicki height set and signed
singleton-pLoF subsets, and signed exome-wide BMI sets.

Important: SCHEMA Phase II carries time-sensitive consortium usage terms. Review the current
terms at <https://schema.broadinstitute.org/terms> before publishing analyses based on its global,
exome-wide, or large gene-set results. The catalog records candidates; it does not constitute a
scientific or publication-policy decision.
