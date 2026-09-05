# G2MH chr22 comparison with the previous pipeline

Comparison performed 2026-08-28 between:

- previous pipeline: `chr22.filtered_annotated.parquet` from the G2MH
  `rv_postprocess` run of 2026-08-12;
- targeted workflow: `chr22-nextflow-clean-89470da-v2`.

The comparison restricted the previous table to final tiered rows and compared
`chromosome, position, REF, ALT, sample, tier` against the targeted workflow's
burden-eligible LoF and missense tables.

| Tier | Previous rows | Targeted rows |
|---|---:|---:|
| `lof_t1` | 4 | 4 |
| `lof_t2` | 34 | 34 |
| `miss_t1` | 3 | 3 |
| `miss_t2` | 16 | 16 |
| `miss_t3` | 82 | 82 |
| `miss_t4` | 1,010 | 1,010 |

All 1,149 allele-sample-tier rows matched exactly: zero rows were exclusive to
either workflow.

One qualifying LoF insertion had an additional Ensembl 115/FastVEP consequence
subcategory, `splice_polypyrimidine_tract_variant`; its primary frameshift and
splice-region consequences, LoF tier, carrier, and burden result were unchanged.

Four multiallelic missense carriers had equivalent dosage but different serialized
GT values. The previous decomposed table used biallelic encodings (`0/1` or `1/1`),
whereas the targeted table retained source allele index 2 (`0/2` or `2/2`). The
extractor now emits biallelic `GT`/`genotype` and retains the original encoding in
`source_GT`. This changes no carrier or burden counts.
