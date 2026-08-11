# archive/

Code that is no longer part of the pipeline, kept rather than deleted so it stays
findable. Nothing here is referenced by `main.nf`, any subworkflow, any module, or
`nextflow.config` — verified by grep, not by assumption.

## modules/

| File | Why it's here |
|---|---|
| `qc_filter.nf` (`APPLY_QC_FILTER`) | Dropped from the DAG in `747e4d6` (2026-05-01) but the process file and its `withName:` selector lingered, emitting a WARN on every run. Selector and `params.qc_filter_script` removed. |
| `reformat_variants.nf` (`REFORMAT_VARIANTS`) | Superseded by `PREPARE_VARIANTS`, which does the same job without multi-allelic resolution (sites arrive biallelic from `bcftools norm -m-`). |

## scripts/

| File | Why it's here |
|---|---|
| `reformat_variants.py` | Superseded by `prepare_variants.py`. |
| `join_filter.py` | Wrote FILTER onto merged output after the fact. Unnecessary: `SPLIT_VEP` emits `%FILTER` and it survives to the parquet (verified — column 16 of 58). |
| `qc_filter.py` | Driven only by the archived `APPLY_QC_FILTER`. |
| `family_query.py` | Never wired in; `CARRIER_QUERY` does carrier extraction in pure bcftools. |
| `alphagenome_score.py`, `annotate_regions.py`, `annotate_variants.py`, `count_per_sample.py`, `filter_lof.py`, `filter_missense.py`, `reshape_trios.py`, `verify_and_gather.py` | Standalone analysis utilities, not pipeline stages. Kept in case they are still wanted; they were never invoked by the workflow. |

## Restoring something

    git mv archive/scripts/foo.py scripts/

and re-add whatever `params.*` entry it needs. Check `git log --follow` for its history.
