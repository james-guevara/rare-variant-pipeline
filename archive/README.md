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
| `filter_pass_and_rescue.py` | Removed from the POSTPROCESS chain deliberately. It did two things, and both belong elsewhere: (1) a `FILTER==PASS` cut — FILTER is now carried as a column so the cut can be applied downstream where it is visible; (2) a de-novo rescue of non-PASS rows — that can only fire in cohorts with trios, so the same variant would be filtered differently between cohorts with and without them, biasing cross-cohort comparison. De novo status belongs in a later annotation step, not a filter. `dnm_paths` is left in `resources.json` as a pointer to the callsets that step will need. |
| `family_query.py` | Never wired in; `CARRIER_QUERY` does carrier extraction in pure bcftools. |
| `annotate_variants.py`, `annotate_regions.py`, `filter_lof.py`, `filter_missense.py`, `count_per_sample.py` | **Superseded POSTPROCESSING, not dead code.** An earlier, partial generation of what now lives in `scripts/postprocess/` and runs as the `POSTPROCESS` subworkflow — e.g. `annotate_variants.py` joins REVEL/popEVE/pext where `join_scores.py` joins full dbNSFP. Kept for reference only. |
| `alphagenome_score.py`, `reshape_trios.py`, `verify_and_gather.py` | Standalone utilities, never invoked by the workflow. |

## Restoring something

    git mv archive/scripts/foo.py scripts/

and re-add whatever `params.*` entry it needs. Check `git log --follow` for its history.
