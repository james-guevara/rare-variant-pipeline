# Candidate-only missense scoring and HC LoF tier annotation

`candidates.nf` is a separate entrypoint after the completed annotation stage. It
reads published picked and LOFTEE TSVs, creates small candidate Parquets, and never
invokes sites preparation, FastVEP, LOFTEE, or genotype extraction. The existing
exact-HC carrier definition and workflow are unchanged. No external gene T1/T2
lists or gene sets are applied.

```text
published picked.tsv → missense consequence token → exact dbNSFP join → n_flag ≥ 1 → missense.parquet
published loftee.tsv → LoF == HC → GeneBayes join → annotate LoF tier → lof_hc.parquet
```

One selected block is one Nextflow task. Start only `chr22_block12`; no NBDC run
has been launched by Codex. The reported completed chr22 annotation/carrier counts
are operator evidence; new missense/score-match counts must be measured by this pilot.

## Inspected code/resources and scientific behavior

The implementation imports `T_STARS`, `LOF_T1`, and `LOF_T2` from the unchanged
[`tier_variants.py`](../../../scripts/postprocess/tier_variants.py), and `SCORES`
and `extract_expr` from unchanged
[`join_scores.py`](../../../scripts/postprocess/join_scores.py).

| Rankscore | Inclusive threshold |
|---|---:|
| ClinPred_rankscore | 0.4298 |
| AlphaMissense_rankscore | 0.9603 |
| popEVE_converted_rankscore | 0.9209 |
| MPC_rankscore | 0.8947 |

`n_flag` counts passing rankscores; missing/unparseable values do not pass. The
existing mapping is 4 → `miss_t1`, 3 → `miss_t2`, 2 → `miss_t3`, 1 → `miss_t4`.
`miss_n_flag` is also supplied for compatibility. Consequence membership is tested
as an exact token in an `&`-separated list, not an arbitrary substring.

The exact resource keys are `#chr`, `pos(1-based)`, `ref`, `alt`. The join checks
all four, allowing only the `chr22`/`22` naming alias, with no allele normalization.
It uses the actual chromosome column from the resource rather than substituting
the requested chromosome. Missing exact matches produce null scores, zero flags,
and explicit nonmatch counts in the receipt. A match with missing scores remains
a resource match, distinguishable from an absent allele.

All eleven existing score outputs are retained for selected missense candidates:
five rankscores (the four above plus REVEL) and ClinPred, AlphaMissense, REVEL, MPC,
popEVE raw scores plus CADD_phred. Rankscores and scalar scores use the existing
TRY_CAST behavior. Raw semicolon lists use MAX, except popEVE uses MIN. Duplicate
resource allele keys aggregate using those same rules, never multiplying output
candidates; the receipt reports duplicate matched keys. These are the existing
variant-level resource reductions, not a new transcript-specific score match.

GeneBayes joins **`Gene == ensg`**, as the existing implementation actually does.
Its old prose mentions symbol matching; that is not the implemented GeneBayes join.
The six existing `genebayes_*` metrics are retained. Inclusive `post_mean >= 0.18`
gives `lof_t1`; `0.03 <= post_mean < 0.18` gives `lof_t2`. Missing and lower scores
remain untiered. Duplicate GeneBayes gene IDs fail rather than multiply rows. There
is no symbol fallback or automatic gene-ID version stripping.

HC allele/gene/transcript identities must agree with the paired picked file. A
wrong file pairing, duplicate candidate keys, or invalid candidate alleles fails.
A variant whose consequence includes both missense and a LoF consequence can be
represented in both outputs; overlap is reported and the outputs must not be
blindly concatenated for carrier counting. The separate masks do not define a
single combined tier. No genotype scans or new carrier implementation are added;
a future missense carrier adapter should feed these exact keys into the existing
indexed exact-allele machinery, preserving the HC extractor's current contract.

## Output retention and spanning deletions

Published under `<outdir>/candidates/<unit_id>/`:

- `missense.parquet`: only rows with `n_flag >= 1`, their picked annotation columns,
  `allele_class`, all score outputs, resource match/multiplicity fields, and tier.
- `lof_hc.parquet`: **all HC rows**, including unmatched/untiered GeneBayes rows,
  their existing annotations, `allele_class`, GeneBayes metrics/match flag, and tier.
- `receipt.json`: input/resource/code identities and hashes, counts by allele class,
  dbNSFP matches/nonmatches, n_flag distribution including zero, HC/GeneBayes match
  counts, LoF tier counts, output row counts/sizes/hashes, DuckDB version and timing.

`allele_class` is `sequence` or `spanning_deletion`. Literal `ALT=*` remains an exact
key and is not converted to another ALT or merged with sequence alleles. HC stars
are retained even when GeneBayes is missing. Assigned LoF tiers on star rows are
annotation labels, not a conclusion that each star is an independent biological
pLoF. No biological event deduplication or burden calculation occurs here.

**Selection tradeoff:** unmatched missense rows and n_flag=0 rows are counted but
not persisted in the selected missense Parquet. They remain recoverable from the
original picked output and immutable resource. All HC LoF rows are retained because
later gene-based analyses may need HC calls outside these two constraint tiers.
No fully scored copy of the complete picked annotation is written.

DuckDB scans the input TSV through views, materializes only missense/HC subsets,
projects the needed resource columns, and joins to candidate keys before score
aggregation. Parquet is not a point-query index: DuckDB may read the projected
columns across the chromosome resource. It does not create another genome-wide
TSV or copy the resource into each work directory. DuckDB uses 70% of task memory
and can spill temporary data within the work directory.

## Resource inspection and identity

Read-only inspection of the canonical Expanse release confirmed:

- Expanded-MANE chr22 Parquet: 1,615,190 rows, documented four key columns and all
  eleven score source columns present, stored as strings.
- GeneBayes header contains `ensg`, `hgnc`, and the six expected metric columns.
- The already-staged pinned LOFTEE SIF contains DuckDB **1.5.5** (confirmed on Expanse).

Canonical v1 file identities, also reported verified by the operator on NBDC:

| File | Bytes | SHA-256 |
|---|---:|---|
| chr22 expanded-MANE Parquet | 298748270 | `7ca84d9bf85627ab66aa86f3fe91a0e28630a7307a92caa3284589ccb21d50a9` |
| GeneBayes TSV | 1572515 | `d5a88129246bb8a1f157c29d6bb566a81234fc752a3b1dd05fce0a422d7e49f3` |

The lock helper verifies these against the repository's canonical inventory and
checks the pinned SIF hash. It records resolved paths, sizes, timestamps, and hashes
once. Subsequent launches/tasks check metadata and use the immutable identities in
Nextflow cache keys. Resource symlinks escaping the read-only root are rejected.
Keep resources immutable; deliberate replacements require a freshly verified lock.

## Proposed block12 commands on NBDC

Run from the updated repository checkout. This consumes existing published
annotation; do not rerun it or the carrier stage. This script creates a **one-block**
manifest and an identity lock (no annotation/genotype processing):

```bash
set -euo pipefail
BASE=/home/ood-guevara-james/abcd-fastvep-smoke
python3 - "$BASE" > abcd-block12-candidates.tsv <<'PY'
from pathlib import Path
import sys
base = Path(sys.argv[1]) / 'annotation-nextflow'
unit = 'chr22_block12'
picked = base / 'fastvep-picker' / unit / 'picked.tsv'
loftee = base / 'loftee' / unit / 'loftee.tsv'
if not picked.is_file() or not loftee.is_file():
    raise SystemExit('Published annotation files are missing')
print('unit_id\tchromosome\tpicked\tloftee')
print(f'{unit}\tchr22\t{picked}\t{loftee}')
PY

# Existing annotation lock identifies the already-validated LOFTEE SIF.
# Adjust only this path if your original lock is in another launch directory.
python3 scripts/lock_candidate_resources.py \
  --resource-root "$BASE/resources" \
  --chromosomes chr22 \
  --annotation-lock "$PWD/abcd-annotation-resources.json" \
  --output abcd-candidate-resources.json
```

The proposed pilot launch is:

```bash
nextflow -C candidates.config run candidates.nf -profile nbdc \
  --candidate_manifest "$PWD/abcd-block12-candidates.tsv" \
  --candidate_resource_lock "$PWD/abcd-candidate-resources.json" \
  --select_units chr22_block12 \
  --outdir "$BASE/candidate-nextflow" \
  -work-dir "$BASE/candidate-nextflow-work" \
  -with-trace "$BASE/candidate-block12.trace.tsv" \
  -resume
```

Defaults: 2 CPUs, 4 GB memory, 2h, medium partition/QoS, queueSize 4; configurable
with `--candidate_cpus`, `--candidate_memory`, `--candidate_time`,
`--candidate_queue`, `--candidate_qos`, `--candidate_queue_size`. `-C` isolates the
legacy configuration. `conf/candidates-nbdc.config` binds the shared resources root
read-only and uses the existing LOFTEE SIF. Host Python 3.9+ is required for the lock
helper. Workflow tasks use the container's Python/DuckDB.

Aggregate inspection only:

```bash
python3 - "$BASE/candidate-nextflow/candidates/chr22_block12/receipt.json" <<'PY'
import json, sys
r = json.load(open(sys.argv[1]))
assert r['status'] == 'passed'
assert r['picked_input_rows'] == 48001 and r['loftee_input_rows'] == 332
assert r['hc_rows'] == 302
assert r['by_allele_class']['sequence']['hc_rows'] == 203
assert r['by_allele_class']['spanning_deletion']['hc_rows'] == 99
keys = ['missense_input_rows','hc_rows','duplicate_dbnsfp_keys',
        'overlapping_missense_hc_keys','by_allele_class','outputs','wall_seconds']
print(json.dumps({k:r[k] for k in keys}, indent=2))
PY
```

New missense/GeneBayes match counts are intentionally not predicted from the old
carrier counts. Inspect this receipt, then repeat block12 with a new trace filename
and `-resume` to confirm caching before considering other blocks. Keep launch/work
paths stable. No all-chr22 candidate job has been launched.

Failed tasks retain a failed receipt and private log in the work directory;
partial Parquets are removed and not published. Keep Parquet data and private
logs on NBDC; only aggregate receipts should be shared. Existing successful files
may remain after a failed rerun, so require a successful Nextflow exit and matching
receipt before using outputs.

## Tests and limits

Synthetic tests exercise actual DuckDB joins and Nextflow tasks, including exact
chromosome/position/REF/ALT matching, missing resource annotations, threshold
boundaries, duplicate resource keys, missing GeneBayes matches, stars, empty inputs,
wrong input pairing, canonical hash verification, caching, and durable publication.
A regression invokes the original `tier_variants.py` and compares its tier labels
against both new products. These are not real NBDC block12 results.

Local validation: **38 tests passed** (8 candidate tests plus the existing exact
carrier, sites annotation, and sites catalog suites):

```bash
uv run --with pytest --with duckdb==1.5.5 --with pysam python -m pytest -q \
  tests/test_scored_candidates.py tests/test_exact_carriers.py \
  tests/test_sites_annotation.py tests/test_sites_catalog.py
```

The NBDC profile also passed `nextflow -C candidates.config config -profile nbdc
-flat`. Real block12 scoring and its new aggregate counts remain pending operator
execution on NBDC.
