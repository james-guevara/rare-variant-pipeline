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

## MANE policy: transcript preference, never variant inclusion (corrected 2026-10-02)

**Use `parquet_expanded`, not `parquet_expanded_mane_select`.** Earlier versions of
this PR incorrectly used a reference already filtered to MANE-containing rows.
The initial brief named that filtered resource, but the current v3 policy and the
operator's clarification require no MANE inclusion filter. Removing a filter on
picked annotations alone cannot recover rows removed during resource preparation.

The earlier April notes and
[historical per-variant analysis](https://github.com/james-guevara/rare-variant-postprocessing/blob/main/steps/missense_counts_pervariant.py)
used the filtered resource. They do not override the later v3 policy. The canonical
v3 notes explicitly remove MANE inclusion filtering and the older SHANK3 exception.
The March investigation also distinguished off-MANE variants in genes that do
have MANE transcripts from genes with no MANE transcript. Neither case is grounds
for exclusion here. No gene whitelist or MANE/canonical fallback filter is applied.

The existing
[resource builder](https://github.com/james-guevara/rare-variant-postprocessing/blob/main/scripts/regen_dbnsfp_expanded.py)
already produces the required unfiltered expanded resource. The score reducers
remain unchanged: use raw-score MAX across transcript lists (MIN for popEVE), and
the supplied rankscores for the four flags. This retains MPC values stored at a
non-MANE transcript position. Historical whole-chr22 MPC coverage was 53.0% at the
MANE position versus 73.5% with list-max; these are prior measurements, not ABCD
results and not the separate cohort-specific coverage denominator.

MANE remains a preference in the existing upstream transcript picker; this stage
neither repicks transcripts nor filters on picked MANE status. HC/GeneBayes has no
MANE restriction either. A gene or allele absent from dbNSFP, or with genuinely
missing predictor values, is not assigned fabricated scores. Unscored missense
rows remain counted but do not acquire a tier. Changing resource coverage can
change candidate counts; old MANE-filtered candidate results are not equivalent.

Existing alternatives were inspected read-only at:

```text
/expanse/projects/sebat1/s3/data/sebat/resources/dbNSFP/5.3.1a/
```

| Product | Existing semantics | Used here? |
|---|---|---|
| `parquet_expanded` | Expanded predictor columns without a MANE row filter | **Yes** |
| `parquet_expanded_mane_select` | Rows with `MANE LIKE '%Select%'`, transcript lists intact | **No; superseded in this PR** |
| `parquet_mpc` | Gene-keyed preferred row, MANE and maximum MPC, transcript IDs, `has_mane_select`; retains a non-MANE row when needed | No |
| `parquet_am` | Analogous MANE/maximum AlphaMissense representation | No |
| `parquet_scores_af` | Earlier scores and AF; builder lacks ClinPred | No |
| `parquet_af` | Allele-frequency columns | No |

`build_dbnsfp_split_parquets.py` beside those resources builds the specialized
MPC/AlphaMissense products. They are not needed for this pilot because the
unfiltered expanded product contains the required score columns. Their earlier
existence is not evidence that a different resource was used in every later run.

`n_scored` counts available rankscores among the four predictors; `n_flag` counts
threshold passes. Receipts include both distributions, each predictor's missing
count, matched rows with no rankscores, and fully scored rows below all thresholds.
Missing scores are not evidence of benignity. Joins remain exact CHROM/POS/REF/ALT;
gene/transcript identities are not added to the score join in this correction.
Duplicate resource rows still reduce independently per score. Scores are therefore
not asserted to be transcript-matched to the picked consequence.

## Output retention and spanning deletions

Published under `<outdir>/candidates/<unit_id>/`:

- `missense.parquet`: only rows with `n_flag >= 1`, their picked annotation columns,
  `allele_class`, all score outputs, `n_scored`, resource match/multiplicity fields, and tier.
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

Read-only Expanse inspection on 2026-10-02 confirmed the unfiltered chr22
Parquet contains **1,799,075 rows** and all eleven score columns. The old filtered
file contains 1,615,190 rows. These are resource row counts, not ABCD candidate counts.
GeneBayes and the pinned LOFTEE SIF (DuckDB 1.5.5) are unchanged.

| File | Bytes | SHA-256 |
|---|---:|---|
| **Unfiltered** chr22 expanded Parquet | 310791513 | `ecb7a6594c7db0c4b8d0ebdf9b1e8a152f633dc3d0c282ccd54b50eb3d6dbfb1` |
| GeneBayes TSV | 1572515 | `d5a88129246bb8a1f157c29d6bb566a81234fc752a3b1dd05fce0a422d7e49f3` |

The [candidate inventory](../../resources/candidate-resource-hashes.json) now
pins unfiltered chr1–22/X/Y and the unchanged GeneBayes identity, imported from
the completed consolidation in PR #14. These files are in the canonical v1
Expanse/S3 roots; historical MANE-filtered files are retained. Execution-host
staging and lock coverage are separate from canonical availability. See the
[autosome readiness guide](AUTOSOME_READINESS.md). There is no chromosome-specific
scientific logic.

The lock helper verifies SHA-256, records metadata and the pinned SIF identity,
and emits schema 2 with `dbnsfp_representation=parquet_expanded`. Old schema-1
MANE-filtered locks are rejected before tasks run. Rebuild the lock; replacing or
renaming the old filtered file will fail its expected hash. Resources stay read-only
and are not copied into individual task work directories.

**One new file must be staged on NBDC.** Codex verified the source on Expanse but
has not transferred it or accessed NBDC:

```text
Source:
/expanse/projects/sebat1/s3/data/sebat/resources/dbNSFP/5.3.1a/parquet_expanded/chr22.parquet

Destination:
/home/ood-guevara-james/abcd-fastvep-smoke/resources/dbNSFP/5.3.1a/parquet_expanded/chr22.parquet
```

Use the established permitted transfer route between your environments. Keep the
old filtered resource separately; do not overwrite it. After transfer, run on NBDC:

```bash
BASE=/home/ood-guevara-james/abcd-fastvep-smoke
printf '%s  %s\n' \
  ecb7a6594c7db0c4b8d0ebdf9b1e8a152f633dc3d0c282ccd54b50eb3d6dbfb1 \
  "$BASE/resources/dbNSFP/5.3.1a/parquet_expanded/chr22.parquet" | sha256sum -c -
```

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

After the initial chromosome runs, review the pipeline with the operator stage
by stage against actual commands and resource builders: inputs, every inclusion
and exclusion rule, transcript/gene assignment, missing-value handling, joins,
and retained outputs. Passing regression tests does not by itself validate those
scientific choices. This review is planned after those runs, not a requirement to
repeat already-completed preparation or annotation now.

Synthetic tests exercise actual DuckDB joins and Nextflow tasks, including exact
chromosome/position/REF/ALT matching, missing resource annotations, threshold
boundaries, duplicate resource keys, missing GeneBayes matches, stars, empty inputs,
wrong input pairing, canonical hash verification, caching, and durable publication.
A regression invokes the original `tier_variants.py` and compares its tier labels
against both new products. Additional regressions cover genes with no MANE,
variants off a gene's MANE transcript, non-MANE MPC scores, genuinely unscored
non-MANE alleles, and rejection of old filtered locks. These are not real NBDC
block12 results.

Local validation: **42 tests passed** (12 candidate tests plus the existing exact
carrier, sites annotation, and sites catalog suites):

```bash
uv run --with pytest --with duckdb==1.5.5 --with pysam python -m pytest -q \
  tests/test_scored_candidates.py tests/test_exact_carriers.py \
  tests/test_sites_annotation.py tests/test_sites_catalog.py
```

The NBDC profile also passed `nextflow -C candidates.config config -profile nbdc
-flat`. Real block12 scoring and its new aggregate counts remain pending operator
execution on NBDC.
