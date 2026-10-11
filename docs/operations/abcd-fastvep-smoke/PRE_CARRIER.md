# Pre-carrier AF and region filtering

`pre_carrier.nf` consumes completed variant-level `missense.parquet` and
`lof_hc.parquet`. It writes new files and never invokes annotation, candidate
selection, or genotype extraction. The pipeline is:

```text
sites → annotation → score/select → pre-carrier filter → carrier extraction
```

## Exact policy

- Cohort: exact allele match to the per-block **sites-only** VCF; compute
  `INFO/AC / INFO/AN < 0.005` for preliminary eligibility (uncorrected). Never use INFO/AF or MLEAF for selection. Original INFO AC/AN/AF strings are retained in `pcf_source_info_*` alongside the AC/AN ratio. Missing exact matches,
  absent/malformed AC/AN, zero AN, negative values, or AC > AN fail the cohort
  criterion and are counted. AC and AN must be single nonnegative integers.
- Population: exact allele join to `parquet_scores_af`, using only
  `gnomAD4.1_joint_POPMAX_AF`. NULL/`.`/empty or absent allele passes; otherwise
  require `< 0.001`. Joint AF is not used. Malformed nonmissing values fail the
  task rather than silently becoming rare. Duplicate identical values collapse;
  conflicting nonmissing values fail rather than invent an aggregation policy.
- Regions: exclude POS overlaps in genomicSuperDups, simpleRepeat, or the
  `Simple_repeat`/`Low_complexity` classes of rmsk. The existing BED convention is
  **start < POS <= end**. This is a POS-only test, not a REF-length interval test.
  BED columns 1–3 are chromosome/start/end; rmsk repClass is **column 7**, as in
  the existing `scripts/postprocess/filter_regions.py`. Other rmsk classes do not
  exclude variants. `chr22`/`22` aliasing applies to all joins/tracks.
- There is no MANE filter, gene filter, genotype-quality filter, or change to the
  functional candidate definitions. Literal `ALT=*` is matched exactly and stays
  classified as `spanning_deletion`; there is no biological event deduplication.

Every candidate is evaluated against all three criteria. Receipt `independent`
counts use **all input candidates** as their denominator; `sequential` counts
show survivors after cohort AF, then POPmax, then regions. Independent failures
and track overlaps can overlap and must not be added to infer total exclusions.

No genotype VCF input exists in the command interface. An eight-column sites-only
header is mandatory and a genotyped header fails before any data record is read
or the file is hashed. The stage streams the small sites VCF once to retrieve
candidate AC/AN, then hashes it for provenance. It never creates a genotype VCF.

## Outputs and execution

One independent tagged task per selected unit publishes copies under:

```text
<outdir>/pre-carrier/<unit_id>/
  missense.filtered.parquet
  lof_hc.filtered.parquet
  filter_audit.parquet
  receipt.json
```

Filtered Parquets preserve original columns and add `pcf_` fields for AC, AN,
cohort AF, POPmax, matches, track overlaps, filter passes and final retention.
`filter_audit.parquet` holds **all input candidate keys**, candidate type,
allele class, and those audit fields, including excluded candidates. It omits
redundant annotation/score columns. This is candidate-sized, not a copy of the
full annotation catalog. Audit rows from the two types are not a disjoint burden.

Receipts report per-type inputs, all filter pass/fail counts, sequential and final
counts, missing cohort/site/POPmax values, each track's overlaps and union,
sequence/spanning-deletion counts, resource identities, input/output SHA-256,
script hash, DuckDB version and wall time. Candidate input hashes are checked
before and after filtering. No individual records are printed in aggregate logs.
Failed tasks retain a failure receipt in work, remove partial new Parquets, and
publish no new products. Require a successful task/receipt; older published files
can remain after a later failed attempt.

The existing candidate files are read-only inputs. `select_scored_candidates.py`
and the HC carrier workflow are unchanged. Shared resources are bound read-only;
they are not copied into each task. Defaults are 2 CPUs/4 GB/2h, Slurm medium
partition and QoS, queueSize 4, configurable via `--filter_cpus`, `--filter_memory`,
`--filter_time`, `--filter_queue`, `--filter_qos`, `--filter_queue_size`.

## NBDC block12 pilot

Use the repository branch containing this workflow. This is a **separate** entrypoint
and separate resource lock. The existing pinned LOFTEE SIF contains DuckDB 1.5.5.
No new dependency installation is required inside the image.

Required shared files on NBDC:

```text
resources/dbNSFP/5.3.1a/parquet_scores_af/chr22.parquet
resources/problematic-regions/genomicSuperDups.bed
resources/problematic-regions/simpleRepeat.bed
resources/problematic-regions/rmsk.bed
```

The operator reports the three BEDs staged. Presence of the `parquet_scores_af`
file has not been checked by Codex on NBDC. The expanded scoring Parquet is **not**
a substitute. If missing, stage the score/AF product documented in the
[resource inventory](../../resources/README.md). The helper below fails if a
required file is absent; it records full hashes of the **actual staged resources**
and verifies the SIF against its previously validated SHA-256. It does not claim
to independently authenticate the provenance of operator-provided BED files.

From the repository checkout on NBDC:

```bash
set -euo pipefail
BASE=/home/ood-guevara-james/abcd-fastvep-smoke
SIF="$BASE/containers/targeted-81e0e00a97f900867808e5f7cf17cdac6f5c5ef7956e2c6641f3188fa335d9ab.sif"

python3 scripts/lock_pre_carrier_resources.py \
  --resource-root "$BASE/resources" --chromosomes chr22 \
  --container "$SIF" --output "$PWD/abcd-pre-carrier-resources.json"

python3 - "$BASE" > abcd-block12-pre-carrier.tsv <<'PY'
from pathlib import Path
import sys
base=Path(sys.argv[1]);unit='chr22_block12'
candidate=base/'candidate-nextflow'/'candidates'/unit
paths=[candidate/'missense.parquet',candidate/'lof_hc.parquet',
       base/'sites-catalog-smoke'/'sites'/(unit+'.sites.vcf.gz')]
if not all(p.is_file() for p in paths):
    raise SystemExit('Missing published candidate/site input; adjust manifest paths')
print('unit_id\tchromosome\tmissense\tlof_hc\tsites')
print('\t'.join([unit,'chr22',*map(str,paths)]))
PY

nextflow -C pre_carrier.config run pre_carrier.nf -profile nbdc \
  --filter_manifest "$PWD/abcd-block12-pre-carrier.tsv" \
  --filter_resource_lock "$PWD/abcd-pre-carrier-resources.json" \
  --select_units chr22_block12 \
  --outdir "$BASE/pre-carrier-nextflow" \
  -work-dir "$BASE/pre-carrier-nextflow-work" \
  -with-trace "$BASE/pre-carrier-block12.trace.tsv" -resume
```

Adapt manifest paths if completed candidates were published elsewhere. For other
blocks, add explicit rows; do not infer block IDs from chromosome alone. After
inspecting the pilot, repeat with a new trace filename and the same work/launch
paths to confirm `CACHED`. Selecting all manifest rows uses `--select_units all`.
No real NBDC task has been launched by Codex.

Aggregate inspection:

```bash
python3 - "$BASE/pre-carrier-nextflow/pre-carrier/chr22_block12/receipt.json" <<'PY'
import json,sys
r=json.load(open(sys.argv[1]));assert r['status']=='passed'
print(json.dumps(r['counts'],indent=2))
PY
```

Operator-provided full-chr22 reference totals: 7,073 missense and 6,200 HC LoF;
cohort AF >=0.001 in 173/416; POPmax >=0.001 in 223/85; region overlap in
1,430/13,273 candidates. These are not local validation results or block12 expected
counts. Their intersections and missing-AC/AN counts are needed to predict final
retained totals; this stage reports those distinctions rather than subtracting
independent totals.

## Local validation

`tests/test_pre_carrier_filter.py` uses synthetic sites-only VCFs, Parquets and
BEDs to exercise boundaries, exact keys, aliases, missing values, every track and
union, both candidate types, preserved inputs, genotype-header rejection, resource
failure cases, empty inputs, resource locking, actual Nextflow execution, durable
publication and `-resume`. Run:

```bash
uv run --with pytest --with duckdb==1.5.5 --with pysam python -m pytest -q \
  tests/test_pre_carrier_filter.py tests/test_scored_candidates.py \
  tests/test_exact_carriers.py tests/test_sites_annotation.py tests/test_sites_catalog.py
```

Local result: **54 tests passed**, including 12 pre-carrier tests and all 42
candidate/annotation/carrier/sites regression tests. NBDC config rendering,
Python compilation, and diff whitespace checks also passed. These are synthetic
validation results; the real NBDC pilot remains to be run by the operator.

The historical chr22 counts above used the earlier 0.001 cohort screen. They are not expected counts for the relaxed 0.005 screen. See [production batching and preliminary screening](CARRIER_BATCHING.md) before rerunning into a new destination. Corrected unrelated rarity is deferred.
