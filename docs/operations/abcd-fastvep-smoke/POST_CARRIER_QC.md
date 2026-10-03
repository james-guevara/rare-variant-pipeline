# Post-extraction QC of filtered carriers

`post_carrier_qc.nf` is a separate stage on completed PR #12 outputs. It reads
`carriers.tsv.gz`, `samples.tsv`, and the extraction `receipt.json` per block.
It never runs annotation, candidate filtering, or extraction and does not read
VCFs, pedigrees, or scientific reference resources. All raw outputs remain intact.

## Fixed row-level policy

A carrier association passes only when every condition is satisfied:

| Condition | Rule |
|---|---|
| Site filter | `site_FILTER == "PASS"` exactly |
| GQ | finite, nonnegative numeric value, `>= 20` |
| DP | finite, nonnegative numeric value, `>= 10` |
| AD | exactly two finite, nonnegative numeric values; finite positive sum |
| AB | `AD_alt / (AD_ref + AD_alt)` |
| Het GT | `0/1`, `1/0`, `0|1`, `1|0`: `0.25 <= AB <= 0.75` |
| Hom ALT GT | `1/1`, `1|1`: `AB >= 0.90` |
| Haploid ALT GT | `1`: `AB >= 0.90`, preserving existing QC behavior |
| Other GT | fail closed, including partial calls and unsupported ploidies |

Missing, malformed, nonfinite and negative quality/depth values fail closed.
Thresholds are fixed in this stage; no chromosome passthrough or new X/Y policy
is introduced. FT is preserved but **not filtered**; `.` is valid for this purpose.
There is no family propagation and no extra gene/MANE/region/frequency filtering.

The existing `scripts/postprocess/qc_genotype.py` establishes GT/AB behavior. A
regression invokes that script on well-formed synthetic values and checks equal
selection. Its permissive SQL parsing does not enforce the newly requested AD
cardinality/nonfinite checks, so this adapter implements those explicitly without
changing the older script or its other callers.

## Failure counts and summaries

Independent reasons are:

```text
site_filter
 gq_invalid / gq_below_min
 dp_invalid / dp_below_min
 gt_unsupported
 ad_invalid
 ab_out_of_range
```

Invalid and below-threshold flags are mutually exclusive for each scalar field.
`ab_out_of_range` applies **only when AD and GT are evaluable**. An unevaluable AB
instead contributes `ad_invalid` and/or `gt_unsupported`. Receipts include
independent counts and mutually exclusive reason combinations, with `PASS` for
passing rows. Combinations partition input rows; non-PASS combinations partition
failed rows. The same accounting is provided by candidate type, allele class and
tier, with explicit reconciliation checks.

Carrier annotations are preserved, including candidate type, tier, Gene, Feature,
SYMBOL and allele class. Untiered HC rows remain in `untiered` strata. Every tier
summary separates sequence and spanning-deletion records. No upstream-deletion
mapping or biological event deduplication is attempted.

Cross-type exact-allele/sample associations retain both annotations. Their genotype
fields must agree; conflicting raw associations fail instead of choosing one.
Within-type duplicate associations fail. Per-type counts are not additive as a
combined burden. Distinct-allele summaries count each surviving exact allele/sample
once **within its allele class**, with separate sequence and star columns.

## New durable outputs

```text
<outdir>/post-carrier-qc/<unit_id>/
  carriers.qc.tsv.gz
  qc_audit.tsv.gz
  samples.tsv
  sample_burden.tsv
  sample_gene_burden.tsv
  gene_burden.tsv
  sample_distinct_alleles.tsv
  receipt.json
```

- `carriers.qc.tsv.gz`: surviving rows, all raw columns plus `qc_AB`.
- `qc_audit.tsv.gz`: all input association keys/type/tier/gene/transcript, AB,
  pass flag and failure combination. This is protected individual-level data.
- `samples.tsv`: the complete source sample universe, validated against the
  extraction receipt; never reconstructed from carrier rows.
- `sample_burden.tsv`: long form by sample/type/class/tier, including zeros for
  every source sample and every raw observed stratum, plus untiered strata.
- `sample_gene_burden.tsv`: sparse surviving sample/gene/type/class/tier counts.
- `gene_burden.tsv`: gene/type/class/tier groups observed in raw carrier rows,
  including groups reduced to zero by QC. This is not a genome-wide gene roster.
- `sample_distinct_alleles.tsv`: all samples, two separate distinct sequence/star
  count columns. Homozygous ALT remains one carrier record with dosage two.
- `receipt.json`: stage/status/policy, all aggregate counts and combinations,
  per-type/class/tier counts, reconciliation, input/output SHA-256 and sizes,
  source extraction lineage, code identity and timing.

The input carrier/sample hashes must match their **passed extraction receipt**,
including its unit/sample/row counts. These are provenance checks, not a rerun of
resource validation. All input identities are checked again at completion. An
output directory that would overwrite input samples is rejected. Failed tasks
remove partial QC products and leave a failed receipt in work; require a successful
Nextflow exit and matching passed receipt because older publications can remain.

## NBDC: block12, resume, three blocks, then all 23

The task uses only the Python standard library, already present in the validated
LOFTEE SIF. Reuse the existing annotation-lock container identity; no new image,
packages, BEDs or dbNSFP files are needed. No NBDC jobs were launched by Codex.

Update without discarding local modifications:

```bash
cd "$HOME/rare-variant-pipeline"
git fetch origin
git switch --detach origin/feat/post-carrier-qc
git rev-parse HEAD
BASE=/home/ood-guevara-james/abcd-fastvep-smoke
```

Build a manifest of already-completed extraction outputs. Listing all blocks here
does not execute them: selection is explicit at launch.

```bash
python3 - "$BASE" > "$PWD/abcd-chr22-post-carrier-qc.tsv" <<'PY'
from pathlib import Path
import sys
root=Path(sys.argv[1])/'filtered-carrier-nextflow'/'filtered-carriers'
print('unit_id\tchromosome\tcarriers\tsamples\tsource_receipt')
for i in range(23):
    unit=f'chr22_block{i}';p=root/unit
    paths=[p/'carriers.tsv.gz',p/'samples.tsv',p/'receipt.json']
    if not all(x.is_file() for x in paths):
        raise SystemExit(f'Missing extraction product for {unit}')
    print('\t'.join([unit,'chr22',*map(str,paths)]))
PY

nextflow -C post_carrier_qc.config run post_carrier_qc.nf -profile nbdc \
  --qc_manifest "$PWD/abcd-chr22-post-carrier-qc.tsv" \
  --resource_lock "$HOME/rare-variant-pipeline/abcd-annotation-resources.json" \
  --select_units chr22_block12 \
  --outdir "$BASE/post-carrier-qc-nextflow" \
  -work-dir "$BASE/post-carrier-qc-nextflow-work" \
  -with-trace "$BASE/post-carrier-qc-block12.trace.tsv" -resume
```

Default NBDC resources remain configurable (`--qc_cpus`, `--qc_memory`, `--qc_time`,
`--qc_queue`, `--qc_qos`, `--qc_queue_size`): 2 CPUs, 4 GB, 4h, medium partition/QoS,
queueSize 4. These execution settings do not expose scientific threshold overrides.

Check the block12 receipt (aggregate output only):

```bash
python3 - "$BASE/post-carrier-qc-nextflow/post-carrier-qc/chr22_block12/receipt.json" <<'PY'
import json,sys
r=json.load(open(sys.argv[1]))
assert r['status']=='passed' and r['stage']=='post_extraction_qc'
assert r['samples_in_source']==8877
assert r['reconciliation']['passed']
assert r['input_rows']==r['pass_rows']+r['fail_rows']
assert sum(r['failure_combinations'].values())==r['input_rows']
assert r['failure_combinations'].get('PASS',0)==r['pass_rows']
assert sum(n for k,n in r['failure_combinations'].items() if k!='PASS')==r['fail_rows']
for field in ['input_rows','pass_rows','fail_rows']:
    assert sum(t[field] for t in r['by_candidate_type'].values())==r[field]
    for t in r['by_candidate_type'].values():
        assert sum(c[field] for c in t['by_allele_class'].values())==t[field]
        for c in t['by_allele_class'].values():
            assert sum(x[field] for x in c['by_tier'].values())==c[field]
print(json.dumps({k:r[k] for k in ['input_rows','pass_rows','fail_rows',
    'independent_failures','failure_combinations','by_candidate_type',
    'post_qc_distinct_alleles_by_class']},indent=2))
PY
```

Then repeat the **same command**, changing only the trace filename to
`$BASE/post-carrier-qc-block12-resume.trace.tsv`; verify `POST_CARRIER_QC
(chr22_block12)` is `CACHED` in the trace.

Next use `--select_units chr22_block0,chr22_block12,chr22_block19` and a new trace
filename with otherwise identical command/work/output paths. Block12 should remain
cached. Check each passed receipt and count reconciliation before proceeding to
`--select_units all` with another trace filename. No upstream entrypoint is invoked.

## Operator-observed whole-chr22 reference counts

These are aggregate **NBDC observations supplied by the operator**, not synthetic
test expectations and not predicted block12 counts:

| Type / allele class / tier | Before | After |
|---|---:|---:|
| HC / sequence / lof_t1 | 1116 | 73 |
| HC / sequence / lof_t2 | 3437 | 474 |
| HC / sequence / untiered | 5930 | 1665 |
| HC / star / lof_t1 | 160 | 34 |
| HC / star / lof_t2 | 1047 | 124 |
| HC / star / untiered | 1047 | 239 |
| Missense / sequence / miss_t1 | 33 | 16 |
| Missense / sequence / miss_t2 | 209 | 110 |
| Missense / sequence / miss_t3 | 1032 | 731 |
| Missense / sequence / miss_t4 | 7812 | 6886 |

Compare whole-chr22 summed stratum counts with these after all 23 receipts pass.
Do not sum candidate types or star/sequence strata into a biological burden.
**Deferred audit note:** the operator observed that most LoF-indel losses were low
AB despite high GQ. Thresholds remain fixed; this implementation does not interpret
those losses or reopen the scientific policy. Differences from the reference
counts should be investigated using failure combinations/provenance, not fixed by
changing thresholds.

## Validation and limits

Focused tests exercise all thresholds, malformed/missing/nonfinite/negative
numbers, AD cardinality and zero/overflow sum, phased/haploid/unsupported GT,
FT independence, star/sequence and tier separation, untiered HC, no propagation,
zero-carrier samples, empty inputs, type overlaps, provenance failures, unchanged
inputs, original-QC parity on valid values, actual PR #12 synthetic extraction
products, and real local Nextflow subset/resume and durable publication.

```bash
uv run --with pytest --with duckdb==1.5.5 --with pysam python -m pytest -q \
  tests/test_post_carrier_qc.py
```

Python runtime needs only the standard library; DuckDB/pysam above are test-only
for comparisons and creating upstream fixtures. Real NBDC data, Slurm execution,
and Apptainer execution of this new stage have **not** been tested by Codex. The
existing pinned image's Python/dependencies were verified during PR #12. The
NBDC config is rendered locally; the operator executes the staged rollout above.

Local result: **62 tests passed**. Python compilation, NBDC profile rendering,
and whitespace/link checks also passed. No existing extraction/QC scripts or
scientific filters were modified by this stage.
