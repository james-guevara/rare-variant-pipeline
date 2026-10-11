# Post-rarity QC with explicit GRCh38 X-only PAR handling

Code pin: `a107a6e920ac18ddc7dee56be72a8d1f9b9b30d8`, branch
`feat/post-rarity-qc`, stacked on `feat/final-rarity` (PR #20).
Entrypoint: `post_rarity_qc.nf`; config: `post_rarity_qc.config`.

The operator reports NBDC final rarity passed 976/976 reconciled receipts:
493,527 → 481,905 alleles and 1,356,838 → 941,371 carrier annotations.
All 11,622 failed alleles had unrelated AF ≥0.001. These are operator-supplied
results, not NBDC measurements by Codex. This QC implementation has local
synthetic validation; real NBDC QC validation remains to be performed.

## Inputs and policy

Each independent block reads only saved `carriers.rare.tsv.gz`, `samples.tsv`,
its final-rarity `receipt.json`, and the existing PSAM. It verifies the passed
receipt, exact final-rarity policy, chromosome, count reconciliation, carrier
and sample SHA-256 hashes, sample universe and observed carrier count. Input
identities are checked again after processing. PSAM validation reuses the
existing generic loader; all source samples require metadata. Representative
and unrelated flags do not limit carrier QC. No frequency recalculation,
annotation, genotype VCF access, or upstream stage is invoked.

Quality rules reuse the existing evaluator: site FILTER exactly PASS, GQ ≥20,
DP ≥10, biallelic AD with two finite nonnegative values and positive sum;
AB = AD_alt / (AD_ref + AD_alt). Heterozygous AB is inclusive 0.25–0.75;
homozygous/haploid ALT requires AB ≥0.90. Missing/malformed/nonfinite quality
values fail closed. FT is preserved, never required to be PASS. No family
propagation is performed.

SEX uses PSAM encoding 1=male, 2=female, other nonempty values=unknown.

| Region | Expected ploidy and handling |
|---|---|
| Autosome | Diploid, regardless of sex; unexpected ploidy excluded |
| X PAR | Diploid for both known sexes |
| X non-PAR | Male haploid, female diploid |
| Y non-PAR | Male haploid, females excluded |
| X/Y unknown sex | Excluded and audited |
| Y PAR | Stage fails clearly under X-only PAR policy |

GRCh38 1-based inclusive PAR intervals are X: 10001–2781479 and
155701383–156030895; Y: 10001–2781479 and 56887903–57217415.
This reuses the extraction-frequency policy implementation. Partial calls and
unsupported ploidy fail. Expected-haploid regions accept literal `1` and
`1/1` or `1|1`; heterozygotes fail even with otherwise valid AB.
The legacy QC entrypoint retains its historical chromosome-blind behavior;
this new entrypoint explicitly excludes haploid autosomal calls.

## Products and counts

Durable products: `<outdir>/post-rarity-qc/<unit_id>/`:

- `carriers.qc.tsv.gz`: passing annotations, preserving raw fields including GT,
  alt_dosage, type, tier, untiered HC, Gene/Feature and allele_class. Adds qc_AB,
  qc_expected_ploidy, qc_frequency_region, qc_sex, qc_effective_alt_dosage.
- `qc_audit.tsv.gz`: all input annotations, genotype/dosage, QC context and
  failure combinations. Effective dosage is missing on failed rows.
- `samples.tsv`, dense `sample_burden.tsv` and `sample_distinct_alleles.tsv`:
  all source samples, including zeros. Sparse `sample_gene_burden.tsv` and
  `gene_burden.tsv` describe surviving associations.
- `receipt.json`: status/policy, input/output hashes, code identities, reconciled
  input/pass/fail counts by type/class/tier, independent failure reasons and
  mutually exclusive combinations, distinct allele/sample counts by class.

`observed_alt_alleles` summaries use effective dosage: one for eligible haploid
calls, raw allele count for diploid calls. Raw `alt_dosage` is never rewritten.
A homozygous carrier still contributes one carrier-variant record per annotation.
Cross-type associations remain separate; distinct allele/sample counts union
those associations. Stars remain separate from sequence in every burden count.
AB-range failure is not assigned to sex/ploidy-ineligible calls, where that
interpretation is inapplicable; independent GQ/DP/AD/site failures still count.
Keep audit and sample-level products local; share aggregate results only.

## Pinned NBDC runbook

Reuse the validated LOFTEE SIF and existing annotation resource lock. No new
container or package installation is required (existing Python/pysam suffice).
Default requests: 1 CPU/task, 4 GB, one hour, maximum 16 tasks, partitions
`small,medium`, no forced QoS, trace overwrite enabled. Memory/time and queue
are configurable with `--post_rarity_qc_memory`, `--post_rarity_qc_time`,
`--post_rarity_qc_queue`; account limits still govern actual concurrency.

Run setup once; manifest creation refuses to overwrite an existing manifest.
Subsequent resumes reuse the same launch directory, work directory and manifest.

```bash
set -euo pipefail
cd "$HOME/rare-variant-pipeline"
git fetch origin feat/post-rarity-qc
git switch --detach a107a6e920ac18ddc7dee56be72a8d1f9b9b30d8
BASE=/home/ood-guevara-james/abcd-fastvep-smoke
REPO=$PWD
RUN="$BASE/post-rarity-qc-nextflow"
PSAM="$BASE/metadata/release-7.0/abcd-frequency-samples.psam"
mkdir -p "$RUN"
python3 "$REPO/scripts/make_post_rarity_qc_manifest.py" \
  --root "$BASE/final-rarity-nextflow/output/final-rarity" \
  --output "$RUN/manifest.tsv"
# Expected inventory: 976 units, 951 autosomal, 25 X/Y.
cd "$RUN"
run_qc() {
  nextflow -C "$REPO/post_rarity_qc.config" run "$REPO/post_rarity_qc.nf" \
    -profile nbdc \
    --post_rarity_qc_manifest "$RUN/manifest.tsv" \
    --psam "$PSAM" --sex_chromosome_policy grch38_x_only_par \
    --resource_lock "$REPO/abcd-annotation-resources.json" \
    --select_units "$1" --outdir "$RUN/output" \
    -work-dir "$RUN/work" -with-trace "$RUN/trace.tsv" -resume
}
run_qc chr22_block12
```

After successful Nextflow completion, reconcile each published unit with the
following aggregate-only checker (also use it after each expansion):

```bash
python3 - "$RUN/manifest.tsv" "$RUN/output/post-rarity-qc" chr22_block12 <<'PY'
import csv, hashlib, json, pathlib, sys
manifest, root, selected = sys.argv[1:]
units = [r['unit_id'] for r in csv.DictReader(open(manifest), delimiter='\t')]
if selected != 'all':
    wanted = selected.split(',')
    assert set(wanted) <= set(units)
    units = wanted
counts = dict(input=0, passed=0, failed=0)
for unit in units:
    directory = pathlib.Path(root) / unit
    r = json.loads((directory / 'receipt.json').read_text())
    assert r['unit_id'] == unit and r['status'] == 'passed'
    assert r['stage'] == 'post_rarity_qc' and r['reconciliation']['passed']
    assert r['policy']['version'] == 'post-rarity-grch38-x-only-par-qc-v1'
    assert r['samples_in_source'] == 8877
    assert r['input_rows'] == r['pass_rows'] + r['fail_rows']
    assert sum(r['failure_combinations'].values()) == r['input_rows']
    source = next(x for x in csv.DictReader(open(manifest), delimiter='\t') if x['unit_id'] == unit)
    upstream = json.loads(pathlib.Path(source['source_receipt']).read_text())
    assert r['input_rows'] == upstream['carrier_annotations']['pass_count']
    for name, identity in r['outputs'].items():
        p = directory / name
        h = hashlib.sha256()
        with p.open('rb') as f:
            for b in iter(lambda: f.read(8 * 1024 * 1024), b''): h.update(b)
        assert p.stat().st_size == identity['bytes'] and h.hexdigest() == identity['sha256']
    for t in r['by_candidate_type'].values():
        assert t['input_rows'] == t['pass_rows'] + t['fail_rows']
        for c in t['by_allele_class'].values():
            assert c['input_rows'] == c['pass_rows'] + c['fail_rows']
            assert sum(v['pass_rows'] for v in c['by_tier'].values()) == c['pass_rows']
    for k, field in [('input','input_rows'),('passed','pass_rows'),('failed','fail_rows')]: counts[k] += r[field]
print(json.dumps(dict(units=len(units), **counts)))
PY
```

Repeat `run_qc chr22_block12` and require CACHED in overwritten trace.tsv.
Then select the three autosomal pilots and one existing unit each on X and Y:

```bash
PILOTS=$(python3 - "$RUN/manifest.tsv" <<'PY'
import csv, sys
rows=list(csv.DictReader(open(sys.argv[1]),delimiter='\t'))
units=['chr22_block0','chr22_block12','chr22_block19']
units += [next(r['unit_id'] for r in rows if r['chromosome']==c) for c in ['chrX','chrY']]
print(','.join(units))
PY
)
run_qc "$PILOTS"
```

Rerun the checker replacing its last argument with `"$PILOTS"`. Inspect aggregate
unknown-sex/female-Y/ploidy/partial/heterozygote failure counts and raw-versus-
effective dosage by class. Repeat pilots to verify caching, then run
`run_qc all` and check with last argument `all`. Require 976 passed receipts.
No real expected QC pass count is asserted before this validation. Keep outputs
and Nextflow cache durable; do not remove the work directory before resuming.

## Downstream gathering compatibility

Use the separate [post-rarity gathering adapter](GATHER_POST_RARITY_QC.md)
for these products. Do not feed them into the legacy `gather_post_qc.nf`.
Its reader requires the old `post_extraction_qc` stage and chromosome-blind policy;
it intentionally rejects these receipts (tested). It also sums raw dosage.
The per-block table keys, carrier filename and receipt strata remain structurally
compatible, but the new gather adapter validates the new stage/policy,
checks consistent PSAM identities, uses `qc_effective_alt_dosage` for allele totals,
and includes ploidy metadata in cross-type genotype consistency checks. Preserve
cross-block duplicate detection, all source samples, distinct allele/sample
unions and separate stars; do not mix old/new QC policies. The QC entrypoint never invokes gathering; its separate runbook provides the
new chromosome gathering commands.

## Local validation

Focused synthetic tests cover PAR boundaries, both sexes/unknown sex, haploid
encodings, partial/unexpected ploidy, exact quality boundaries, invalid fields,
FT behavior, immutable inputs, receipt/hash failures, tiers/stars/untiered HC,
cross-type distinct counts, zeros, empty/all-fail units, actual local Nextflow
execution and cached resume, and legacy gather rejection. The combined suite passed 157 tests, including existing QC, final
rarity and gathering regression tests. Real SIF/Slurm execution,
NBDC output counts and full-cohort resource use remain untested by Codex.
