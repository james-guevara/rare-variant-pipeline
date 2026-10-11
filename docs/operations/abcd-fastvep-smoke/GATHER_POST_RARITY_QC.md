# Chromosome gathering after final rarity and sex/PAR-aware QC

Current validation status: the operator reports completed NBDC validation through
976-block/24-chromosome gathering. See the [finalization record](../../validation/abcd-post-rarity-finalization.md).
Historical code pins and stage-specific instructions below remain unchanged.

Branch: `feat/gather-post-rarity-qc`, stacked on PR #21.
Code pin: `aed4008eab38073110896f5d1e170ae698c9e9ca`.
Entrypoint/config: `gather_post_rarity_qc.nf` / `gather_post_rarity_qc.config`.

The operator reports post-rarity QC validation on NBDC: 976/976 receipts
reconciled; 941,371 input carrier annotations → 454,232 passing, comprising
107,704 sequence HC LoF, 329,667 missense and 16,861 spanning-deletion
associations. These are supplied validation results, not measurements by Codex.
Gathering consumes these saved passing products without rerunning QC, frequency
counting, extraction or other upstream stages.

## Contract and outputs

Manifest columns are `unit_id`, `chromosome`, `carriers`, `samples`,
`source_receipt`, tab-separated. Each row names saved `carriers.qc.tsv.gz`,
`samples.tsv` and the matching passed `post_rarity_qc` receipt. The builder below
inventories that root once; it does not inspect original VCFs. Explicit unique
unit IDs and chromosome grouping come from receipts, not filename guesses.

The new entrypoint uses the existing gather implementation through an adapter:

- Require the exact PR #21 stage/policy, reconciled source counts (including
  nested type/class/tier totals), matching output hashes and source sample order.
- Require PSAM content identities to agree with each QC receipt's input identity
  and across all selected blocks/chromosomes. The saved bytes/SHA-256 identity
  is used; no new PSAM file or genotype source is needed.
- Validate saved GT/ploidy/region/sex/effective-dosage consistency. Include all
  ploidy metadata in cross-type genotype comparison. This checks the saved QC
  contract; it does not recalculate AF or apply new quality thresholds.
- Reject cross-block exact allele overlaps, including overlaps with different
  carrier samples. Do not silently deduplicate. Reject duplicate within-type
  allele/sample associations; preserve distinct candidate-type annotations.
- Sum `qc_effective_alt_dosage` for allele totals. Preserve raw GT/alt_dosage and
  every source annotation column in gathered carriers; add `source_unit_id`.
- Recompute sample sets and distinct allele/sample unions from surviving rows.
  Never add block-level unique-sample counts. Keep sequence/stars separate.

One task per selected chromosome writes durable copies to:
`<outdir>/post-rarity-qc-gather/<chromosome>/<gather_id>/`.

| Product | Meaning |
|---|---|
| `carriers.qc.tsv.gz` | All selected passing annotations, raw fields retained, source unit added |
| `samples.tsv` | All source samples in source order |
| `sample_burden.tsv` | Dense sample × type/class/tier association counts and effective allele totals, including zeros |
| `sample_gene_burden.tsv` | Sparse sample/gene × type/class/tier association counts and effective allele totals |
| `gene_burden.tsv` | Gene × type/class/tier associations, recomputed unique carrier samples and effective allele totals |
| `sample_distinct_alleles.tsv` | Dense per-sample exact-allele unions separately for sequence and stars |
| `sample_gene_distinct_alleles.tsv` | Sparse per-sample/gene exact-allele unions separately by allele class |
| `receipt.json` | Policy, PSAM/input/output/code identities, source lineage, strata, class totals, duplicate checks and reconciliations |

Sparse missing cells mean zero; dense sample summaries and `samples.tsv` retain
every source sample, including all 8,877 in this dataset. Tiers and untiered HC
remain separate. Association-based dosage can count the same allele/sample in
multiple candidate types. Receipt `distinct_allele_sample_observed_alt_alleles`
is the effective-dosage union across types within each allele class;
`annotation_observed_alt_alleles` is the association-based total. Neither is a
sum of sequence and spanning-deletion biological burden.

This remains chromosome-wide gathering. Do not sum chromosome unique-sample
counts to claim a genome-wide unique-sample count. No genome-wide gather or
change to the legacy `gather_post_qc.nf` entrypoint is required for this stage.
The legacy entrypoint retains its old policy and raw dosage semantics.

## NBDC setup and block12 pilot

The NBDC profile reuses the existing LOFTEE SIF and annotation lock; no new
container, resource download or PSAM staging is required. Defaults: 1 CPU/task,
4 GB, 2 hours, 16-task cap, `small,medium` partitions, no forced QoS, trace
overwrite enabled. Scheduler account limits govern actual concurrency. Adjust
`--gather_memory`/`--gather_time` only if observed gathering resource use requires
it. This aggregates in memory per chromosome; NBDC peak memory remains untested.

Run setup once. Manifest creation refuses to overwrite an existing manifest.

```bash
set -euo pipefail
cd "$HOME/rare-variant-pipeline"
git fetch origin feat/gather-post-rarity-qc
git switch --detach aed4008eab38073110896f5d1e170ae698c9e9ca
BASE=/home/ood-guevara-james/abcd-fastvep-smoke
REPO=$PWD
RUN="$BASE/post-rarity-qc-gather-nextflow"
mkdir -p "$RUN"
python3 "$REPO/scripts/make_post_rarity_gather_manifest.py" \
  --root "$BASE/post-rarity-qc-nextflow/output/post-rarity-qc" \
  --output "$RUN/manifest.tsv"
# Expected inventory: 976 units, 951 autosomal and 25 X/Y.
cd "$RUN"
run_gather() {
  nextflow -C "$REPO/gather_post_rarity_qc.config" run "$REPO/gather_post_rarity_qc.nf" \
    -profile nbdc --gather_manifest "$RUN/manifest.tsv" \
    --resource_lock "$REPO/abcd-annotation-resources.json" \
    --select_units "$1" --gather_id "$2" --outdir "$RUN/output" \
    -work-dir "$RUN/work" -with-trace "$RUN/trace.tsv" -resume
}
run_gather chr22_block12 block12
```

Require successful Nextflow completion, then use this aggregate-only checker.
It verifies selected unit coverage, source-count reconciliation, 8,877 samples,
PSAM consistency, no cross-block duplicates, and published hashes. Keep
individual-level files local.

```bash
check_gather() {
  python3 - "$RUN/manifest.tsv" "$RUN/output/post-rarity-qc-gather" "$1" "$2" <<'PY'
from collections import Counter, defaultdict
import csv, hashlib, json, pathlib, sys
manifest, root, selection, gid = sys.argv[1:]
all_rows = list(csv.DictReader(open(manifest), delimiter='\t'))
wanted = set(selection.split(',')) if selection != 'all' else {r['unit_id'] for r in all_rows}
assert wanted <= {r['unit_id'] for r in all_rows}
groups = defaultdict(list)
for row in all_rows:
    if row['unit_id'] in wanted: groups[row['chromosome']].append(row)
totals = Counter(); strata = Counter(); psams = set()
for chrom, units in groups.items():
    directory = pathlib.Path(root) / chrom / gid
    r = json.loads((directory / 'receipt.json').read_text())
    assert r['stage'] == 'post_rarity_qc_gather' and r['status'] == 'passed'
    assert r['chromosome'] == chrom and r['gather_id'] == gid
    assert r['reconciliation']['passed'] and r['samples_in_source'] == 8877
    assert r['policy']['version'] == 'post-rarity-grch38-x-only-par-qc-v1'
    assert r['summary_dosage_field'] == 'qc_effective_alt_dosage'
    assert set(r['selected_units']) == {u['unit_id'] for u in units}
    assert not any(r['duplicate_checks'].values())
    psams.add((r['psam_identity']['bytes'], r['psam_identity']['sha256']))
    expected = Counter()
    for unit in units:
        source = json.loads(pathlib.Path(unit['source_receipt']).read_text())
        assert source['psam_identity'] == r['psam_identity']
        expected.update({k: source[k] for k in ['input_rows','pass_rows','fail_rows']})
    assert dict(expected) == r['source_counts']
    assert r['carrier_annotation_records'] == expected['pass_rows']
    assert sum(s['carrier_annotation_records'] for s in r['by_stratum']) == expected['pass_rows']
    assert sum(c['carrier_annotation_records'] for c in r['by_allele_class'].values()) == expected['pass_rows']
    assert sum(s['observed_alt_alleles'] for s in r['by_stratum']) == sum(c['annotation_observed_alt_alleles'] for c in r['by_allele_class'].values())
    if selection == 'all': assert r['full_manifest_selected']
    for name, identity in r['outputs'].items():
        p = directory / name; h = hashlib.sha256()
        with p.open('rb') as f:
            for b in iter(lambda: f.read(8 * 1024 * 1024), b''): h.update(b)
        assert p.stat().st_size == identity['bytes'] and h.hexdigest() == identity['sha256']
    totals.update(expected)
    for s in r['by_stratum']:
        strata[s['candidate_type'], s['allele_class']] += s['carrier_annotation_records']
assert len(psams) == 1
if selection == 'all':
    assert len(wanted) == 976 and len(groups) == 24
    assert totals['input_rows'] == 941371 and totals['pass_rows'] == 454232
    assert strata['lof_hc','sequence'] == 107704
    assert strata['missense','sequence'] == 329667
    assert sum(n for (t,c),n in strata.items() if c == 'spanning_deletion') == 16861
print(json.dumps(dict(chromosomes=len(groups), units=len(wanted), **totals)))
PY
}
check_gather chr22_block12 block12
```

These full-rollout expectations are the supplied NBDC observations, not synthetic
test expectations. If they disagree, inspect receipts before proceeding; do not
change the thresholds or drop records to force agreement.

## Resume, mixed pilots and full rollout

```bash
run_gather chr22_block12 block12
check_gather chr22_block12 block12
# Require CACHED in the overwritten trace.tsv.
PILOTS=$(python3 - "$RUN/manifest.tsv" <<'PY'
import csv, sys
rows=list(csv.DictReader(open(sys.argv[1]),delimiter='\t'))
units=['chr22_block0','chr22_block12','chr22_block19']
units += [next(r['unit_id'] for r in rows if r['chromosome']==c) for c in ['chrX','chrY']]
print(','.join(units))
PY
)
run_gather "$PILOTS" mixed-pilots
check_gather "$PILOTS" mixed-pilots
run_gather "$PILOTS" mixed-pilots
check_gather "$PILOTS" mixed-pilots
# Require three CACHED chromosome tasks in trace.tsv.
```

After pilots reconcile, inspect effective versus raw dosage preservation and the
separate star/untiered strata locally, then expand:

```bash
run_gather all all
check_gather all all
# Expect 24 chromosome receipts covering all 976 blocks.
```

An identical selection, gather ID, code and inputs can resume cached tasks.
Changing selection from one block to three or all blocks necessarily recomputes
that chromosome's gather. Distinct gather IDs keep pilot and complete outputs
separate. No upstream work runs. Preserve the launch/work directories for resume.

## Validation scope

The combined local suite passed **175 tests**.

Synthetic tests exercise effective haploid dosage and preserved raw calls,
PSAM mismatch within/across chromosomes, saved ploidy consistency, cross-type
conflicts/unions, cross-block duplicate alleles even with different samples,
input hashes, fixed policy/counts, tiers/untiered HC/stars, all 8,877 samples with
empty inputs, recomputed unique-sample counts and immutable inputs. Local
Nextflow tests run block selection, mixed chromosomes, full selection and cached
resume using an overwritten trace, then verify durable products after removing
test work files. Existing gathering, QC and final-rarity regression suites are
also run. No NBDC job, real SIF execution or real-data gathering was performed
by Codex.
