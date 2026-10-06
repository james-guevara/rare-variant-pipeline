# Chromosome gathering after carrier QC

This is a separate `gather_post_qc.nf` entrypoint stacked on PR #13 at
`1e1731960ca750f396d02bac65899050fafb3273`. The old `GATHER_CARRIERS` module
only concatenates older chunk TSVs and cannot supply these summaries or checks.
No existing workflow is modified. No VCF, annotation, filtering, extraction or
QC task runs here. The stage reads only published `carriers.qc.tsv.gz`,
`samples.tsv`, and `receipt.json` for selected units.

The operator reports that PR #13 passed block12/resume, blocks 0/12/19, and all
23 blocks on NBDC, with all receipts reconciled: 21,823 input associations →
10,352 passing (7,743 missense, 2,212 sequence HC LoF, 397 star associations).
These are operator-provided observations, not synthetic test expectations or
Codex-observed NBDC results. PR #13 was still open when this stage was started.

## Counting and integrity rules

- Verify each product against its passed post-QC receipt, unit, chromosome,
  fixed QC policy, sample count, passed-row count and type/class/tier counts.
- Require the identical full sample list/order across blocks. Preserve every
  source sample in `samples.tsv`, `sample_burden.tsv` and
  `sample_distinct_alleles.tsv`, even if all its counts are zero.
- Normalize only the `chr` prefix for exact CHROM/POS/REF/ALT identity.
- **Fail** if any passing exact allele appears in more than one selected block,
  even when those blocks list different samples or annotations. Also report
  duplicate cross-block allele/sample/type associations. No arbitrary winner or
  silent deduplication is used. Only post-QC carrier alleles are inspected; this
  does not claim to detect overlaps among sites with no passing carriers.
- Reject within-unit duplicate allele/sample/type associations. Preserve both
  candidate types for a within-unit overlap; their genotype fields must agree.
- Recompute unique samples using sets. Never add block-level unique-sample counts.
- Annotation-association burdens retain type, tier, gene, and allele class.
  Untiered HC is retained. Homozygous ALT is one association with dosage two.
- Distinct allele/sample burdens union across candidate types **within allele
  class**. Per-gene distinct burdens union across types within the same gene.
  Distinct counts across different genes are not additive as an overall burden.
- Star counts remain separate in every summary. They are representation records,
  not inferred independent deletion events. No upstream-deletion mapping occurs.
- Per-sample/per-gene tables are sparse: absent cells are zero. Join to `samples.tsv`
  to retain zero-carrier samples in downstream analyses; no arbitrary gene universe
  or huge sample × gene zero cross-product is generated.
  Zero means no retained post-QC candidate carrier in this selected scope; it
  is not evidence of complete callable coverage or absence of every variant.

## Products

```text
<outdir>/post-qc-gather/<chromosome>/<gather_id>/
  carriers.qc.tsv.gz              # preserves all columns plus source_unit_id
  samples.tsv                    # complete sample universe
  sample_burden.tsv              # sample/type/class/tier association counts + dosage
  sample_gene_burden.tsv          # sparse sample/gene/type/class/tier associations
  gene_burden.tsv                 # gene/type/class/tier associations + unique samples
  sample_distinct_alleles.tsv     # per-sample sequence/star distinct counts
  sample_gene_distinct_alleles.tsv # sparse per-sample/gene/class distinct counts
  receipt.json                   # aggregate counts, reconciliation, hashes, lineage
```

The receipt reports counts and unique samples by type/class/tier and by allele
class, source QC input/pass/fail sums, duplicate checks, selected units, full
manifest units, container lock identity, input/output hashes and script identities.
`full_manifest_selected` means every chromosome unit in the supplied manifest was
selected; it cannot prove that the operator supplied every real block. For ABCD
chr22, explicitly check the expected block0–block22 set as below.

Use a distinct `gather_id` for each pilot scope. It prevents a one-block result
from overwriting a chromosome-wide publication. Resume requires the same manifest,
selection, gather ID, work directory, inputs and code. Expanding selection reruns
only this small gather task; it does not run upstream stages. Failed tasks remove
partial products and leave a failed receipt in work. Always require a successful
Nextflow exit and matching receipt; an older publication can remain after failure.

## Pinned NBDC rollout

Implementation commit: `282a0862cb479deeb284e54b0469f5e1c7e45450`.
Keep local modifications; do not force checkout or reset them.

```bash
cd "$HOME/rare-variant-pipeline"
git fetch origin
git switch --detach 282a0862cb479deeb284e54b0469f5e1c7e45450
BASE=/home/ood-guevara-james/abcd-fastvep-smoke
python3 - "$BASE" > "$PWD/abcd-chr22-gather-post-qc.tsv" <<'PY'
from pathlib import Path
import sys
root=Path(sys.argv[1])/'post-carrier-qc-nextflow'/'post-carrier-qc'
print('unit_id\tchromosome\tcarriers\tsamples\tsource_receipt')
for i in range(23):
    unit=f'chr22_block{i}';p=root/unit
    paths=[p/'carriers.qc.tsv.gz',p/'samples.tsv',p/'receipt.json']
    if not all(x.is_file() for x in paths):
        raise SystemExit(f'Missing post-QC product for {unit}')
    print('\t'.join([unit,'chr22',*map(str,paths)]))
PY

nextflow -C gather_post_qc.config run gather_post_qc.nf -profile nbdc \
  --gather_manifest "$PWD/abcd-chr22-gather-post-qc.tsv" \
  --resource_lock "$PWD/abcd-annotation-resources.json" \
  --select_units chr22_block12 --gather_id block12 \
  --outdir "$BASE/post-qc-gather-nextflow" \
  -work-dir "$BASE/post-qc-gather-nextflow-work" \
  -with-trace "$BASE/gather-block12.trace.tsv" -resume
```

The existing validated LOFTEE SIF is reused for standard-library Python only;
LOFTEE itself is not invoked. NBDC defaults are 2 CPUs, 4 GB, 2 hours, medium
partition/QoS, queue size 4. Override `gather_cpus`, `gather_memory`, `gather_time`,
`gather_queue`, `gather_qos`, or `gather_queue_size` for execution needs.
No scientific thresholds are configurable in this gather stage.

Check the receipt using the command below with `block12`. Then repeat the exact
Nextflow command, changing only the trace filename to `gather-block12-resume.trace.tsv`;
verify `GATHER_POST_QC (chr22:block12)` is `CACHED`.

Next change selection to `chr22_block0,chr22_block12,chr22_block19`, gather ID to
`three-blocks`, and trace filename to `gather-three-blocks.trace.tsv`. Check that
receipt before using `--select_units all --gather_id chr22-all` and trace
`gather-chr22-all.trace.tsv`. Keep work/output paths fixed. Each new selection is
one new gather task, not a rerun of any block-level stage.

Aggregate-only checker (change the final argument for each stage):

```bash
python3 - "$BASE/post-qc-gather-nextflow/post-qc-gather/chr22" block12 <<'PY'
import json,sys
from pathlib import Path
root=Path(sys.argv[1]);scope=sys.argv[2]
r=json.loads((root/scope/'receipt.json').read_text())
expected={'block12':{'chr22_block12'},
          'three-blocks':{'chr22_block0','chr22_block12','chr22_block19'},
          'chr22-all':{f'chr22_block{i}' for i in range(23)}}[scope]
assert r['status']=='passed' and r['stage']=='post_qc_gather'
assert set(r['selected_units'])==expected and r['samples_in_source']==8877
assert r['reconciliation']['passed'] and not any(r['duplicate_checks'].values())
assert r['source_counts']['pass_rows']==r['carrier_annotation_records']
assert r['source_counts']['input_rows']==r['source_counts']['pass_rows']+r['source_counts']['fail_rows']
assert sum(x['carrier_annotation_records'] for x in r['by_stratum'])==r['carrier_annotation_records']
if scope=='chr22-all':
    assert r['full_manifest_selected']
    assert r['source_counts']['input_rows']==21823
    assert r['carrier_annotation_records']==10352
    assert sum(x['carrier_annotation_records'] for x in r['by_stratum'] if x['candidate_type']=='missense')==7743
    assert sum(x['carrier_annotation_records'] for x in r['by_stratum'] if x['candidate_type']=='lof_hc' and x['allele_class']=='sequence')==2212
    assert r['by_allele_class']['spanning_deletion']['carrier_annotation_records']==397
print(json.dumps({k:r[k] for k in ['status','samples_in_source','carrier_annotation_records',
    'source_counts','by_stratum','by_allele_class','duplicate_checks','reconciliation']},indent=2))
PY
```

These expected counts refer to annotation associations. Unique allele/sample and
unique-sample totals are recomputed outputs, not assumed equal to association
counts. A discrepancy or duplicate failure should be investigated from local
protected products and aggregate receipts, not fixed by changing QC thresholds.
Do not print or share individual carriers or sample-level burden rows.

## Local validation and limits

```bash
uv run --with pytest --with duckdb==1.5.5 --with pysam python -m pytest -q \
  tests/test_gather_post_qc.py tests/test_post_carrier_qc.py
nextflow -C gather_post_qc.config config -profile nbdc -flat
```

Tests cover receipt/product identity, sample-set/order consistency, exact-allele
cross-block duplicates including chr aliases and different samples, within-unit
duplicates, cross-type gene union, untiered/star separation, empty inputs, all
8,877 zero-carrier samples, input preservation, distinct vs association counts,
recomputed unique samples, and actual local Nextflow subset/resume/multi-block
execution plus durable publications after work deletion.

Local result: **77 tests passed**, including actual Nextflow resume and multi-block
execution. NBDC profile rendering, Python compilation and diff checks passed.

Real NBDC/Slurm/Apptainer execution and real-data gather counts remain untested
by Codex. No NBDC jobs were launched. Deferred scientific auditing, QC-policy
changes, upstream deletion deduplication, and sex-chromosome policy are out of scope.
