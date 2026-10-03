# Carrier extraction from filtered candidate Parquets

Use `filtered_carriers.nf` after the completed pre-carrier stage. It consumes
`missense.filtered.parquet` and `lof_hc.filtered.parquet`, plus the original block
VCF and its explicit TBI/CSI index. It does not rerun annotation, scoring, filtering,
or reference-resource validation. The broader scientific audit remains deferred.

The operator reports all 23 chr22 pre-carrier blocks passed and resume cached:
6,093 retained missense, 4,504 sequence HC-LoF, and 611 star records. These are input
candidate counts supplied by the operator, not results measured by Codex.

## Implementation and compatibility

The new manifest, process and Parquet adapter are separate from `carriers.nf` and
its HC TSV interface. The indexed lookup and genotype decoding were factored into
shared helpers in `extract_exact_carriers.py`; both entrypoints use those helpers.
The existing HC-only CLI, selection, output columns, and validation remain intact.
Updating code can invalidate old Nextflow task caches; do not rerun that entrypoint
as part of this new run. Its published results remain independent.

The new adapter requires exact keys, Gene/Feature identity, allele_class, tier,
Consequence and `pcf_retained=true`. HC input additionally requires `LoF=HC`.
Duplicate exact alleles **within a candidate type** fail. Chromosome aliases
`chr22`/`22` are supported; source chromosome aliases must be unambiguous. Existing
source checks still reject multiallelic records at candidate positions, duplicate
exact source records, missing GT, invalid genotype allele indexes, or a bad index.

There is no new genotype QC. A called ALT index 1 is a carrier, including partial
and haploid ALT calls. GT/GQ/DP/AD/FT, site FILTER and ALT dosage are retained.
Homozygous ALT contributes one carrier variant record, with dosage two.

### Overlap, tiers and star records

The same exact allele can be present in both candidate types. The workflow queries
it once and writes **two annotation associations per carrier sample**, preserving
each type's Gene, Feature, SYMBOL and tier even when those differ. It never chooses
one annotation to overwrite the other. Counts across candidate types are therefore
**not additive**. The receipt reports overlap and distinct-allele counts explicitly.

All tier summaries are grouped by candidate type **and allele_class**. Null/empty
HC tiers become the explicit `untiered` group; these candidates are retained.
Star records stay in `spanning_deletion` audit strata; they are never added to a
sequence biological burden. No upstream-deletion mapping or event deduplication
is attempted. Raw carrier-annotation totals are labeled as such, not as a combined
biological burden. Within each type a homozygote counts once, not twice.

## Outputs

A new durable destination is used:

```text
<outdir>/filtered-carriers/<unit_id>/
  carriers.tsv.gz
  candidates.tsv
  unmatched.tsv
  samples.tsv
  sample_burden.tsv
  sample_gene_burden.tsv
  gene_burden.tsv
  sample_distinct_alleles.tsv
  receipt.json
```

- Carrier rows retain genotype fields and append `candidate_type` and `tier`.
- Candidate/unmatched audit rows retain type, tier, gene/transcript and allele class.
- `sample_burden.tsv` is long-form, grouped by sample/type/class/tier. **Every VCF
  sample is present for every observed stratum**, plus explicit untiered strata
  for both types/classes even with empty inputs. Zero-carrier samples remain zero.
- `sample_gene_burden.tsv` is sparse over nonzero sample/gene/type/class/tier cells.
- `gene_burden.tsv` includes candidate gene/type/class/tier groups with zero carriers.
- `sample_distinct_alleles.tsv` includes every sample and separately counts distinct
  sequence alleles and star records. A cross-type overlapping allele is counted
  once in its class; annotations in `carriers.tsv.gz` remain fully preserved.
- The receipt reports candidates/matched/unmatched/carriers and partial-call/dosage
  counts by type/class/tier, zero-carrier matched candidates, overlapping types,
  sample count, input paths/identities, code hashes and output hashes.

Parquet inputs and index are SHA-256 hashed. The large genotype VCF uses path,
size and mtime (no whole-VCF hashing/scanning). Input identities are checked before
and after extraction. Failed tasks remove partial products and retain a failed
receipt in work. Require a successful receipt and Nextflow exit, since older
published outputs can remain after a failed rerun. Only aggregate receipts should
leave the protected environment; other files contain individual-level data.

## Container check already completed

A synthetic Parquet write/read probe passed in the **existing pinned LOFTEE SIF**
on Expanse: DuckDB 1.5.5 and pysam 0.23.3. This was a dependency check, not NBDC
execution or a new resource audit. The NBDC profile reuses that image and existing
annotation-lock identity metadata. No new container/dependency installation is
needed. The pre-carrier resource lock is not needed to extract genotypes: the
filtered input hashes identify exactly which selected variants were consumed.

## Exact NBDC block12 commands

Update to this branch without discarding local changes. If checkout reports local
conflicts, keep those changes and resolve them; do not reset/clean the checkout.

```bash
cd "$HOME/rare-variant-pipeline"
git fetch origin
git switch --detach origin/feat/filtered-candidate-carriers
git rev-parse HEAD

BASE=/home/ood-guevara-james/abcd-fastvep-smoke
python3 - "$BASE" > "$PWD/abcd-block12-filtered-carriers.tsv" <<'PY'
from pathlib import Path
import sys
base=Path(sys.argv[1]);unit='chr22_block12'
p=base/'pre-carrier-nextflow'/'pre-carrier'/unit
vcf=Path('/shared/release/abcd/abcd/concatenated/genetics/sequencing/snv_indel/population_vcf/abcd_cohort_chr22_block12.vcf.gz')
paths=[p/'missense.filtered.parquet',p/'lof_hc.filtered.parquet',vcf,Path(str(vcf)+'.tbi')]
if not all(p.is_file() for p in paths):
    raise SystemExit('A filtered candidate, VCF or index input is missing')
print('unit_id\tchromosome\tmissense\tlof_hc\tvcf\tindex')
print('\t'.join([unit,'chr22',*map(str,paths)]))
PY

nextflow -C filtered_carriers.config run filtered_carriers.nf -profile nbdc \
  --carrier_manifest "$PWD/abcd-block12-filtered-carriers.tsv" \
  --resource_lock "$HOME/rare-variant-pipeline/abcd-annotation-resources.json" \
  --select_units chr22_block12 \
  --outdir "$BASE/filtered-carrier-nextflow" \
  -work-dir "$BASE/filtered-carrier-nextflow-work" \
  -with-trace "$BASE/filtered-carrier-block12.trace.tsv" -resume
```

The existing NBDC carrier profile defaults to 2 CPUs, 4 GB, 4h, medium partition
and QoS, queueSize 4. Existing `--carrier_*` resource overrides apply.

## Receipt gate before expansion

These assertions check **candidate counts**, not predicted carrier counts:

```bash
python3 - "$BASE/filtered-carrier-nextflow/filtered-carriers/chr22_block12/receipt.json" <<'PY'
import json,sys
r=json.load(open(sys.argv[1]));assert r['status']=='passed'
m=r['by_candidate_type']['missense'];h=r['by_candidate_type']['lof_hc']
assert m['candidate_records']==296 and h['candidate_records']==196
assert m['by_allele_class']['sequence']['candidate_records']==296
assert m['by_allele_class']['spanning_deletion']['candidate_records']==0
assert h['by_allele_class']['sequence']['candidate_records']==181
assert h['by_allele_class']['spanning_deletion']['candidate_records']==15
assert m['unmatched_candidates']==h['unmatched_candidates']==0
assert m['matched_candidates']==296 and h['matched_candidates']==196
assert r['samples_in_source']==8877
assert r['distinct_candidate_alleles']==492-r['overlapping_type_alleles']
for t in (m,h):
    for c in t['by_allele_class'].values():
        for field in ['candidate_records','matched_candidates','unmatched_candidates','carrier_records']:
            assert sum(x[field] for x in c['by_tier'].values())==c[field]
print(json.dumps({k:r[k] for k in ['samples_in_source','overlapping_type_alleles',
    'candidate_annotation_records','distinct_candidate_alleles','carrier_annotation_records',
    'distinct_allele_audit_by_class','by_candidate_type']},indent=2))
PY
```

Also inspect untiered HC counts, matched candidates with no carriers, the type
intersection, sequence/star tier separation, and source/input provenance. No exact
carrier count is asserted in advance. Repeat the command with a different trace
filename, same work/launch paths and `-resume`; block12 should be `CACHED`.

Then add explicit manifest rows for blocks 0 and 19 and run `--select_units
chr22_block0,chr22_block12,chr22_block19`. Apply the block12 count assertions only
to block12. Keep output/work paths and execution parameters stable so block12
remains cached. Optional `--expected_missense`/`--expected_hc` flags are per-task
guards and change task cache keys; the commands above deliberately use receipt
assertions instead, allowing subset expansion without changing those parameters. After the three-block gate passes, extend the
manifest to all 23 blocks and use `--select_units all`. Inspect every receipt;
do not sum candidate-type carrier totals into a combined burden. No such NBDC
runs have been launched by Codex.

## Local validation

Synthetic tests cover filtered input selection, tier/gene/transcript propagation,
untiered HC, both index types, exact REF/ALT/star matching, quality and dosage
preservation, zero-carrier samples, empty inputs, invalid/unfiltered inputs,
cross-type overlaps with different annotations, and real Nextflow subset/resume
and durable publication. Existing HC-only regression tests also pass.

```bash
uv run --with pytest --with duckdb==1.5.5 --with pysam python -m pytest -q \
  tests/test_filtered_carriers.py tests/test_exact_carriers.py
```

Result: **30 tests passed** (12 new and 18 existing). NBDC config rendering and
Python compilation are checked locally. Real NBDC carrier counts remain pending.
