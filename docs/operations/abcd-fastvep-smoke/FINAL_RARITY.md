# Final rarity from saved extraction products

Code pin: `3edbec9613cb05ad564961b3a022f55f4d8f9d6b` on `feat/final-rarity`,
stacked on PR #19. Entrypoint: `final_rarity.nf`, config: `final_rarity.config`.

The operator reports 951 autosomal and 25 X/Y extraction/frequency receipts passed
on NBDC, the X/Y pilot resumed from cache, and X/Y reported no zero-AN alleles.
These are operator-supplied results, not new NBDC validation by Codex.

## Fixed eligibility and products

For each exact CHROM/POS/REF/ALT, retain only a matched allele with saved
`unrelated_an > 0` and **saved `unrelated_af < 0.001`**. Exactly 0.001 fails.
Missing, invalid, nonfinite, negative, or out-of-range AF fails closed; missing,
invalid or zero AN fails closed. Unmatched alleles fail. No source INFO AF,
cohort AF, or recalculated genotype frequency is substituted. `22` and `chr22`
are equivalent for joining; output source fields remain preserved.

One decision applies to every candidate annotation and carrier row of an allele.
Raw genotype/quality/dosage fields are untouched. Types, tiers, untiered HC and
stars are preserved; star records remain separate in all burden/audit counts.
Passing variants with no observed carriers remain in the retained candidate and
frequency tables. All source samples, including zero-carrier samples, remain in
sample summaries. No genotype QC or downstream gathering is invoked.

Inputs per unit are only `carriers.tsv.gz`, `candidates.tsv`,
`variant_frequencies.tsv`, `samples.tsv`, and the passed extraction `receipt.json`.
All four product hashes are checked against that receipt; sample, candidate and
carrier counts must reconcile. Candidate and frequency exact-allele sets must
agree. Duplicate alleles/associations, annotation inconsistencies, conflicting
cross-type genotypes and changed inputs fail clearly instead of silently losing
rows. Genotype VCFs, annotation resources and PSAM are not read by this stage.

Durable products under `<outdir>/final-rarity/<unit_id>/`:

- `carriers.rare.tsv.gz`: passing raw carrier annotations, same columns/values.
- `candidates.rare.tsv`, `variant_frequencies.rare.tsv`: retained variant products.
- `rarity_audit.tsv`: all distinct frequency alleles with pass/fail and one reason.
- `samples.tsv`, `sample_burden.tsv`, `sample_gene_burden.tsv`, `gene_burden.tsv`,
  `sample_distinct_alleles.tsv`: recomputed summaries; distinct allele/sample
  counts avoid cross-type double counting, separately by allele class.
- `receipt.json`: input identities/source frequency policy, reconciled allele,
  candidate-annotation and carrier-annotation pass/fail counts; breakdowns by
  type/class/tier; per-class reason counts and distinct allele/sample counts;
  code/output identities and stage status.

Keep variant/sample tables local. Share aggregate receipts only. Raw extraction
outputs remain unchanged. Require a successful task exit and receipt; use the
new destination below so old products are not mistaken for this stage's outputs.

## Pinned NBDC block12, resume, then expansion

The NBDC profile requests **1 CPU/task**, at most **16 concurrent tasks**, and
eligible Slurm partitions `small,medium`. Trace overwrite is enabled. No `medium`
QoS is forced; use `--rarity_qos` only if required by your established allocation.
The scheduler's account/QoS limits still govern actual running concurrency.
Defaults are 4 GB/task and one hour, configurable by `--rarity_memory` and
`--rarity_time`. The existing validated LOFTEE SIF/resource lock is reused; the
stage needs only Python's standard library and no new container.

Generate the manifest once from the two supplied roots. This inventories saved
files/receipts and does not invoke upstream workflows. Existing manifests are
not overwritten; use the same manifest for resume.

```bash
set -euo pipefail
cd "$HOME/rare-variant-pipeline"
git fetch origin feat/final-rarity
git switch --detach 3edbec9613cb05ad564961b3a022f55f4d8f9d6b
BASE=/home/ood-guevara-james/abcd-fastvep-smoke
RUN="$BASE/final-rarity-nextflow"
mkdir -p "$RUN"
python3 scripts/make_final_rarity_manifest.py \
  --root "$BASE/carrier-af005-with-frequencies/output/filtered-carriers" \
  --root "$BASE/carrier-xy-grch38-x-only-par/output/filtered-carriers" \
  --output "$RUN/manifest.tsv"
python3 - "$RUN/manifest.tsv" <<'PY'
import csv, sys
r=list(csv.DictReader(open(sys.argv[1]),delimiter='\t'))
assert len(r)==976
assert sum(x['chromosome'] not in ('chrX','chrY') for x in r)==951
assert sum(x['chromosome'] in ('chrX','chrY') for x in r)==25
print('Manifest: 951 autosomal + 25 X/Y saved units')
PY
run_rarity() {
  nextflow -C final_rarity.config run final_rarity.nf -profile nbdc \
    --rarity_manifest "$RUN/manifest.tsv" \
    --resource_lock "$HOME/rare-variant-pipeline/abcd-annotation-resources.json" \
    --select_units "$1" --outdir "$RUN/output" \
    -work-dir "$RUN/work" -with-trace "$RUN/trace.tsv" -resume
}
run_rarity chr22_block12
```

Inspect its aggregate receipt and reconcile it with the saved extraction counts.
No fixed retained count is predicted:

```bash
python3 - "$RUN/output/final-rarity/chr22_block12/receipt.json" <<'PY'
import json, sys
r=json.load(open(sys.argv[1]))
assert r['status']=='passed' and r['stage']=='final_rarity'
assert r['samples_in_source']==8877 and r['reconciliation']['passed']
for k in ['distinct_alleles','candidate_annotations','carrier_annotations']:
    c=r[k]; assert c['input']==c['pass_count']+c['fail_count']
for c in r['by_allele_class'].values():
    assert sum(c['reasons'].values())==c['distinct_alleles']['input']
print(json.dumps({k:r[k] for k in ['samples_in_source','distinct_alleles',
    'candidate_annotations','carrier_annotations','by_candidate_type','by_allele_class']},indent=2))
PY
run_rarity chr22_block12
```

The repeated command overwrites `trace.tsv`; require `CACHED` for block12. Then
validate blocks 0/12/19 and one existing X and Y unit from the manifest:

```bash
PILOTS=$(python3 - "$RUN/manifest.tsv" <<'PY'
import csv, sys
rows=list(csv.DictReader(open(sys.argv[1]),delimiter='\t'))
units=['chr22_block0','chr22_block12','chr22_block19']
units += [next(r['unit_id'] for r in rows if r['chromosome']==c) for c in ['chrX','chrY']]
print(','.join(units))
PY
)
run_rarity "$PILOTS"
```

Inspect those receipts with the same reconciliation rules, especially separate
star counts and zero-carrier samples. After the pilot passes, use `run_rarity all`
to apply **only final rarity** to the 976 saved units. Previously completed final
rarity tasks should cache. No extraction, annotation, frequency counting or QC
stage is launched by these commands. Inspect every receipt after successful
Nextflow completion; the manifest contains all expected units but does not imply
that all final-rarity tasks have already completed.

## X/Y downstream QC

The [post-rarity QC adapter](POST_RARITY_QC.md) now implements the requirements
below as a separate entrypoint. The description below records the legacy
entrypoint limitations; do not wire its chromosome-blind QC directly to X/Y.

### Requirements identified during final-rarity implementation

The current `scripts/qc_filtered_carriers.py` is chromosome-blind. Its `HET` set
accepts `0/1`/phased heterozygotes, and its `HOM` set accepts both diploid ALT
homozygotes and literal `1`, regardless of sex or region. It has no PSAM input.
Final rarity does **not** fix this, because frequency exclusions never removed
raw carrier rows for otherwise rare alleles.

A subsequent QC change must:

1. Join source samples to the existing PSAM and assign expected ploidy by the same
   explicit GRCh38 X-only PAR policy. Keep all source samples in summaries.
2. Exclude/audit unknown-sex X/Y calls and female Y calls. These can currently pass
   ordinary quality/AB checks, and remain present after final rarity.
3. For known-sex X PAR and female X non-PAR, require diploid calls. Literal `1`
   currently passes QC but is unexpected ploidy in these regions.
4. For male X non-PAR/Y non-PAR, accept native haploid ALT and homozygous diploid
   encodings, exclude/audit heterozygous calls even if AB is 0.25–0.75, and retain
   the agreed partial-call exclusion. Apply the existing ALT AB >=0.90 rule to
   the eligible haploid ALT calls; do not alter thresholds here.
5. Preserve raw GT/dosage for provenance and explicitly distinguish effective
   haploid dosage one from the raw `1/1` dosage two in post-QC allele totals.
   Homozygous ALT still counts as one carrier variant. Keep stars separate.
6. Add an explicit input adapter for `carriers.rare.tsv.gz` and the final-rarity
   receipt. Current QC expects extraction output names/hashes and counts; simply
   renaming files or passing this receipt is insufficient. Update subsequent QC
   provenance/summary handling for the ploidy policy before gathering.

The existing GQ >=20, DP >=10, site PASS, AD/AB and FT behavior is otherwise
unchanged. No QC code or policy was modified by this PR. This stage reuses only
its table-writing and summary helpers; it never calls genotype evaluation.

## Validation

28 local tests passed: strict threshold boundaries, missing/nonfinite AF, zero
and malformed AN, exact allele consistency, cross-type propagation, input hash
checks, immutable inputs, stars/tiers/untiered HC, empty/all-fail inputs, all
samples including zeros, X/Y raw-call preservation, and actual local Nextflow
execution/resume with trace overwrite after deleting genotype VCF inputs.
NBDC execution and large-cohort runtime/memory remain untested by Codex.
