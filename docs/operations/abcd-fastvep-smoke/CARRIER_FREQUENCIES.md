# Corrected autosomal frequencies during carrier extraction

Code pin: `89a7dfd72a0c1c0657ad25e5251e23afa3ea6594`.
Branch: `feat/carrier-corrected-frequencies`, stacked on the batching/preliminary
screen update. This adds counting; it does not apply the final rarity filter.

## Agreed policy

Enable explicitly with `--compute_frequencies true --psam /path/to/samples.psam`.
The default remains raw-only extraction. The existing production SIF and one-CPU,
10 kb pysam batching are unchanged. Counting takes place on each exact matched
record already fetched for extraction; there is no second genotype VCF scan.

- Cohort: `frequency_representative=1`.
- Unrelated: `frequency_representative=1 AND unrelated=1`.
- At most one selected representative per participant among source VCF samples;
  duplicate representatives fail. Participants without a selected representative
  are excluded from frequency sets and counted in the receipt. Representatives
  are never automatically chosen. Extra PSAM rows are allowed and not counted.
- All source samples remain in raw carriers and sample/burden summaries, regardless
  of flags. Sample order and zero-carrier samples are preserved.
- **No additional genotype QC:** no GQ, DP, AD, AB, FT, or site FILTER exclusions
  from frequency counts. Label: `representative-autosomal-no-genotype-qc-v1`.
  Existing downstream carrier QC is unchanged and separate.
- Autosomal expected ploidy is two. Complete and partial diploid calls count
  called alleles, including reference calls: `0/0` gives AC=0/AN=2, `0/.` gives
  AC=0/AN=1, `1/.` gives AC=1/AN=1. Phasing does not affect counts. Fully missing
  calls contribute zero; unexpected ploidy is excluded and audited.
- Unknown SEX values do not exclude autosomal calls. SEX is not interpreted by
  this autosomal implementation.
- AN=0 gives missing AF (`.` in TSV), not zero. Unmatched candidates have missing
  AC, AN and AF, distinct from matched sites with AN=0.
- Source INFO AC/AN/AF and their AC/AN ratio are preserved separately. Corrected
  counts are from the selected samples' GTs, not from INFO or carrier rows.
- Each exact allele is counted once even if it has both missense and HC-LoF
  annotations. Sequence and star records remain separate; no deletion-event
  inference, tier changes, or carrier deduplication policy changes are introduced.

The pure counting helper also tests the agreed expected-haploid rules: literal
`0` **or** `1` counts one allele; `0/0` and `1/1` collapse to one allele;
heterozygous `0/1` and partial calls are excluded and audited. These rules are
**not yet applied to X/Y records** because their expected ploidy/PAR assignment
is unresolved. Invalid allele indices are excluded by the counting primitive;
VCF decoding otherwise follows the existing pysam behavior and exact-source
validation, including rejection of multiallelic source records.

## Explicit X/Y boundary

With counting enabled, Python and the Nextflow entrypoint reject X/Y before
extraction. No PAR intervals, X/Y duplicate-locus treatment or sex-code mapping
are guessed. Confirm ABCD genome build, PAR location on X/Y, and whether equivalent
PAR records occur on both before enabling sex-aware counting. The agreed unknown
sex rule will exclude/audit X/Y when that path is added. Raw X/Y extraction remains
available with `--compute_frequencies false`.

The user agreed to discuss the PAR representation later. This restriction does
not block autosome counting and does not certify any X/Y frequency results.

## Products and reconciliation

All existing raw products retain their schemas and values. New durable products:

- `variant_frequencies.tsv`: one row per distinct input candidate allele, with
  exact key, allele class, matched status, original source INFO fields/ratio,
  and `cohort_ac/an/af`, `unrelated_ac/an/af`. Each set also has
  `counted_genotypes`, `reference_genotypes`, and `excluded_genotypes`.
  A reference genotype here means a contributing call with no called ALT,
  including `0/.`; partial calls with a called ALT are not reference genotypes.
- `frequency_audit.tsv`: counts per exact allele, class, sample set, and mutually
  exclusive counting/exclusion reason. Contains variant identities; keep local.
- `receipt.json`: `frequencies.policy`, sample-set sizes/selection counts, and
  matched distinct alleles, reason counts and zero-AN counts separately by allele
  class and sample set. Counts must reconcile to matched alleles × selected samples.
  It retains PSAM/input/code/output hashes and reports explicit frequency status.

`sample_metadata.tsv` retains the PSAM linkage. No sample-level or variant-level
output should be pasted into shared logs. Aggregate receipts can be shared.

## Final rarity is a separate PR

The intended downstream rule is **corrected unrelated AF < 0.001**, with missing
AF unable to pass. It is not applied here: raw outputs still include candidates
that fail that final threshold. That future filter can join `variant_frequencies.tsv`
to candidates/carriers by exact allele without repeating extraction or reading
reference genotypes again. It must recompute downstream summaries after selection.

This revision supersedes the earlier deferral of autosomal counting in
[CARRIER_BATCHING.md](CARRIER_BATCHING.md). PAR/X/Y counting remains deferred.

## Pinned NBDC block12 command

Reuse completed annotation/scoring and the **new source-AF <0.005** pre-carrier
outputs from the [batching runbook](CARRIER_BATCHING.md). Do not substitute the
older 0.001-screened files. The existing `$BASE/carrier-source-af005/carriers.tsv`
manifest names those filtered Parquets and indexed source VCF. PSAM must be
prepared on NBDC; its actual path is supplied explicitly below, not guessed.

```bash
set -euo pipefail
cd "$HOME/rare-variant-pipeline"
git fetch origin feat/carrier-corrected-frequencies
git switch --detach 89a7dfd72a0c1c0657ad25e5251e23afa3ea6594
BASE=/home/ood-guevara-james/abcd-fastvep-smoke
# Set this to the actual prepared NBDC PSAM before running the remaining commands.
: "${PSAM:?Set PSAM to the absolute path of the prepared tab-delimited PSAM}"
test -f "$PSAM"
RUN="$BASE/carrier-af005-with-frequencies"
mkdir -p "$RUN"
nextflow -C filtered_carriers.config run filtered_carriers.nf -profile nbdc \
  --carrier_manifest "$BASE/carrier-source-af005/carriers.tsv" \
  --resource_lock "$HOME/rare-variant-pipeline/abcd-annotation-resources.json" \
  --psam "$PSAM" --compute_frequencies true \
  --carrier_cpus 1 --carrier_batch_bp 10000 --select_units chr22_block12 \
  --outdir "$RUN/output" -work-dir "$RUN/work" \
  -with-trace "$RUN/block12.trace.tsv" -resume
```

The first run computes frequencies and raw carriers together in this new output
location; earlier raw-only tasks cannot supply the reference genotypes needed
for AN. Existing annotation, scoring and pre-carrier products are reused. No
production image rebuild or fresh resource lock is needed.

Check aggregate receipts after a successful Nextflow exit:

```bash
python3 - "$RUN/output/filtered-carriers/chr22_block12/receipt.json" \
  "$BASE/carrier-source-af005/pre-carrier-output/pre-carrier/chr22_block12/receipt.json" <<'PY'
import json, sys
r, pre = [json.load(open(p)) for p in sys.argv[1:]]
assert r['status'] == pre['status'] == 'passed'
assert '< 0.005' in pre['policy']['cohort']
assert r['samples_in_source'] == 8877 and r['query_batch_bp'] == 10000
assert r['frequency_status'].startswith('computed_autosomes_no_genotype_qc')
f = r['frequencies']; assert f['status'] == 'passed'
assert f['sample_selection']['unrelated'] <= f['sample_selection']['cohort'] <= 8877
assert f['policy']['unrelated'] == 'frequency_representative == 1 AND unrelated == 1'
for kind, c in r['by_candidate_type'].items():
    assert c['candidate_records'] == pre['counts'][kind]['final_retained']
    assert c['unmatched_candidates'] == 0
for cls in f['by_allele_class'].values():
    for name, counts in cls['sample_sets'].items():
        assert sum(counts['by_reason'].values()) == counts['evaluated_genotypes']
        assert counts['evaluated_genotypes'] == counts['expected_evaluated_genotypes']
        assert counts['expected_evaluated_genotypes'] == cls['matched_distinct_alleles'] * f['sample_selection'][name]
print(json.dumps({'frequency_status': r['frequency_status'], 'frequencies': f}, indent=2))
PY
```

Inspect representative counts, participants without representatives, unexpected
ploidy exclusions, zero-AN counts, and star separation. A successful receipt
alone does not establish that a user-prepared PSAM encodes the intended sample set.
No real candidate, carrier, representative or unrelated count is predicted here.

Repeat with a different trace filename and identical parameters/work directory
plus `-resume`; require CACHED. Changing PSAM must rerun counting. After block12
and resume pass, expand the explicit manifests to blocks 0/12/19, then all desired
autosomal blocks. Keep X/Y out of the counting selection until PAR is resolved.

## Validation limits

Synthetic tests cover called reference/ALT alleles, partial/missing calls,
phasing, sample-set intersection, duplicate representatives, no quality-field
access by the counter, zero AN, unmatched/empty candidates, original carrier
parity, overlapping annotations, separate stars, both index formats, X/Y blocking,
and real local Nextflow resume/PSAM cache invalidation/durable publication.
Haploid policy is unit-tested without assigning source X/Y ploidy.

No NBDC job or real-data frequency validation was performed by Codex. The existing
container dependencies suffice, but this revision has not been executed in the
production SIF on NBDC. Final rarity filtering and X/Y/PAR counting are not included.
