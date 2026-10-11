# Pysam batching and preliminary frequency screen

Current validation status: the operator reports completed NBDC validation through
976-block/24-chromosome gathering. See the [finalization record](../../validation/abcd-post-rarity-finalization.md).
Historical code pins and stage-specific instructions below remain unchanged.

For the subsequent opt-in autosomal frequency-counting implementation, see
[CARRIER_FREQUENCIES.md](CARRIER_FREQUENCIES.md). The pinned commands and deferral
status below describe the earlier batching revision.

## Scope and policy

Production code pin: `ea937a70125292d35a6d5a8f50fe5cb361c151c7`.
Branch: `feat/carrier-batching-frequency-policy`.

The operator reported block12 full-output parity and wall time of 321.6 s versus
68.1 s for pysam with 10 kb batching. These are operator-provided observations,
not timings obtained by Codex or a guarantee for this production revision.

This update keeps pysam and the existing validated LOFTEE SIF. Filtered extraction
uses one CPU per block and combines candidate positions into nonoverlapping query
spans of at most 10,000 bases, anchored at the first candidate in each span.
Exact CHROM/POS/REF/ALT matching and duplicate detection remain mandatory. Other
positions read within a batch cannot become candidates. The older HC-only
entrypoint retains its default single-position query behavior.

Execution order:

1. Existing scored candidate Parquets and sites-only VCFs remain inputs.
2. Preliminary filter: **uncorrected source INFO/AC / INFO/AN < 0.005**; gnomAD
   `gnomAD4.1_joint_POPMAX_AF < 0.001` or missing; existing region exclusions.
3. Batched raw carrier extraction. All source samples remain represented.
4. Corrected cohort/unrelated frequencies and final unrelated AF < 0.001 are
   **deferred**, following the user's clarification. Do not label these outputs
   as passing final corrected rarity. No new frequency genotype QC is imposed.

X and Y are accepted for preliminary filtering and raw extraction where inputs
and chromosome resources exist. This does not imply that sex-aware frequency
counting or sex-chromosome burdens have been validated.

## Preserved provenance and optional PSAM

The pre-carrier Parquets retain original AC/AN/AF strings in
`pcf_source_info_ac`, `pcf_source_info_an`, `pcf_source_info_af`; `pcf_cohort_af`
continues to mean their valid AC/AN ratio. INFO/AF and MLEAF never drive selection.

Extraction adds `source_frequencies.tsv`: one row per distinct candidate allele,
with matched status, allele class, source INFO AC/AN/AF and their AC/AN ratio.
Missing/unmatched values are `.`. These values come from the genotype VCF's INFO,
not a genotype recount. Pysam formats parsed numeric values; original lexical
INFO strings are preserved in the pre-carrier Parquets. Existing carrier and
burden table schemas stay unchanged, including separate star/tier counts.

Optional `--psam /absolute/path/samples.psam` accepts a tab-delimited table:

```text
#IID    SEX    participant_id    frequency_representative    unrelated
```

Use tabs, unique exact VCF sample IDs, nonempty participant IDs, and explicit
`0`/`1` flags. SEX is preserved verbatim, including an explicit unknown code;
it is not interpreted. Every VCF sample must have a row. Extra PSAM rows and
columns are allowed. `sample_metadata.tsv` retains metadata in VCF sample order,
including zero-carrier samples. No rows are excluded according to either flag,
no intersection is applied, and representative uniqueness is not yet enforced.
The receipt records the PSAM SHA-256 and aggregate metadata counts.

## Decisions deferred before corrected frequency implementation

- The representative and unrelated sets must be defined independently; do not
  assume their intersection. Establish duplicate-participant validation for each.
- Define genome-build-specific X/Y PAR intervals and how PAR records on X and Y
  are represented so the same biological locus is not counted twice.
- Define sex/ploidy expectations and conversion of diploid-encoded haploid calls.
- Decide unknown-sex and unexpected heterozygous-call treatment.
- Decide whether partially called GTs contribute their called alleles or no
  alleles. Raw extraction currently retains any called ALT, unchanged.
- No GQ/DP/AB filter is implicitly part of frequency counting. Any future rule
  must be explicit and distinguish reference and ALT calls. Existing downstream
  carrier QC remains a separate stage.
- Define zero-AN handling before applying final corrected unrelated rarity.

Corrected AN cannot be reconstructed from carrier-only tables. Because counting
is deferred, that later step must read eligible reference and ALT genotypes from
the source VCF at retained sites (or use an explicitly designed sufficient-count
product). Source INFO AN is not a corrected denominator.

## Pinned NBDC block12 pilot

Do not rerun annotation or scoring. Reuse the existing resource locks and SIF.
Use new output/work directories: older filtered candidates were already censored
at 0.001, so extraction alone cannot recover variants newly eligible below 0.005.
No old 296/196 candidate-count assertion applies to this relaxed screen.

```bash
set -euo pipefail
cd "$HOME/rare-variant-pipeline"
git fetch origin feat/carrier-batching-frequency-policy
git switch --detach ea937a70125292d35a6d5a8f50fe5cb361c151c7
BASE=/home/ood-guevara-james/abcd-fastvep-smoke
RUN="$BASE/carrier-source-af005"
mkdir -p "$RUN"
python3 - "$BASE" "$RUN" <<'PY'
import csv, sys
from pathlib import Path
base, run = map(Path, sys.argv[1:])
u = 'chr22_block12'
c = base/'candidate-nextflow/candidates'/u
sites = base/'sites-catalog-smoke/sites'/f'{u}.sites.vcf.gz'
paths = [c/'missense.parquet', c/'lof_hc.parquet', sites]
if not all(p.is_file() for p in paths):
    raise SystemExit('Adjust manifest paths to the existing published candidates/sites')
with (run/'pre-carrier.tsv').open('w') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['unit_id','chromosome','missense','lof_hc','sites'])
    w.writerow([u,'chr22',*map(str,paths)])
PY
nextflow -C pre_carrier.config run pre_carrier.nf -profile nbdc \
  --filter_manifest "$RUN/pre-carrier.tsv" \
  --filter_resource_lock "$HOME/rare-variant-pipeline/abcd-pre-carrier-resources.json" \
  --select_units chr22_block12 --outdir "$RUN/pre-carrier-output" \
  -work-dir "$RUN/pre-carrier-work" -with-trace "$RUN/pre-carrier.trace.tsv" -resume
```

Inspect the successful preliminary receipt, then create the extraction manifest:

```bash
python3 - "$RUN" <<'PY'
import csv, json, sys
from pathlib import Path
run = Path(sys.argv[1]); u = 'chr22_block12'
p = run/'pre-carrier-output/pre-carrier'/u
r = json.loads((p/'receipt.json').read_text())
assert r['status'] == 'passed' and '< 0.005' in r['policy']['cohort']
print(json.dumps(r['counts'], indent=2))
v = Path('/shared/release/abcd/abcd/concatenated/genetics/sequencing/snv_indel/population_vcf')/f'abcd_cohort_{u}.vcf.gz'
paths = [p/'missense.filtered.parquet', p/'lof_hc.filtered.parquet', v, Path(str(v)+'.tbi')]
assert all(x.is_file() for x in paths), 'Missing carrier input'
with (run/'carriers.tsv').open('w') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['unit_id','chromosome','missense','lof_hc','vcf','index'])
    w.writerow([u,'chr22',*map(str,paths)])
PY
nextflow -C filtered_carriers.config run filtered_carriers.nf -profile nbdc \
  --carrier_manifest "$RUN/carriers.tsv" \
  --resource_lock "$HOME/rare-variant-pipeline/abcd-annotation-resources.json" \
  --carrier_cpus 1 --carrier_batch_bp 10000 \
  --select_units chr22_block12 --outdir "$RUN/carrier-output" \
  -work-dir "$RUN/carrier-work" -with-trace "$RUN/carrier.trace.tsv" -resume
```

If PSAM is ready, append `--psam /absolute/path/to/the/prepared.psam` to the
extraction command before its first run. Omit it otherwise; corrected frequencies
are not enabled by supplying it. Adding/changing PSAM changes the task cache.

## Validation and expansion

- Require passed receipts and zero unmatched candidates. Reconcile extraction
  candidate counts with this run's pre-carrier `final_retained` counts by type.
- Confirm `query_batch_bp=10000`, `frequency_status` reports deferred, and all
  8,877 source samples remain. Inspect separated sequence/star/tier counts.
- Repeat each command with a new trace filename and the same code, inputs,
  parameters, work directory and `-resume`; require CACHED.
- Then add explicit block0/19 manifest rows and select those plus block12.
  After this gate, extend the manifests and select `all`. This update does not
  guess unpublished manifests or launch NBDC jobs.
- Completed annotation/scoring outputs remain reusable. New filtered candidates
  require new extraction outputs; old results remain in their original locations.
- Keep all sample/variant tables local. Share aggregate receipts and timings only.

Local synthetic validation covers full-output batching parity, query boundaries,
duplicate/multiallelic safeguards, threshold boundaries, INFO provenance, X/Y raw
extraction, PSAM coverage/order/no-subsetting, and real local Nextflow subset,
resume, and durable publication with and without PSAM. The existing SIF has the
required pysam/DuckDB dependencies; no image was changed. Production NBDC execution
of this revision and its performance remain untested.
