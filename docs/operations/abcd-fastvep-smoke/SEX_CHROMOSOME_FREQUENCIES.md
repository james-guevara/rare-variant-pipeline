# GRCh38 X-only PAR frequency counting

Current validation status: the operator reports completed NBDC validation through
976-block/24-chromosome gathering. See the [finalization record](../../validation/abcd-post-rarity-finalization.md).
Historical code pins and stage-specific instructions below remain unchanged.

Code pin: `a1cec12a6404afe120dc79563be77146277d69b6` on
`feat/carrier-corrected-frequencies` (PR #19, including PR #18).

Enable with `--compute_frequencies true --sex_chromosome_policy grch38_x_only_par`
and the **existing PSAM**. The flag declares the supplied source representation
and genome build; the workflow does not infer or discover PAR policy.

## Policy

| Region (1-based inclusive, variant POS) | Male | Female | Unknown sex |
|---|---|---|---|
| X PAR1: 10001–2781479 | diploid | diploid | exclude/audit |
| X PAR2: 155701383–156030895 | diploid | diploid | exclude/audit |
| X outside PAR | haploid | diploid | exclude/audit |
| Y outside PAR | haploid | exclude/audit | exclude/audit |
| Y PAR1: 10001–2781479; PAR2: 56887903–57217415 | **fail the task for any candidate** | same | same |

PSAM `SEX=1` means male and `SEX=2` means female. Other nonempty values mean
unknown, including `0`, `-9`, `NA`, and `.`. Text `M`/`F` is not silently interpreted
as numeric PSAM encoding. Inspect the aggregate `sample_selection.by_sex` counts
against the existing PSAM's intended encoding before expanding.

For expected haploidy, literal GT `0` or GT `1` contributes one allele; `0/0` and
`1/1` (also phased) collapse to one allele. Heterozygous and partial calls are
excluded with distinct audit reasons. Fully missing calls contribute nothing.
Expected diploid calls retain the prior rule: count called alleles from partial
calls. There are **no GQ/DP/AD/AB, FT, or site FILTER exclusions** from frequencies.

Cohort remains `frequency_representative=1`; unrelated remains representative
**AND** unrelated. Reference genotypes contribute AN. Raw carriers, raw dosage,
source INFO frequencies/ratio, tiers, untiered HC and separate star counts remain
unchanged. Female/unknown-sex carriers are not removed from raw tables. Final AF
filtering is not applied. AN=0 still produces missing AF.

Assignment is by POS, not the reference allele's full span. Candidate validation
rejects Y-PAR before opening the genotype VCF, even when the candidate would be
unmatched. It does not scan noncandidate Y records for a separate whole-callset
PAR audit. Raw extraction without frequency counting is unaffected by this guard.

## Evidence supplied by the operator

NBDC checks reported on 2026-10-10:

- PAR1: 399,739 records on X; zero on Y.
- PAR2: 29,044 records on X; zero on Y.
- Twenty SNPs sampled per X PAR: complete diploid calls in both sexes, including
  heterozygotes in both; male/female median DP 38/38 (PAR1), 36/36 (PAR2).
- All 951 autosomal blocks completed and passed receipts.

These are operator-reported observations, not new checks run by Codex. Completed
autosomal products should remain as they are. The autosomal frequency policy,
counts and output schemas are unchanged. Updated staged code can change Nextflow
cache keys: **do not submit the completed autosomal manifest on this revision**.

## Outputs and audit

For X/Y, `variant_frequencies.tsv` adds `frequency_region` (`X_PAR`, `X_nonPAR`,
`Y_nonPAR`). Autosomal table schemas stay unchanged. Per-class receipt counts add
`matched_by_region`, and selection metadata adds representative counts by sex.
The new policy version is `representative-grch38-x-only-par-no-genotype-qc-v1`.

Audit reasons include `excluded_unknown_sex`, `excluded_female_Y`,
`excluded_haploid_partial`, `excluded_haploid_heterozygous`,
`counted_haploid`, and `counted_diploid_encoded_haploid`, alongside existing
reasons. Receipt `eligible_samples` means representatives selected before
per-site sex/ploidy exclusions. All reasons reconcile to matched distinct alleles
× selected representatives, separately by sample set and sequence/star class.

## NBDC: prepare an X/Y-only manifest

Reuse existing published X/Y candidates screened with source AC/AN <0.005 and
existing source VCF/index paths. No annotation/scoring/autosomal extraction is
launched here. Set `PSAM` to your existing PSAM and `SOURCE_MANIFEST` to a carrier
manifest containing X/Y rows (it may also contain autosomes). If it has no X/Y
rows, prepare those rows with the same six-column carrier schema first.

```bash
set -euo pipefail
cd "$HOME/rare-variant-pipeline"
git fetch origin feat/carrier-corrected-frequencies
git switch --detach a1cec12a6404afe120dc79563be77146277d69b6
BASE=/home/ood-guevara-james/abcd-fastvep-smoke
: "${PSAM:?Set PSAM to the existing absolute PSAM path}"
: "${SOURCE_MANIFEST:?Set SOURCE_MANIFEST to a carrier manifest containing prepared X/Y inputs}"
test -f "$PSAM"
RUN="$BASE/carrier-xy-grch38-x-only-par"
mkdir -p "$RUN"
python3 - "$SOURCE_MANIFEST" "$RUN" <<'PY'
import csv, sys
from pathlib import Path
source, out = Path(sys.argv[1]).resolve(), Path(sys.argv[2])
with source.open() as f:
    reader = csv.DictReader(f, delimiter='\t')
    header = reader.fieldnames
    required = ['unit_id','chromosome','missense','lof_hc','vcf','index']
    assert header and set(required) <= set(header), 'Invalid carrier manifest header'
    selected = []
    for row in reader:
        c = row['chromosome'].removeprefix('chr')
        if c not in ('X','Y'): continue
        row['chromosome'] = 'chr'+c
        for k in ['missense','lof_hc','vcf','index']:
            p = Path(row[k]); p = p if p.is_absolute() else source.parent/p
            assert p.is_file(), 'Missing X/Y input; prepare the corresponding manifest row'
            row[k] = str(p.resolve())
        selected.append(row)
assert selected, 'No X/Y rows; completed autosomes must not be rerun'
with (out/'xy-carriers.tsv').open('w') as f:
    writer = csv.DictWriter(f, fieldnames=header, delimiter='\t', lineterminator='\n')
    writer.writeheader(); writer.writerows(selected)
pilots = [next(r['unit_id'] for r in selected if r['chromosome']==c)
          for c in ['chrX','chrY'] if any(r['chromosome']==c for r in selected)]
(out/'pilot-units.txt').write_text(','.join(pilots)+'\n')
print({c: sum(r['chromosome']==c for r in selected) for c in ['chrX','chrY']})
PY
PILOT_UNITS=$(cat "$RUN/pilot-units.txt")
nextflow -C filtered_carriers.config run filtered_carriers.nf -profile nbdc \
  --carrier_manifest "$RUN/xy-carriers.tsv" \
  --resource_lock "$HOME/rare-variant-pipeline/abcd-annotation-resources.json" \
  --psam "$PSAM" --compute_frequencies true \
  --sex_chromosome_policy grch38_x_only_par \
  --carrier_cpus 1 --carrier_batch_bp 10000 --select_units "$PILOT_UNITS" \
  --outdir "$RUN/output" -work-dir "$RUN/work" \
  -with-trace "$RUN/pilot.trace.tsv" -resume
```

The helper picks the first listed block of each available sex chromosome for the
pilot; it does not inspect variants to choose blocks. If you know a block spanning
X PAR, include it in `PILOT_UNITS` for a real PAR check. The source manifest, PSAM
and completed autosomal outputs are read-only inputs to these instructions.
Existing validated SIF/resource lock is reused; no container rebuild is required.

## Pilot, resume, then X/Y expansion

Inspect aggregate receipts only:

```bash
python3 - "$RUN/output/filtered-carriers" <<'PY'
import json, sys
from pathlib import Path
paths = sorted(Path(sys.argv[1]).glob('*/receipt.json'))
assert paths, 'No published receipts'
for p in paths:
    r = json.loads(p.read_text()); f = r['frequencies']
    assert r['status'] == f['status'] == 'passed'
    assert r['chromosome'] in ('chrX','chrY')
    assert f['policy']['sex_chromosomes'] == 'grch38_x_only_par'
    assert f['policy']['genome_build'] == 'GRCh38'
    assert all(c['unmatched_candidates']==0 for c in r['by_candidate_type'].values())
    for cls in f['by_allele_class'].values():
        n = cls['matched_distinct_alleles']
        assert sum(cls['matched_by_region'].values()) == n
        for group, c in cls['sample_sets'].items():
            assert sum(c['by_reason'].values()) == c['evaluated_genotypes'] == n*f['sample_selection'][group]
    print(json.dumps({'unit_id':r['unit_id'],'source_samples':r['samples_in_source'],
                      'frequencies':f}, indent=2))
PY
```

Confirm source/representative counts agree with the existing PSAM and inspect
unknown-sex, female-Y, heterozygosity and partial-call exclusions. Star audit
counts stay separate. A run with no matched X PAR variants does not validate
real PAR counting; inspect `matched_by_region` and include a PAR block when ready.

Repeat the identical pilot command with a new trace filename and `-resume`;
require CACHED. Then use `--select_units all` **only with `xy-carriers.tsv`**, keeping
PSAM, parameters and work/output directories unchanged. This expands to X/Y only.
Do not reattach the 951 autosomal rows or change their output/work locations.

## Local validation and limits

95 targeted tests passed, followed by four additional exclusion tests: PAR start,
end and adjacent positions on both chromosomes; actual indexed X/Y extraction;
PSAM/sample-set rules; native and diploid-encoded haploidy; unknown/female-Y
exclusions without GT access; AN=0; raw/source output parity; separate stars;
unmatched Y-PAR rejection before VCF access; unchanged autosomal tables/receipts;
and local Nextflow explicit-policy enforcement, X/Y wiring, resume and durable
publication. Both TBI and CSI fixtures were exercised.

No NBDC jobs were launched. This X/Y revision still needs the operator's real
pilot validation. Final corrected unrelated AF <0.001 filtering remains separate.
