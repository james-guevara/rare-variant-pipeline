# Sex-chromosome handling

## Current status

The G2MH chrX and chrY sharded VCZ stores are complete. Ensembl 115 FastVEP,
LOFTEE, and transcript-priority resources are staged for both chromosomes, and
the chrX LoF-only branch has completed its first end-to-end validation. The
missense dbNSFP tables and combined target inputs still need to be prepared;
chrY analysis policy and postprocessing also require validation before a
production sex-chromosome run.

The VCZ genotype arrays have fixed ploidy width two. Effective ploidy is encoded
by `call_genotype_mask`, not by array shape or the rendered GT string. In G2MH,
haploid chrX non-PAR alternate calls render as `1/.`. Other callers or conversion
paths may render the same state as `1` or `./1`.

Carrier extraction therefore emits these mask-derived fields:

- `called_ploidy`
- `alt_dosage`
- `is_heterozygous`
- `is_homozygous_alt`
- `is_hemizygous_alt`

Genotype QC uses `called_ploidy` and `alt_dosage`. `GT` is retained only for
display and provenance.

## G2MH ploidy audit

Run:

```bash
python scripts/audit_sex_chromosome_ploidy.py \
  --config config/sex-chromosomes/g2mh.json \
  --x-zarr /path/to/chrX.sharded-v3.zarr \
  --y-zarr /path/to/chrY.sharded-v3.zarr \
  --autosome-zarr /path/to/chr22.sharded-v3.zarr \
  --output g2mh.inferred-karyotype.tsv
```

The 2026-08-27 G2MH audit classified 1,065 samples as:

- 549 `XX-like`
- 510 `XY-like`
- 6 `ambiguous`

Fourteen otherwise classified samples had discordant chrX/autosome depth. The
saved FSx audit is checksum-pinned at:

```text
ec8268ccdeec4e888a8c1d5dc9a2bbde728555fbab0c848cb708eccd90784721
```

This is an inferred karyotype QC result, not a declaration of biological sex.
It should be compared with authoritative cohort metadata when that becomes
available. Ambiguous and depth-discordant samples must remain visible rather
than being silently forced into a binary category.

For non-diagnostic triage of ambiguous or discordant samples, run:

```bash
python scripts/interpret_sex_chromosome_qc.py \
  --audit g2mh.inferred-karyotype.tsv \
  --config config/sex-chromosomes/g2mh.json \
  --output g2mh.sex-chromosome-qc.tsv
```

The added copy-number evidence patterns are screening flags, not clinical
karyotypes. Ambiguous samples receive `sex_chromosome_analysis_eligible=0` by
default. In particular, excess raw chrY depth is only a prompt for orthogonal
review until a uniquely mappable chrY mask has been validated.

Carrier rows are retained and annotated before burden collapse:

```bash
python scripts/annotate_sex_chromosome_carriers.py \
  --input chrX.carriers.parquet \
  --sample-qc g2mh.sex-chromosome-qc.tsv \
  --regions resources/grch38-sex-chromosome-regions.json \
  --output chrX.carriers.sex-qc.parquet
```

`burden_count_available` preserves ambiguous non-PAR counts for QC and
sensitivity analysis. `primary_analysis_eligible` and
`frequency_denominator_eligible` remain false for those rows. This separation
prevents uncertainty about karyotype from destroying observed burden data or
silently entering the primary test. Noncanonical chrY PAR records are retained
for provenance but have `burden_count_available=false` to prevent double counts.

The genotype-summary sidecar likewise retains whole-cohort `genotype_*` fields
and adds `primary_genotype_ac`, `primary_genotype_an`, and
`primary_genotype_af`. Sex-chromosome cohort-AF filtering explicitly selects the
`primary_genotype` prefix. Thus ambiguous carriers remain inspectable without
contributing alleles to the primary frequency numerator or denominator.

The chrX/autosome median-DP ratio is a strong independent signal in G2MH:
approximately 1.0 for XX-like samples and 0.5 for XY-like samples. Raw chrY depth
is inflated by repetitive/multicopy mapping and is diagnostic only; it is not a
classification threshold without a validated uniquely mappable chrY mask.

The chrY 11-13 Mb block has extensive off-target calls in XX-like samples and
must not be used for inference. The audit uses chrY 3-5 Mb, which cleanly
separates Y-supported and Y-absent G2MH samples.

## Proposed analysis policy

1. Treat GRCh38 chrX PAR1 and PAR2 as diploid and apply autosomal genotype QC.
2. In chrX non-PAR, accept diploid heterozygous/homozygous and haploid
   hemizygous calls using mask-derived state.
3. In chrY non-PAR, include only samples eligible under the cohort karyotype
   resource. Exclude XX-like and ambiguous samples from both carriers and AN.
4. Preserve combined and karyotype-stratified AC, AN, and AF. Cohort eligibility
   thresholds must use the explicitly configured measure.
5. Prefer the chrX representation for PAR variants. If a future input contains
   homologous chrY PAR records, exclude them to prevent double counting.
6. Fail preflight when a sex-chromosome run lacks a checksum-pinned sample
   karyotype/sex resource and GRCh38 PAR definition.

The current G2MH chrY primary-contig store contains no PAR1 or PAR2 variants, so
there is no duplicate PAR burden in the existing data. This is an observed
dataset property, not a portable assumption.

Reference-build PAR coordinates are stored separately in
`resources/grch38-sex-chromosome-regions.json`. Inference windows and thresholds
are cohort calibration in `config/sex-chromosomes/g2mh.json`; they are not global
defaults and must be recalibrated or validated for each additional cohort.

## First chrX LoF validation

The checksum-pinned G2MH chrX LoF-only run completed on AWS on 2026-08-28.
Targeted extraction covered 6,180 intervals and 40,340 observed alleles;
FastVEP took 9 seconds and standalone LOFTEE took 8 seconds. Of 94 tiered LoF
alleles, problematic-region filtering reduced 536 raw carrier rows to 109 and
genotype QC retained 67 rows across 38 alleles and 65 samples. Sixty-six rows
are primary eligible and one ambiguous-sample `ABCB7` carrier is retained for
sensitivity analysis only. Exact counts and canonical hashes are pinned in
`resources/g2mh-chrX-lof-regression.json`.

The corresponding chrY LoF-only validation targeted 254 intervals and 606
observed alleles. Eight tiered LoF alleles produced 11 raw carrier rows: two
XY-like primary rows, one ambiguous sensitivity row, and eight XX-like Y calls
retained for QC provenance only. Region filtering retained eight rows; genotype
QC retained one high-quality, haploid `KDM5D` LoF carrier in an XY-like sample.
No ambiguous chrY burden row survived genotype QC. Exact chrY counts and hashes
are pinned in `resources/g2mh-chrY-lof-regression.json`.
