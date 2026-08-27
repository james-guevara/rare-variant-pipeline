# Sex-chromosome handling

## Current status

The G2MH chrX and chrY sharded VCZ stores are complete, but the production
targeted workflow does not yet run them. Chromosome-specific FastVEP, LOFTEE,
dbNSFP, target BED, problematic-region, and manifest resources still need to be
staged and validated.

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
