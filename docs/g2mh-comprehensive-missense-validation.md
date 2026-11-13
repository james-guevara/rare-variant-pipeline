# Comprehensive G2MH missense validation

The corrected all-observed G2MH workflow was validated across chromosomes 1–22,
X, and Y against the previous rare-variant pipeline and by independently
recomputing every missense score flag and tier.

## Tier definitions

The workflow uses percentile-normalized dbNSFP `*_rankscore` columns and retains
the maximum available value rather than restricting score extraction to MANE:

| Predictor | Passing threshold |
|---|---:|
| `ClinPred_rankscore` | >= 0.4298 |
| `AlphaMissense_rankscore` | >= 0.9603 |
| `popEVE_converted_rankscore` | >= 0.9209 |
| `MPC_rankscore` | >= 0.8947 |

The direction is intentionally:

- `miss_t1`: four of four predictors pass;
- `miss_t2`: three of four pass;
- `miss_t3`: two of four pass;
- `miss_t4`: one of four passes.

Missing predictor values do not pass their threshold; they do not cause an allele
to be discarded when another predictor supports a valid tier.

## Selected-allele audit

| Tier | Selected alleles | MPC present | MPC passes |
|---|---:|---:|---:|
| `miss_t1` | 162 | 162 | 162 |
| `miss_t2` | 908 | 853 | 706 |
| `miss_t3` | 5,253 | 4,410 | 3,141 |
| `miss_t4` | 43,382 | 35,314 | 4,776 |
| **Total** | **49,705** | **40,739** | **8,785** |

Independent recomputation found zero `miss_n_flag` errors and zero tier-label
errors. Overall selected-allele score coverage was 48,144 for ClinPred, 49,328
for AlphaMissense, 43,355 for popEVE, and 40,739 (82.0%) for MPC.

FastVEP transcript picking is used to require a selected `missense_variant` in a
candidate gene. MANE participates in transcript priority but is not an
inclusion/exclusion filter. The dbNSFP maximum score is retained even when the
score-producing transcript is not MANE, avoiding the known loss of MPC coverage
from MANE-strict extraction.

## Legacy comparison

The corrected workflow restores nine multiallelic autosomal `miss_t4` carrier
rows that were absent from an earlier targeted extraction. Those nine rows pass
candidate selection, region filtering, genotype QC, population AF filtering, and
cohort AF filtering. With that normalized multiallelic matching fix, autosomal
missense carrier rows match the previous pipeline exactly.

Primary-analysis carrier totals are:

| Tier | Corrected G2MH carrier rows |
|---|---:|
| `miss_t1` | 113 |
| `miss_t2` | 611 |
| `miss_t3` | 3,693 |
| `miss_t4` | 41,148 |

The only differences from the previous totals (`113`, `612`, `3,695`, and
`41,149`) are four chrX carriers in samples with ambiguous sex-chromosome QC:
one `miss_t2`, two `miss_t3`, and one `miss_t4`. They are excluded only from the
primary analysis and are retained in the sex-chromosome sensitivity outputs.

## Conclusion

The comprehensive G2MH missense branch is validated without another annotation
run. Its score thresholds, tier direction, transcript-selection behavior,
multiallelic carrier recovery, and autosomal legacy parity are internally
consistent. Sex-chromosome primary-versus-sensitivity handling explains all
remaining differences.
