# ABCD frequency, final rarity, QC and chromosome gathering validation

NBDC validation is complete for the pipeline implemented in PRs #18–#22.
The operator supplied the results below; Codex did not access NBDC or independently
inspect protected data/receipts. No further NBDC reruns are needed to finalize
this validation. Saved outputs and pinned runbooks remain the execution record.

## Validated stages

| Stage | Operator-reported result |
|---|---|
| Batched extraction and corrected frequencies | 951 autosomal plus 25 X/Y receipts passed; X/Y pilot resume cached; no zero-AN X/Y alleles reported |
| Final unrelated rarity | 976/976 receipts reconciled; 493,527 → 481,905 alleles; 1,356,838 → 941,371 carrier annotations; all 11,622 failed alleles had unrelated AF ≥0.001 |
| Post-rarity QC | 976/976 receipts reconciled; 941,371 → 454,232 carrier annotations |
| Chromosome gathering | 976 blocks → 24 chromosome outputs, reconciled; all 8,877 source samples retained |

Passing QC annotations: 107,704 sequence HC LoF, 329,667 missense and 16,861
spanning-deletion records (sum 454,232). Stars remain separate from sequence
burdens. These are annotation associations, not necessarily distinct allele/sample
counts or genome-wide unique-carrier sample counts. No additional unique-count
results are inferred from this report.

Final chromosome products on NBDC:

```text
/home/ood-guevara-james/abcd-fastvep-smoke/post-rarity-qc-gather-nextflow/output/post-rarity-qc-gather/<chromosome>/all/
```

Each chromosome directory contains gathered carriers, the full sample universe,
sample, sample/gene and gene summaries, distinct allele/sample summaries and a
reconciliation receipt. Individual-level products remain on NBDC; this document
contains aggregate evidence only.

## Implementation provenance and runbooks

| PR | Implementation/documentation head before finalization | Runbook |
|---|---|---|
| [#18](https://github.com/james-guevara/rare-variant-pipeline/pull/18) | `c822ddb422660a097962d927b793b38b0730ca7e` | [Batching and preliminary screen](../operations/abcd-fastvep-smoke/CARRIER_BATCHING.md) |
| [#19](https://github.com/james-guevara/rare-variant-pipeline/pull/19) | `e26e0d48e1ffbbab81da2084df180cd9809636c8` | [Corrected frequencies](../operations/abcd-fastvep-smoke/CARRIER_FREQUENCIES.md), [X/Y policy](../operations/abcd-fastvep-smoke/SEX_CHROMOSOME_FREQUENCIES.md) |
| [#20](https://github.com/james-guevara/rare-variant-pipeline/pull/20) | `5bee30e148ede0d88af8e11539e6e78c12aea29b` | [Final rarity](../operations/abcd-fastvep-smoke/FINAL_RARITY.md) |
| [#21](https://github.com/james-guevara/rare-variant-pipeline/pull/21) | `8470b610f81e3f6d484dbe087902a4894d7d1260` | [Post-rarity QC](../operations/abcd-fastvep-smoke/POST_RARITY_QC.md) |
| [#22](https://github.com/james-guevara/rare-variant-pipeline/pull/22) | `ee68dd5b258babb160d0e589e03271c6ff9baa35` | [Chromosome gathering](../operations/abcd-fastvep-smoke/GATHER_POST_RARITY_QC.md) |

Finalization adds documentation only to the final PR. Merge commits preserve
these implementation histories and their code pins. Earlier runbook statements
about deferred functionality or pending NBDC validation describe the stage at
its original implementation time; this record supersedes their validation
status. The later frequency, rarity, QC and gather stages resolve those earlier
deferrals without changing historical commands/pins.

## Current policy and scope

The preliminary filter uses uncorrected source AC/AN <0.005, gnomAD joint
POPMAX <0.001 or missing, and the established region exclusions. Extraction uses
10 kb pysam batching. Corrected cohort frequencies use one representative per
participant; unrelated frequencies use representative AND unrelated. Eligible
reference genotypes contribute to AN. Frequency counting adds no genotype QC;
the explicit GRCh38 X-only PAR policy controls X/Y ploidy.

Final rarity requires saved unrelated AF <0.001 with positive AN. Subsequent QC
uses site PASS, GQ ≥20, DP ≥10, heterozygous AB 0.25–0.75 and homozygous/haploid
ALT AB ≥0.90, with the agreed sex/PAR exclusions. FT is not a filter. Gathering
validates QC policy/hashes/PSAM identity, preserves raw GT/dosage, uses effective
dosage, rejects cross-block overlaps, and recomputes sample unions. Zero-carrier
samples, tiers and untiered HC are retained.

This closes implementation validation through chromosome gathering. It does
not change the scientific policy, resolve the deferred biological interpretation
of star records or low-AB LoF losses, or add a genome-wide analysis stage.
