# ABCD chr22 post-QC gathering — operator-reported NBDC validation

Implementation: `282a0862cb479deeb284e54b0469f5e1c7e45450` (PR #15).
The operator reported the following completed validation in this session:

- Block12 passed.
- Repeated block12 execution with `-resume` was cached.
- Blocks 0/12/19 passed.
- All 23 chr22 blocks passed; all receipts reconciled.
- Zero cross-block duplicate alleles or carrier associations.
- All 8,877 source samples retained.

| Allele class | Distinct alleles | Distinct allele/sample records | Samples with carriers |
|---|---:|---:|---:|
| Sequence | 7,023 | 9,955 | 5,839 |
| Spanning deletion (`*`) | 231 | 397 | 247 |

Passing annotation associations:

| Candidate class | Associations |
|---|---:|
| Missense | 7,743 |
| Sequence HC LoF | 2,212 |
| Spanning deletion HC LoF | 397 |
| Total annotation associations (bookkeeping only) | 10,352 |

The sequence association sum is 7,743 + 2,212 = 9,955. Stars remain separate.
Do not add the two unique-sample counts: the same sample can appear in both.
The total above is not a combined biological burden.

Provenance limitation: this is a record of the operator's aggregate report, not
an independent inspection of NBDC receipts, protected data, traces or files by
Codex. No original NBDC receipt hashes or job IDs were supplied. The prior local
synthetic tests are documented in the [gather runbook](../operations/abcd-fastvep-smoke/GATHER_POST_QC.md).

This establishes the reported chr22 execution gate. It does not resolve the
separately deferred scientific audit, deletion-event interpretation, callable
coverage, or sex-chromosome policy. Remaining autosomes use the existing policy.
