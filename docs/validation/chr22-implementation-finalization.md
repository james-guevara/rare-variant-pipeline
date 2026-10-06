# Chr22 implementation finalization

The implementation stack is finalized under the existing scientific policy.
PRs #10 (scored candidates), #11 (pre-carrier filtering), #12 (filtered carriers),
#13 (post-extraction QC), #15 (chromosome gathering), and #16 (generic autosome
readiness) were merged into `main` in dependency order using merge commits.
Resource consolidation/documentation PR #14 integrates the same runtime tree;
its overlapping resource README was reconciled to retain both access/download
instructions and the autosome guidance.

## Preserved validation and commits

| Stage | Original implementation commit |
|---|---|
| Candidate scoring | `cddc9e8a599b99952d919694889c147b73c5b1e0` |
| Pre-carrier filtering | `bca3ca91b8a2a00eef0db48d4209150f32424c08` |
| Filtered carriers | `f4adc07c08a116aef15892bd53f331ea4231cd89` |
| Post-extraction QC | `1e1731960ca750f396d02bac65899050fafb3273` |
| Gathering implementation | `282a0862cb479deeb284e54b0469f5e1c7e45450` |
| Autosome readiness | `4eedac8cb53b5b71eb4222b26ae65483f00207bc` |

The runtime tree (scripts, modules, libraries, tests, entrypoints and configs)
is identical to the reviewed autosome-readiness head. Documentation reconciliation
does not require another real chr22 run.

Before merging, the combined local regression run passed **131 tests**:

```bash
uv run --with pytest --with duckdb==1.5.5 --with pysam python -m pytest -q \
  tests/test_scored_candidates.py tests/test_pre_carrier_filter.py \
  tests/test_exact_carriers.py tests/test_filtered_carriers.py \
  tests/test_post_carrier_qc.py tests/test_gather_post_qc.py
```

This includes exact-allele and star handling, original HC-only compatibility,
functional scoring/AF/region boundaries, QC behavior, zero samples, association
versus distinct counts, and actual local Nextflow execution/resume. Resource
manifest hashes and documentation links were checked after reconciliation.

Real NBDC validation is [operator-reported](abcd-chr22-post-qc-gather-nbdc.md):
23 gathered blocks, 8,877 samples, zero cross-block duplicates, 9,955 sequence
allele/sample records and 397 separately retained star records. Codex did not
run NBDC jobs or independently inspect protected data there.

## Handoff

- Use a pinned final `main` merge commit for subsequent code deployments. The
  merge SHA is recorded in GitHub PR #14 and the session handoff, rather than
  embedding a self-referential commit hash in this file.
- Keep existing publications, `.nextflow` cache metadata, launch directories,
  work directories and resource locks. A new checkout alone does not require
  rerunning completed stages. Do not include chr22 in new execution manifests.
- Follow the [generic autosome readiness guide](../operations/abcd-fastvep-smoke/AUTOSOME_READINESS.md)
  for input contracts and the [resource guide](../resources/README.md) for
  canonical Expanse/S3 and the shared ddp195 copy.
- The operator handles actual NBDC paths, populated manifests, staging and jobs;
  repository work handles generic code, contracts, inventories and tests.
- Remaining autosomes proceed under the unchanged policy. LoF-indel low-AB losses,
  biological interpretation of spanning-deletion records, broader scientific
  auditing and sex-chromosome policy remain separate work. Finalization does not
  claim those questions are resolved.

No NBDC jobs, scientific resource revalidation, or real chr22 reruns were performed
during finalization. Existing branch names are retained for historical runbook links.
