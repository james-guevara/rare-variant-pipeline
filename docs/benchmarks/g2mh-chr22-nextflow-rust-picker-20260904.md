# G2MH chr22 Nextflow regression with Rust picker

Date: 2026-09-04

## Execution contract

- Integrated workflow commit: `1dee29159b6b543cb8f170a3c2c87da4470bc4b8`
- Rare-variant component: `ca95e891aaee35e02b814a1a4632e71410bba868`
- FastVEP: fork commit `cb8113d7bab2db42cb06bb2b2a40c57b60ea2561`
- Container digest:
  `sha256:7d5b76a28e2427ca97af6ebec1d5e38aec63b93609419c0cc77c063e48e2917d`
- Nextflow: 26.04.6
- Execution: AWS Batch with FSx for Lustre
- Scope: all 996,414 chr22 records and 1,148,541 observed ALT alleles

The run used the operator-facing command:

```bash
nextflow run rare_variant.nf \
  -profile aws_batch \
  -params-file config/runs/g2mh.chr22-test.rare.json \
  -resume
```

## Timing

| Task or stage | Time |
|---|---:|
| Nextflow preparation task | 18.3 seconds duration; 153 ms realtime |
| AWS runtime preflight | 19.1 seconds duration; 1.1 seconds realtime |
| chr22 scientific task | 2:10 duration; 1:47 realtime |
| Cohort gather | 19.6 seconds duration; 86 ms realtime |
| Zarr all-observed extraction | 51 seconds |
| Sites VCF emission | 4 seconds |
| `bcftools norm` | 1 second |
| FastVEP streamed into Rust picker | 24 seconds |
| Standalone LOFTEE | 6 seconds |

Batch task duration includes scheduling and staging; stage timings come from the
chromosome worker log.

## Regression result

The run produced 53 qualifying LoF alleles: 10 `lof_t1` and 43 `lof_t2`.
Genotype recovery found 1,341 raw carrier rows across 1,055 samples. Final
eligibility retained 38 carrier rows across 19 alleles and 26 samples.

The formal row/schema comparison passed for the tiered variants and carrier
Parquets. The FastVEP-picked table, tiered LoF table, region-filtered table,
genotype-QC table, population-AF tables, cohort-AF table, final eligibility
table, and four per-sample burden outputs were byte-identical to the prior
validated all-observed chr22 run.

The FastVEP-picked SHA-256 was:

```text
37ddb67c9c4908aebb7785b5e750c9acbecfa8e685991916d5cfc2f430d2c5d4
```

Standalone LOFTEE intentionally differs in storage representation: it now
contains only the 399 actual LoF candidates and is approximately 108 KB instead
of retaining blank annotations for more than one million non-LoF rows.

## Outputs

- FSx run root:
  `/fsx/rare-variant-production/g2mh/annotation-v2-test/chr22`
- Published results:
  `s3://sebat-genomics-work/results/integrated-genomics/g2mh-rare-production-v2-test/`
- Regression fixture:
  `resources/g2mh-chr22-all-observed-lof-regression.json`
