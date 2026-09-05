# Rust FastVEP picker validation: G2MH chr22

Date: 2026-09-04

## Implementation

The standalone `fastvep-picker` Rust binary reproduces the production subset of
`scripts/pick_fastvep_consequences.py`. It reads FastVEP VCF from a file or
standard input, applies the same versioned transcript-priority and consequence-
rank tables, and emits the same selected-transcript TSV schema. The Python
implementation remains the development oracle.

Source commit: `1671c5a76c369da50c64320d7dc3c719ac2ab95a`

Validated image:
`640838474376.dkr.ecr.us-east-1.amazonaws.com/rare-variant-pipeline-targeted@sha256:7d5b76a28e2427ca97af6ebec1d5e38aec63b93609419c0cc77c063e48e2917d`

## Full-chromosome parity

Input was the normalized all-observed G2MH chr22 VCF containing 1,148,541 allele
records. The Rust picker emitted 1,148,541 data rows. Both its materialized-input
and streamed outputs were byte-identical to the established Python-picked
production checkpoint:

```text
SHA-256  37ddb67c9c4908aebb7785b5e750c9acbecfa8e685991916d5cfc2f430d2c5d4
```

## Performance

Measurements used the same warm persistent 8-vCPU AWS host and FSx resources as
the split Python benchmark.

| Path | Wall time |
|---|---:|
| Python picker, materialized FastVEP VCF | 2:47.64 |
| Rust picker, materialized FastVEP VCF | 19.20 seconds |
| FastVEP streamed into Rust picker | 43.87--48.62 seconds |

The standalone picker improved by 8.7-fold. The streamed production path also
avoids the 3.1-GB raw FastVEP intermediate. The previous all-observed chr22
FastVEP-plus-Python-picker checkpoint took approximately 124 seconds, so the
measured combined annotation/picking stage improved by approximately 2.5--2.8-fold
on this workload. The final immutable image produced the exact expected checksum
in 48.62 seconds.

Small timing records are retained under:

```text
/fsx/loftee-parity/workflows/benchmarks/rust-fastvep-picker-052be33/
```

The temporary 3.1-GB raw VCF should be deleted after validation.
