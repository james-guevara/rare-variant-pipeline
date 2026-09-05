# FastVEP and picker split benchmark: G2MH chr22

Date: 2026-09-04

## Contract

- FastVEP fork: `james-guevara/fastVEP`
- FastVEP commit: `cb8113d7bab2db42cb06bb2b2a40c57b60ea2561`
- Container digest:
  `sha256:a43ef2306bd3961ebf9c328f58c946fae4e857de2584843c5ed9eaa76b3cfd07`
- Host: warm persistent AWS orchestration instance, 8 vCPUs
- Input: normalized all-observed G2MH chr22 sites VCF
- Input allele records: 1,148,541

The fork includes actual FastVEP source changes for indel `HGVS_OFFSET` parity
and terminal shifted-insertion consequence handling. VEP-compatible transcript
picking remains a separate downstream pipeline implementation.

## Results

| Stage | Wall time | Output size |
|---|---:|---:|
| FastVEP only | 52.45 seconds | 3.1 GB raw VCF |
| Picker, raw VCF on FSx | 2:47.64 | 136 MB TSV |
| Picker, raw VCF in memory | 2:46.63 | 136 MB TSV |

The memory-backed result demonstrates that the separated picker's runtime is
primarily Python parsing and ranking rather than FSx input latency. Materializing
FastVEP output also creates a temporary 3.1-GB file, so it is not the desired
production path.

Both picker measurements emitted 1,148,541 data rows and exactly matched the
existing streamed production checkpoint:

```text
SHA-256  37ddb67c9c4908aebb7785b5e750c9acbecfa8e685991916d5cfc2f430d2c5d4
```

## Decision

Keep FastVEP piped directly into the picker in production. If further annotation
speed is needed, optimize the picker first—preferably by implementing the
validated hierarchy in Rust or FastVEP itself—and require byte-identical chr22
regression output before adoption.

The benchmark artifacts are under:

```text
/fsx/loftee-parity/workflows/benchmarks/fastvep-picker-cb8113d/
```
