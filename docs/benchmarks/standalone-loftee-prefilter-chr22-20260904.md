# Standalone LOFTEE candidate-prefilter benchmark: G2MH chr22

Date: 2026-09-04

## Change

Standalone LOFTEE now inspects consequence and biotype before transcript lookup
or resource-dependent context construction. By default it emits only classified
protein-coding stop-gained, frameshift, splice-donor, and splice-acceptor rows.
`--include-non-lof` retains the legacy full-table representation for debugging.

## Inputs

- FastVEP/picker table: G2MH chr22, 1,148,541 allele rows
- Original LOFTEE runtime: 129 seconds
- Optimized source commit: `61c21e0b7f63`
- Image digest:
  `sha256:a43ef2306bd3961ebf9c328f58c946fae4e857de2584843c5ed9eaa76b3cfd07`

## Results

| Measurement | Original | Optimized |
|---|---:|---:|
| Rows scanned | 1,148,541 | 1,148,541 |
| Rows classified/emitted | 399 | 399 |
| Warm runtime | 129 s | 7.33 s |
| Output size | 143 MB | 108 KB |

The warm runtime improved by approximately 17.6-fold. The initial optimized run,
including the one-time container pull, took 19.82 seconds.

## Parity validation

The optimized output was compared byte-for-byte with the header plus all `HC` and
`LC` rows selected from the original full-table output. The files matched exactly:

```text
rows including header: 400
sha256: 66d1ac09daa6dab066db3baf93965ff0f8b3f786eb4edab14223682113cec74c
comparison: EXACT_MATCH
```

This establishes parity for every retained column, including `LoF`, `LoF_filter`,
`LoF_flags`, and `LoF_info`.
