# ABCD FastVEP 115 + Rust picker smoke test

No rare-pipeline changes or branch merge are required. [`benchmark.py`](benchmark.py) is an operational command wrapper, not another annotation implementation. Run it on NBDC using the existing generic binaries. It performs no filtering, VEP execution, LOFTEE execution, or genotype processing.

## Exact image and resources

Validated image from the [Rust picker benchmark](https://github.com/james-guevara/rare-variant-pipeline/blob/b686ef796ed6706bcf0ab15ff3731021441f1cd1/docs/benchmarks/rust-fastvep-picker-chr22-20260904.md):

```text
640838474376.dkr.ecr.us-east-1.amazonaws.com/rare-variant-pipeline-targeted@sha256:7d5b76a28e2427ca97af6ebec1d5e38aec63b93609419c0cc77c063e48e2917d
```

FastVEP source: `james-guevara/fastVEP@cb8113d7bab2db42cb06bb2b2a40c57b60ea2561`.
Rust picker: version `0.1.0`, documented source `1671c5a76c369da50c64320d7dc3c719ac2ab95a` in `rare-variant-pipeline`.

Stage the matching Ensembl 115 resource bundle with these six files in one directory:

```text
Homo_sapiens.GRCh38.115.chr22.gff3
Homo_sapiens.GRCh38.115.chr22.gff3.fastvep.cache
Homo_sapiens.GRCh38.dna.primary_assembly.fa
Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
vep115.transcript-priority.tsv
vep115.consequence-ranks.tsv
```

These exact six consolidated resources are in `s3://sebat-genomics-work/resources/rare-variant-pipeline/v1/targeted-annotation/ensembl-115/`. Their expected hashes are supplied in [SHA256SUMS.chr22](SHA256SUMS.chr22), checked against the existing GitHub runtime manifest, the S3 checksum index, and S3 object SHA-256 metadata on 2026-09-29. The restored Expanse copy is `/expanse/projects/sebat1/resources/rare-variant-pipeline/releases/v1/targeted-annotation/ensembl-115/`; use this versioned path. Historical chromosome-specific filenames are retained only for compatibility. See the [full resource inventory](../../resources/README.md) and the [existing consolidation documentation](https://github.com/james-guevara/integrated_genomics_pipeline/blob/main/docs/operations/rare-variant-resource-locations.md).

The ECR registry requires authorized access. If NBDC cannot pull it, build/export a SIF from this exact digest on an authorized machine and transfer that SIF with its SHA-256. Keep its build provenance. A SIF checksum is not the OCI digest. Do not assume an older `41de024` SIF has the validated Rust picker. No AWS/Expanse resource paths are assumed here.

On an authorized machine, the pull command is:

```bash
apptainer pull targeted-rust-picker.sif \
  docker://640838474376.dkr.ecr.us-east-1.amazonaws.com/rare-variant-pipeline-targeted@sha256:7d5b76a28e2427ca97af6ebec1d5e38aec63b93609419c0cc77c063e48e2917d
sha256sum targeted-rust-picker.sif > targeted-rust-picker.sif.sha256
```

## Already staged on Expanse

The repaired Expanse registry now supplies the exact resources and pinned smoke-test SIF:

```bash
RES=/expanse/projects/sebat1/resources/rare-variant-pipeline/releases/v1/targeted-annotation/ensembl-115
SIF=/expanse/projects/sebat1/resources/rare-variant-pipeline/containers/v1/targeted-7d5b76a28e2427ca97af6ebec1d5e38aec63b93609419c0cc77c063e48e2917d.sif
```

The SIF has a portable `.sha256` sidecar. Transfer that image, its sidecar, and the six annotation resources to NBDC, preserving timestamps where possible. Verify the resource and SIF checksums at the destination and enforce the cache-newer-than-GFF rule below before binding the resources read-only. See the [repair report](../../resources/repair-20260929/README.md) for validation and exact image hashes. Codex has not run the ABCD benchmark on NBDC.

## Download the six verified resources

On a machine with authorized S3 access, use an empty staging directory and the provided checksum file. Transfer the resulting files and checksum file to NBDC if NBDC has no AWS access. No presigned URLs or credentials need to be shared with ChatGPT.

```bash
set -euo pipefail
export RES=/your/staging/path/ensembl-115
export HANDOFF=/your/path/abcd-fastvep-smoke
mkdir -p "$RES"
S3_ANNOTATION=s3://sebat-genomics-work/resources/rare-variant-pipeline/v1/targeted-annotation/ensembl-115
for name in \
  Homo_sapiens.GRCh38.115.chr22.gff3 \
  Homo_sapiens.GRCh38.115.chr22.gff3.fastvep.cache \
  Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai \
  vep115.transcript-priority.tsv \
  vep115.consequence-ranks.tsv; do
  aws s3 cp "$S3_ANNOTATION/$name" "$RES/$name" --only-show-errors
done
(cd "$RES" && sha256sum -c "$HANDOFF/SHA256SUMS.chr22")
```

FastVEP decides cache freshness using **strict filesystem modification-time ordering** (`cache_mtime > gff_mtime`). S3 copies can give the pair equal or reversed timestamps. After verifying downloaded bytes, set the cache timestamp newer than its matching GFF3; this changes metadata only, not the pinned SHA-256. Then bind resources read-only as shown below. Otherwise FastVEP may rebuild a valid cache or waste time reparsing GFF3.

```bash
python3 - "$RES" <<'PY_CACHE_TIME'
import os, pathlib, sys
root = pathlib.Path(sys.argv[1])
gff = root / 'Homo_sapiens.GRCh38.115.chr22.gff3'
cache = pathlib.Path(str(gff) + '.fastvep.cache')
ns = max(cache.stat().st_mtime_ns, gff.stat().st_mtime_ns + 1_000_000_000)
os.utime(cache, ns=(ns, ns))
PY_CACHE_TIME
```

These six files total approximately 3.20 GB. Run the checksum check again after transfer to NBDC. Downloading/hashing resources will warm filesystem caches; record whether the benchmark is a first or repeated invocation.

## NBDC presence checks

Set these to your actual NBDC locations. Copy `benchmark.py` there too. `apptainer` can be replaced by `singularity` if that is the installed runtime.

```bash
export SIF=/your/nbdc/path/targeted-rust-picker.sif
export RES=/your/nbdc/path/ensembl-115
export INPUT=/your/nbdc/path/chr22_block12.sites.vcf.gz
export HANDOFF=/your/nbdc/path/abcd-fastvep-smoke-handoff
export OUT=/your/nbdc/path/block12-fastvep-smoke-01

command -v apptainer
test -s "$SIF" && test -s "$INPUT" && test -s "$HANDOFF/benchmark.py"
for name in \
  Homo_sapiens.GRCh38.115.chr22.gff3 \
  Homo_sapiens.GRCh38.115.chr22.gff3.fastvep.cache \
  Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai \
  vep115.transcript-priority.tsv \
  vep115.consequence-ranks.tsv; do
  if test -s "$RES/$name"; then echo "FOUND $name"; else echo "MISSING $name"; fi
done
apptainer inspect "$SIF"
apptainer exec --cleanenv "$SIF" sh -c \
  'command -v fastvep; command -v fastvep-picker; command -v bcftools; command -v python'
apptainer exec --cleanenv "$SIF" fastvep annotate --help
apptainer exec --cleanenv "$SIF" fastvep-picker --help
```

If you need to locate existing copies, restrict the search to your resource/container storage root:

```bash
find /your/nbdc/resource-storage-root -type f \( \
  -name '*.sif' -o -name 'Homo_sapiens.GRCh38.115.chr22.gff3*' -o \
  -name 'vep115.transcript-priority.tsv' -o -name 'vep115.consequence-ranks.tsv' -o \
  -name 'Homo_sapiens.GRCh38.dna.primary_assembly.fa*' \) -print
```

## Run the benchmark

Use a fresh output directory. Run within your normal NBDC compute allocation, not necessarily a login node. All paths above should be absolute.

```bash
mkdir -p "$OUT"
SIF_SHA=$(sha256sum "$SIF" | cut -d ' ' -f 1)
apptainer exec --cleanenv \
  --bind "$RES:/resources:ro" \
  --bind "$INPUT:/input.sites.vcf.gz:ro" \
  --bind "$HANDOFF:/smoke:ro" \
  --bind "$OUT:/output" \
  "$SIF" python /smoke/benchmark.py \
    --input /input.sites.vcf.gz --resources /resources \
    --chromosome chr22 --outdir /output --sif-sha256 "$SIF_SHA"
```

Success leaves only `picked.tsv` and `receipt.json`. Share the receipt/aggregate console summary, not picked rows. Expected input and picked counts for block12 are **48,001**; the wrapper checks their equality. It captures subprocess errors privately on failure and deletes partial picked output.

The wrapper uses a temporary uncompressed VCF because the pinned FastVEP reader consumes plain VCF. It removes only the old `INFO/CSQ`: pinned FastVEP appends fresh CSQ while the picker selects the first CSQ field, so leaving VEP112 CSQ would select the wrong annotation. It preserves all records and the durable source. It aligns `chr22`/`22` aliases only when necessary to match both GFF3 and FASTA, recording the mapping; incompatible contigs fail.

The scientific command is the existing production pipeline:

```bash
fastvep annotate --input temporary-clean.vcf \
  --gff3 "$GFF3" --fasta "$FASTA" --transcript-cache "$CACHE" \
  --hgvs --symbol --canonical --output-format vcf --output - |
fastvep-picker --fastvep - --transcript-priority "$PRIORITY" \
  --consequence-ranks "$RANKS" --output picked.tsv
```

Reported wall time includes FastVEP/picker startup and resource loading, excludes temporary-input preparation, container startup and subsequent hashing. CPU is summed child user+system time. RSS is the largest individual child high-water mark on Linux, not the sum of simultaneously resident processes. Hashing occurs afterward to avoid warming resources before the measurement. Record whether this was the first or a repeated run; neither command flushes filesystem caches.

No real ABCD/NBDC benchmark has been performed by Codex. Inspect runtime/counts first; standalone LOFTEE is a separate next experiment.

## Local validation

The command wrapper was checked with synthetic data using real bcftools and the repository Rust picker, with FastVEP stubbed. Checks covered old CSQ removal, contig alias conversion, two-record/two-row parity, receipt generation, and temporary-file cleanup. Python and shell snippets passed syntax checks. This is not an ABCD annotation or container validation.

Updated production extraction: [10 kb pysam batching and preliminary source AF <0.005](CARRIER_BATCHING.md). Corrected frequencies remain deferred.
