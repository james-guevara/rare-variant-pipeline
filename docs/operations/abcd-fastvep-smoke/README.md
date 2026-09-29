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

Use the validated cache and tables, not freshly generated approximations. The current [production runner](https://github.com/james-guevara/rare-variant-pipeline/blob/b686ef796ed6706bcf0ab15ff3731021441f1cd1/scripts/run_targeted_chromosome.sh) supplies these filenames. The documentation does not establish expected SHA-256 values for all six resources. Obtain a source checksum manifest when staging; this wrapper records actual NBDC checksums but cannot independently certify their provenance.

The ECR registry requires authorized access. If NBDC cannot pull it, build/export a SIF from this exact digest on an authorized machine and transfer that SIF with its SHA-256. Keep its build provenance. A SIF checksum is not the OCI digest. Do not assume an older `41de024` SIF has the validated Rust picker. No AWS/Expanse resource paths are assumed here.

On an authorized machine, the pull command is:

```bash
apptainer pull targeted-rust-picker.sif \
  docker://640838474376.dkr.ecr.us-east-1.amazonaws.com/rare-variant-pipeline-targeted@sha256:7d5b76a28e2427ca97af6ebec1d5e38aec63b93609419c0cc77c063e48e2917d
sha256sum targeted-rust-picker.sif > targeted-rust-picker.sif.sha256
```

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
