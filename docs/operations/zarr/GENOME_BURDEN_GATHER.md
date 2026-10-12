# Genome-wide post-rarity burden gathering

`gather_genome_burdens.nf` consumes completed chromosome products from
`gather_post_rarity_qc.nf`. It does not repeat annotation, genotype extraction,
frequency correction, or carrier QC. The input TSV has exactly two columns:
`chromosome` and `directory`, where each directory contains a passed receipt,
`carriers.qc.tsv.gz`, and `samples.tsv`.

Specify the expected chromosome set explicitly (the default is 1–22):

```bash
nextflow run gather_genome_burdens.nf \
  --chromosome_manifest chromosomes.tsv \
  --expected_chromosomes 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
  --gather_container /path/to/pinned-python.sif \
  --gather_publish_mode link --outdir /path/to/results
```

Use a site config to enable Slurm/Singularity, pin/localize the container and
set TMPDIR to job-local storage. Default task resources are 1 CPU / 4 GB / 2 h.
SQLite bounds memory and holds disposable intermediate tables in TMPDIR.

Checks require exactly one full, passed gather for every expected chromosome,
identical sample order, PSAM provenance, QC policy and strata. Carrier/sample
checksums must match receipts. Carrier rows must belong to the declared
chromosome, sample universe and source units. Duplicate annotation associations
and conflicting effective dosages fail closed. Per-chromosome carrier counts,
stratum counts and dosages reconcile before any final receipt is issued.

`genome/` contains sample, sample/gene and gene burdens, distinct allele/sample
summaries, a dense zero-inclusive sample universe, stratum summaries, chromosome
source manifest and a checksum-bearing reconciliation receipt. Unique samples
are recomputed as unions, not added across chromosomes. Sequence and spanning
-deletion records remain separate, as do tiers, untiered HC and candidate types.
Annotation associations and distinct allele/sample counts are separate outputs.
Sample/gene tables are sparse: missing cells are zero.

Chromosome carrier files remain the durable source; the genome stage does not
make a redundant concatenated carrier file. Source receipts and file hashes are
recorded in the genome receipt. Publication mode defaults to copy; link is
available for shared filesystems. Nextflow caches the declared chromosome inputs,
and the wrapper/resume path is covered by synthetic tests.
