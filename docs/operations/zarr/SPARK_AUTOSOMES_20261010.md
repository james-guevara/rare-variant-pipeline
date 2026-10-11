# SPARK WGS 2026_08 autosomal rare workflow

Status: batch preparation/launch; not a completed release.

Source code: james-guevara/rare-variant-pipeline, feat/zarr-policy-backend,
commit a75e97f05c30257cab1cbc84d959357378a89377. Nextflow 26.04.6.
Inputs: 1,125 nonempty published Zarr blocks across chr1–22. X/Y excluded here.
The prior single-block and legacy chr22 experiments remain separate; this run
includes the pilot block again so all autosomal outputs share one policy/root.

Layout:
- inputs/: complete block manifest, chromosome order, PSAM, checksum resource locks.
- code/: frozen snapshot of the tested repository revision.
- runs/chrN/: stage input manifests, Nextflow launch/cache/work and execution reports.
- results/sites, annotation, candidates, pre-carrier, filtered-carriers,
  final-rarity, post-rarity-qc: stage outputs keyed by unique source block ID.
- results/post-rarity-qc/<unit_id>/: final block carrier table, sample and gene
  burden tables, distinct-allele counts, QC audit, sample list and receipt.
- inventory/: expected-versus-completed block inventories with final receipt and
  published-output size/hash checks. A complete block inventory does NOT mean
  chromosome/genome-wide burden gathering has been performed.
- logs/: resource-lock job and chromosome-controller logs.

Scientific policy: unfiltered dbNSFP parquet_expanded, no MANE inclusion filter;
ABCD picked annotation/tiering, untiered HC LoF retained, stars kept separate;
preliminary all-source genotype-derived AF <0.005, gnomAD joint POPMAX <0.001
or missing, existing problematic-region rules; corrected unrelated AF <0.001
with positive AN; explicit representative/family-proxy sample selection;
post-rarity PASS-or-dot site filter, GQ >=20, DP >=10, heterozygous AB 0.25–0.75,
homozygous ALT AB >=0.90. Multiallelic candidate projection supported. Source
GT/AD retained in carriers. Source Zarr stores are read-only.

Run order: sites -> annotation -> candidates -> Zarr preliminary frequencies /
pre-carrier filter -> Zarr carriers / corrected frequencies -> final rarity ->
post-rarity QC -> published block inventory. Three chromosome controllers,
ten tasks/controller. Each stage is a native Nextflow invocation. No global
barrier across all chromosomes. Controller allocation: 1 CPU, 2 GB, 24 hours.
Worker allocations use the tested stage settings (carriers 1 CPU/16 GB,
pre-carrier 2 CPU/8 GB, annotation 2 CPU/8 GB, candidates 2 CPU/4 GB,
final rarity and QC 1 CPU/1 GB). No upstream conversion or PLINK reruns.

Recovery: identify the failed stage and its UUID from that chromosome's
launch-<stage>/.nextflow.log. Confirm no surviving controller/orphan submission,
then resume that exact stage UUID and continue later stages. Do not blindly
restart the entire chromosome driver. Preserve previous execution reports.

Inventory checks never delete source data or intermediates. After completion,
review inventory/all-autosomes.json and address missing/failed blocks before
calling these block outputs complete. Chromosome/genome burden gathering must
reconcile duplicate allele/sample associations and unique samples separately.

Submitted resource-lock job 54829415, chromosome array 54829416, and final
inventory job 54829417. Production root:
`/expanse/projects/sebat1/j3guevar/SPARK-WGS-2026_08/rare/autosomes-20261010/`.
