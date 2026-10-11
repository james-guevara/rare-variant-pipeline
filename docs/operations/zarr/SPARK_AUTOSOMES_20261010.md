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

## chr7 / chr11 recovery — 2026-10-11

- chr7 failed before QC while Python imported gzip/zlib: `libz.so.1` returned an I/O error. Controller **54851784** resumed session `6b1b99a2-ee2c-4de7-841b-557627f7b988`, using a SHA-256-verified node-local copy of the identical SIF. It reused 28 tasks and completed 36. `inventory/chr7.json` passed **64/64**. The shared-filesystem failure itself was not reproducible enough to assign a lower-level cause.
- chr11 failed reproducibly in FastVEP's frameshift HGVSp code. Incomplete-CDS coordinates included two virtual bases but the edited CDS+UTR sequence omitted them. Baseline `cb8113d7` reproduced the panic on a single private record and on a synthetic regression; the shared phase-aware helper in FastVEP **b1b79335af455a4f300c9bd6fa68274f2bfb94ae** fixes both CLI and library paths. See https://github.com/james-guevara/fastVEP/pull/1 (draft, not merged).
- Validation: 40 library tests; CLI phase 0/1/2 regression (baseline fails, patched passes); six annotation-adapter tests. The complete failing block passed with **437,566 input records and picked output rows**, taking 31.16 seconds. No variant suppression or filtering-policy changes.
- Pinned rebuilt SIF: `/expanse/projects/sebat1/j3guevar/operations/fastvep-fix-20261011/fastvep-b1b7933-r2.sif`; SHA-256 `74a265bc86774aa1a7d6d1343c0a19a549427ff18e4691664fdbd69c569a707a`. Build provenance inside `/opt/rvp/fastvep-build.json` is checked against the executable and included in receipts.
- chr11 controller **54851842** uses code snapshot `code-fixes-95a22b9`, explicit annotation session `d41121e9-9b68-4912-9b35-7511fee0fd55`, and `inputs/resource-locks/annotation.fastvep-b1b7933.json`. It reuses all sites exports, reruns annotation with the corrected pin, then continues downstream. All containers are localized with checksums before execution. The chromosome is not complete until its final inventory passes.
- Recovery scripts, configs, and pins: `/expanse/projects/sebat1/j3guevar/operations/rare-fixes-20261011/`. Diagnostic/private evidence: sibling `fastvep-fix-20261011/`. Earlier failed packaging and config attempts remain as diagnostic logs; only the tested `r2` image and controller above are active.

The FastVEP fix changes frameshift HGVSp construction for incomplete CDS starts. Previously completed chromosomes retain their original annotation pin; they have not been reannotated by this targeted recovery. Phase-zero construction is unchanged.
