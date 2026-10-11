# Zarr inputs with the current rare-variant policy

This branch is stacked on `feat/post-rarity-qc` at
`8470b610f81e3f6d484dbe087902a4894d7d1260` (PR #21). It does not change
ABCD thresholds or replace the VCF workflow. It adds source adapters, keeping
candidate selection, preliminary filtering, final rarity, and post-rarity QC shared.

## Stage contract

1. `zarr_sites.nf`: export a small **genotype-free** VCF of the original site keys,
   FILTER and source INFO AC/AN/AF. No genotype reads, filtering, normalization,
   or AF recalculation. Every source ALT is retained; multiallelic rows become separate site-only allele rows, with original row and ALT-index pointers. This is not a
   reconstruction of the large genotyped VCF.
2. `annotation.nf`: consume that sites file using the existing locked FastVEP,
   picker, and LOFTEE resources. Build a sites manifest with `unit_id,chromosome,vcf`.
3. `candidates.nf`: use **unfiltered** `parquet_expanded`, not
   `parquet_expanded_mane_select`. MANE is a preference, not an inclusion rule.
4. `pre_carrier.nf`: shared current policy: source INFO AC/AN **<0.005**;
   gnomAD joint **POPMAX <0.001 or missing**; existing region exclusions.
5. `zarr_carriers.nf`: shared carrier-output writer with Zarr-backed exact allele
   lookup, sparse carrier materialization and corrected frequencies. The
   corrected-frequency population is representative samples; the unrelated
   population is representative AND unrelated. Count reference and partial
   genotypes under the existing policy; no GQ/DP/AD/AB/FT/FILTER gating here.
6. `final_rarity.nf`: unchanged; saved corrected unrelated AF **<0.001**, positive
   AN and valid/matched frequency required. No source genotype reread.
7. `post_rarity_qc.nf`: unchanged current QC and receipt validation. All sample
   IDs, untiered HC and separate sequence/star counts remain in the shared outputs.

Do not feed old 0.01/MANE-filtered test outputs into this chain. No production
SPARK scientific results are established by the synthetic tests in this PR.
The existing legacy gather rejects post-rarity-QC receipts; the adapter described
in PR #21 is still needed before chromosome-wide gathering.

## Inputs and commands

`zarr_sites.nf` manifest: `unit_id\tchromosome\tzarr`.
`zarr_carriers.nf` manifest: `unit_id\tchromosome\tzarr\tmissense\tlof_hc`, where
missense/lof_hc are the **new pre-carrier filtered Parquets**. Paths can be
absolute or relative to the manifest. Units must be unique. Supply just the
selected pilot units; neither entrypoint scans results directories for inputs.

```bash
nextflow -C zarr_sites.config,conf/site.config run zarr_sites.nf \
  --zarr_manifest manifests/pilot-zarr.tsv --outdir runs/pilot/sites \
  -work-dir runs/pilot/work-sites

# Run the existing annotation/candidates/pre_carrier entrypoints with their
# resource locks and explicit manifests; see the linked runbooks below.

nextflow -C zarr_carriers.config,conf/site.config run zarr_carriers.nf \
  --carrier_manifest manifests/pilot-filtered.tsv --psam samples.frequency.psam \
  --outdir runs/pilot/extraction -work-dir runs/pilot/work-extraction
```

Set `params.carrier_container` to an immutable SIF containing Python, Zarr 3,
NumPy, pysam and DuckDB; record its hash in the site run provenance. The older
ABCD LOFTEE SIF is **not assumed** to contain Zarr. pysam is used for BGZF output
and optional site-only VCF serialization, **not genotype decoding or AF counting**.
Zarr genotype inputs are opened read-only. Do not mutate source stores while
running or resuming. Receipts hash array metadata; they do not claim exhaustive
source chunk integrity validation. Keep source release validation evidence.

The backend groups candidate records by genotype variant chunk and reads each
needed FORMAT chunk once. Memory scales with chunk dimensions, sample count and
FORMAT arrays (particularly AD); the conservative carrier default is 16 GB.
Autosomal frequency counting is vectorized. X/Y use the same explicit scalar
PAR/ploidy helper as the VCF backend and require `grch38_x_only_par` explicitly;
confirm the dataset's PAR representation before enabling it.

## Pedigree-based unrelated proxy

`make_family_frequency_psam.py` accepts column mappings rather than cohort names.
It retains source sample order, chooses the lexicographically first sample per
participant, then the representative of the lexicographically first participant
per family. Missing family identities and conflicting participant families fail.
The JSON receipt records counts, rules, limitations and input/output hashes.
It is **not** a genetically verified unrelated subset and does not claim to
exclude cross-family relatedness. All source samples remain in carrier outputs.

SPARK metadata matched 45,178 source samples, 45,178 unique participants and
21,003 families. The authorized family-based rule selects 21,003 samples for the
unrelated-frequency proxy. Sample IDs/PSAMs remain private and are not committed.

## Validation and rollout boundary

Synthetic parity compares VCF and Zarr carrier rows, all frequency products and
frequency receipts, including reference/partial/haploid calls, overlapping
annotations, untiered HC, stars, unmatched candidates, and all zero-carrier samples.
Randomized vector/scalar counting matches. Empty products remain readable.
Real local Nextflow execution/resume and sites export are tested. The existing
VCF carrier/frequency/final-rarity/post-rarity-QC tests remain unchanged and pass.

**Not yet a validated full SPARK chr22 run.** First run a small new-policy pilot,
reconcile its candidate/extraction/frequency/rarity/QC receipts, and measure memory
before chromosome fan-out. Do not interpret the earlier successful SPARK targeted
run as scientific-policy parity; it used 0.01 and MANE-filtered scoring resources.

Multiallelic source records are supported. Each candidate ALT has its own AC,
source INFO AC/AF and carrier associations; AN includes every eligible called
allele. Other called ALTs project to zero in the candidate-specific GT. AD is
reference depth plus the selected ALT depth. Original GT/AD, zero-based source
row and one-based ALT index are retained as `source_*` carrier fields. This is
an explicit allele-specific QC representation, not a claim that another ALT is
biologically reference. Existing AB thresholds operate on reference plus target
ALT depths. Source genotypes are never rewritten. Duplicate exact alleles across
source rows fail rather than being double counted. X/Y counting checks original
multiallelic heterozygosity before applying the haploid policy.

The adapter supports FORMAT/AD; LAD/LAA-only stores need an explicit allele-depth
adapter before use. PL is not required. Sites retain the original allele spelling;
if upstream normalization changes keys, an explicit source-pointer mapping is
required. The provided sites/annotation route preserves source keys.

## Explicit genotype-derived preliminary frequencies

When source INFO AC/AN is absent, use `zarr_pre_carrier.nf` (with
`zarr_pre_carrier.config`) instead of `pre_carrier.nf`. Its filter manifest adds
an absolute `zarr` path to `unit_id,chromosome,missense,lof_hc,sites`; `--psam` is
required. Resource locks and population/region filters remain the same.

This opt-in route computes candidate-ALT AC/AN using **all source samples**, ignoring
participant-representative and unrelated selection flags for this preliminary
count only. It reads GT/mask chunks, with no genotype-quality filtering. The shared
sex/PAR counting rules require an explicit policy for X/Y. Reference and partial
calls follow the shared counting rules; zero AN and unmatched alleles cannot pass.
Counts are saved as `preliminary-frequencies.tsv` with a provenance/audit JSON.
The filter verifies the table checksum and candidate input hashes before use.
`pcf_source_info_*` remains the original source INFO (including missing values);
`pcf_cohort_*` contains the derived counts and AF. The threshold remains strictly
AF <0.005. Later representative/family-based frequencies and final unrelated
AF <0.001 are unchanged. This is an explicit frequency-source choice, never a
silent missing-INFO fallback.

The first real SPARK pilot using source INFO completed technically but retained
zero candidates because source AC/AN was missing for all candidate alleles. Those
empty outputs do not validate real carrier extraction. The genotype-derived rerun
uses a separate downstream output/work root and reuses annotation/candidates.

### Bounded FORMAT reads

Carrier extraction reads GT/mask for frequency counts across all samples. AD,
DP, GQ, FT and phasing are selected only for candidate rows and carrier samples,
with reads grouped by each FORMAT array's sample-chunk boundaries. These fields
are not materialized as whole variant-block-by-cohort arrays. Zarr still decodes
the physical chunks intersecting a selection; selected output alone is not the
peak memory requirement. Tests verify carrier/frequency parity and sample-chunk
bounded AD selections. The first whole-FORMAT pilot exceeded 16 GB; the bounded
reader is benchmarked at the same allocation, retaining preliminary counts.
