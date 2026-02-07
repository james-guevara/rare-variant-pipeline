# Rare Variant Pipeline

A modular Nextflow DSL2 pipeline for annotating and analyzing rare genetic variants in family-based sequencing studies (WGS/WES).

## Overview

This pipeline processes VCF files through variant annotation (VEP), filters for consequential variants (HIGH/MODERATE impact), and queries family genotypes to identify rare variants segregating in families.

## Pipeline Structure

```
┌─────────────────────────────────────────────────────────────────────┐
│                         VCF_PROCESSING                              │
│  BCFTOOLS_FILTER → VEP_ANNOTATE → SPLIT_VEP → REFORMAT_VARIANTS    │
│                                                    ↓                │
│                                         consequential.bed/.tsv      │
└─────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────┐
│                       FAMILY_PROCESSING                             │
│  SCATTER_BED → FAMILY_QUERY → GATHER_GENOTYPES → RESOLVE_GENOTYPES │
│       ↓                                                             │
│  (parallel chunks)                                                  │
└─────────────────────────────────────────────────────────────────────┘
                                    ↓
┌─────────────────────────────────────────────────────────────────────┐
│                         MERGE_INDEX                                 │
│              MERGE_ANNOTATIONS → SORT_INDEX                         │
│                                      ↓                              │
│                              final.tsv.gz + .tbi                    │
└─────────────────────────────────────────────────────────────────────┘
```

## Directory Structure

```
├── main.nf              # Main workflow with entry points
├── nextflow.config      # Configuration (profiles, resources, containers)
├── modules/             # Individual process definitions
│   ├── bcftools_filter.nf
│   ├── vep_annotate.nf
│   ├── split_vep.nf
│   ├── reformat_variants.nf
│   ├── scatter_bed.nf
│   ├── family_query.nf
│   ├── gather_genotypes.nf
│   ├── resolve_genotypes.nf
│   ├── merge_annotations.nf
│   └── sort_index.nf
├── subworkflows/        # Grouped process chains
│   ├── vcf_processing.nf
│   ├── family_processing.nf
│   └── merge_index.nf
└── scripts/             # Python scripts for data processing
    ├── reformat_variants.py
    ├── family_query.py
    ├── resolve_family_genotypes.py
    ├── merge_genotypes_annotations.py
    └── verify_and_gather.py
```

## Requirements

- Nextflow >= 25.10.0
- Singularity (for containerized processes)
- Python 3.12 with: polars, cyvcf2, numpy

### Containers
- `bcftools:1.22` - for filtering and split-vep
- `ensembl-vep:115.2` - for variant annotation

## Usage

### Full Pipeline
```bash
nextflow run main.nf -profile <cohort> --chroms <chromosomes> -resume
```

### Individual Subworkflows
```bash
# VCF processing only (annotation)
nextflow run main.nf -profile ssc -entry RUN_VCF_PROCESSING --chroms chr22

# Family processing only (requires VCF_PROCESSING outputs)
nextflow run main.nf -profile ssc -entry RUN_FAMILY_PROCESSING --chroms chr22

# Merge and index only (requires both previous outputs)
nextflow run main.nf -profile ssc -entry RUN_MERGE_INDEX --chroms chr22
```

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--chroms` | Comma-separated chromosomes | `chr22` |
| `--outdir` | Output directory | `output` |
| `--regions_per_chunk` | Regions per scatter chunk | `1000` |
| `--trace_prefix` | Prefix for report files | `""` |

### Profiles

- `spark_wgs` - SPARK WGS data
- `spark_iwes` - SPARK iWES v3 data  
- `ssc` - SSC data

## Output

Final output is a gzipped, tabix-indexed TSV with:
- Variant annotations (VEP consequences, impact, gene, etc.)
- Family genotypes (GT, GQ, DP, AD per sample)
- Family relationships

```
output/
├── sites/          # Filtered VCF sites
├── vep/            # VEP annotated VCFs
├── variants/       # Split VEP TSV
├── reformat/       # Reformatted + consequential variants
├── scatter/        # BED chunks
├── family_query/   # Per-chunk genotypes
├── gather/         # Gathered genotypes per chromosome
├── resolve/        # Resolved family genotypes
├── merged/         # Merged annotations + genotypes
└── indexed/        # Final gzipped + tabix indexed output
```

### Output Columns

Each row in the final TSV represents one **sample x ALT allele x variant** — multi-allelic sites and multi-member families expand into multiple rows. Columns come from three sources, joined on `(#CHROM, POS0, END, REF, ALT_specific)`.

#### Genotype columns (from `family_query.py` + `resolve_family_genotypes.py`)

| Column | Description |
|--------|-------------|
| `#CHROM` | Chromosome |
| `POS0` | 0-based start position |
| `END` | End position (exclusive) |
| `POS` | 1-based position |
| `REF` | Reference allele |
| `ALT` | All ALT alleles (comma-separated) |
| `SAMPLE` | Sample ID |
| `GT` | Genotype (e.g., `0/1`) |
| `GQ` | Genotype quality |
| `DP` | Read depth |
| `AD` | Allelic depths (comma-separated) |
| `FAMILY` | Family ID from PED file |
| `allele_num` | Which ALT allele this row represents (1-based index) |
| `ALT_specific` | The specific ALT allele for this row |
| `carrier` | 1 if sample carries this allele, 0 if not, null if missing GT |
| `AD_ref` | Reference allele depth |
| `AD_alt` | ALT allele depth (allele-specific) |
| `sample_is_multiallelic_het` | True if sample is het for >1 different ALT allele |
| `family_has_multiallelic_het` | True if any family member is multiallelic het |
| `alleles_in_gt` | Comma-separated ALT allele indices in this sample's GT |
| `is_missing_gt` | True if genotype contains `.` |

#### VEP annotation columns (from `split_vep.nf` + `reformat_variants.py`)

VCF-level fields:

| Column | Description |
|--------|-------------|
| `ID` | Variant ID |
| `QUAL` | Variant quality score |
| `INFO` | INFO field (CSQ blob stripped) |

Standard VEP fields:

| Column | Description |
|--------|-------------|
| `Allele` | VEP allele |
| `Consequence` | VEP consequence type(s) |
| `IMPACT` | HIGH, MODERATE, LOW, or MODIFIER |
| `SYMBOL` | Gene symbol |
| `Gene` | Ensembl gene ID |
| `Feature_type` | Type of feature (Transcript, RegulatoryFeature, etc.) |
| `Feature` | Ensembl transcript/feature ID |
| `BIOTYPE` | Transcript biotype |
| `EXON` | Exon number/total |
| `INTRON` | Intron number/total |
| `HGVSc` | HGVS coding notation |
| `HGVSp` | HGVS protein notation |
| `cDNA_position` | Position in cDNA |
| `CDS_position` | Position in CDS |
| `Protein_position` | Position in protein |
| `Amino_acids` | Reference/variant amino acids |
| `Codons` | Reference/variant codons |
| `Existing_variation` | Known variant IDs (e.g., rsIDs) |
| `ALLELE_NUM` | Allele number from VCF |
| `DISTANCE` | Distance to transcript |
| `STRAND` | Transcript strand |
| `FLAGS` | Transcript quality flags |
| `SYMBOL_SOURCE` | Source of gene symbol |
| `HGNC_ID` | HGNC identifier |
| `CANONICAL` | YES if canonical transcript |
| `MANE_SELECT` | MANE Select transcript ID |
| `MANE_PLUS_CLINICAL` | MANE Plus Clinical transcript ID |
| `NEAREST` | Nearest gene symbol |
| `DOMAINS` | Protein domains |
| `PUBMED` | PubMed IDs |
| `Uploaded_allele` | Allele as uploaded |

LOFTEE plugin:

| Column | Description |
|--------|-------------|
| `LoF` | LOFTEE verdict: HC, LC, or empty |
| `LoF_filter` | LOFTEE filter reason |
| `LoF_flags` | LOFTEE warning flags |
| `LoF_info` | LOFTEE detailed info |

SpliceAI plugin:

| Column | Description |
|--------|-------------|
| `SpliceAI_pred` | SpliceAI prediction scores |

MaxEntScan plugin (with SWA, NCSS):

| Column | Description |
|--------|-------------|
| `MaxEntScan_ref` | MaxEntScan reference score |
| `MaxEntScan_alt` | MaxEntScan alternate score |
| `MaxEntScan_diff` | MaxEntScan score difference |

dbNSFP plugin — pathogenicity scores:

| Column | Description |
|--------|-------------|
| `MPC_score` | MPC missense constraint score |
| `MPC_rankscore` | MPC rank score |
| `PrimateAI_score` | PrimateAI pathogenicity score |
| `PrimateAI_rankscore` | PrimateAI rank score |
| `PrimateAI_pred` | PrimateAI prediction |
| `ClinPred_score` | ClinPred pathogenicity score |
| `ClinPred_rankscore` | ClinPred rank score |
| `ClinPred_pred` | ClinPred prediction |
| `AlphaMissense_score` | AlphaMissense pathogenicity score |
| `AlphaMissense_rankscore` | AlphaMissense rank score |
| `AlphaMissense_pred` | AlphaMissense prediction |
| `CADD_raw` | CADD raw score |
| `CADD_raw_rankscore` | CADD raw rank score |
| `CADD_phred` | CADD phred-scaled score |

dbNSFP plugin — population frequencies:

| Column | Description |
|--------|-------------|
| `1000Gp3_AC` | 1000 Genomes allele count |
| `1000Gp3_AF` | 1000 Genomes allele frequency |
| `AllofUs_ALL_AF` | All of Us overall AF |
| `AllofUs_POPMAX_AF` | All of Us population max AF |
| `AllofUs_POPMAX_POP` | All of Us population max population |
| `RegeneronME_ALL_AF` | Regeneron ME overall AF |
| `gnomAD4.1_joint_flag` | gnomAD v4.1 joint flag |
| `gnomAD4.1_joint_AF` | gnomAD v4.1 joint AF |
| `gnomAD4.1_joint_nhomalt` | gnomAD v4.1 joint homozygous alt count |
| `gnomAD4.1_joint_POPMAX_AF` | gnomAD v4.1 population max AF |
| `gnomAD4.1_joint_POPMAX_nhomalt` | gnomAD v4.1 population max homozygous alt count |
| `ALFA_Total_AF` | ALFA total allele frequency |
| `dbNSFP_POPMAX_AF` | dbNSFP cross-database population max AF |
| `dbNSFP_POPMAX_AC` | dbNSFP cross-database population max AC |
| `dbNSFP_POPMAX_POP` | dbNSFP population max population |

dbNSFP plugin — ClinVar:

| Column | Description |
|--------|-------------|
| `clinvar_id` | ClinVar variant ID |
| `clinvar_clnsig` | ClinVar clinical significance |
| `clinvar_trait` | ClinVar trait name |
| `clinvar_review` | ClinVar review status |
| `clinvar_hgvs` | ClinVar HGVS notation |
| `clinvar_var_source` | ClinVar variant source |
| `clinvar_MedGen_id` | ClinVar MedGen ID |
| `clinvar_OMIM_id` | ClinVar OMIM ID |
| `clinvar_Orphanet_id` | ClinVar Orphanet ID |

dbNSFP plugin — conservation scores:

| Column | Description |
|--------|-------------|
| `GERP++_NR` | GERP neutral rate |
| `GERP++_RS` | GERP rejected substitutions |
| `GERP++_RS_rankscore` | GERP RS rank score |
| `GERP_92_mammals` | GERP score across 92 mammals |
| `GERP_92_mammals_rankscore` | GERP 92 mammals rank score |
| `phyloP100way_vertebrate` | phyloP 100-way vertebrate conservation |
| `phyloP100way_vertebrate_rankscore` | phyloP 100-way rank score |
| `phyloP470way_mammalian` | phyloP 470-way mammalian conservation |
| `phyloP470way_mammalian_rankscore` | phyloP 470-way rank score |
| `phyloP17way_primate` | phyloP 17-way primate conservation |
| `phyloP17way_primate_rankscore` | phyloP 17-way primate rank score |
| `phastCons100way_vertebrate` | phastCons 100-way vertebrate |
| `phastCons100way_vertebrate_rankscore` | phastCons 100-way rank score |
| `phastCons470way_mammalian` | phastCons 470-way mammalian |
| `phastCons470way_mammalian_rankscore` | phastCons 470-way rank score |
| `phastCons17way_primate` | phastCons 17-way primate |
| `phastCons17way_primate_rankscore` | phastCons 17-way primate rank score |

Per-allele INFO fields (extracted by `reformat_variants.py`):

| Column | Description |
|--------|-------------|
| `AC` | Allele count (allele-specific) |
| `AF` | Allele frequency (allele-specific) |
| `AQ` | Allele quality (allele-specific) |
| `MLEAC` | Maximum likelihood allele count |
| `MLEAF` | Maximum likelihood allele frequency |

#### Annotation overlay columns (from `reformat_variants.py`)

| Column | Description |
|--------|-------------|
| `segDups` | Overlap count with segmental duplications |
| `simpleRepeats` | Overlap count with simple repeats |
| `constraint_z_1kb_mean` | Mean genomic constraint z-score (1kb windows) |
| `constraint_oe_1kb_mean` | Mean genomic constraint observed/expected (1kb windows) |
| `genebayes_obs_lof` | GeneBayes observed LoF count |
| `genebayes_exp_lof` | GeneBayes expected LoF count |
| `genebayes_prior_mean` | GeneBayes prior mean |
| `genebayes_post_mean` | GeneBayes posterior mean |
| `genebayes_post_lower_95` | GeneBayes posterior 95% lower bound |
| `genebayes_post_upper_95` | GeneBayes posterior 95% upper bound |
| `gnomad_lof_pLI_tx` | gnomAD pLI for this specific transcript |
| `gnomad_lof_oe_ci_upper_tx` | gnomAD LOEUF for this specific transcript |
| `gnomad_lof_pLI_max` | gnomAD pLI max across all transcripts for gene |
| `gnomad_lof_oe_ci_upper_min` | gnomAD LOEUF min across all transcripts for gene |
| `gnomad_lof_pLI_canonical` | gnomAD pLI for canonical transcript |
| `gnomad_lof_oe_ci_upper_canonical` | gnomAD LOEUF for canonical transcript |

## Post-Pipeline Region Annotation

After the pipeline produces per-chromosome parquet files, you can add region-based annotations (repeat overlaps, blacklist flags, regulatory elements, etc.) using `annotate_regions.py`.

### Setup

1. Define tracks in the manifest (`resources/region_tracks.tsv`)
2. Download the BED files:
   ```bash
   bash scripts/download_region_tracks.sh
   ```

### Usage

```bash
# Annotate with all downloaded tracks
python scripts/annotate_regions.py output/indexed/chr22.parquet \
    -o output/annotated/chr22.parquet \
    --manifest resources/region_tracks.tsv \
    --regions-dir resources/regions/

# Annotate with specific BEDs only
python scripts/annotate_regions.py output/indexed/chr22.parquet \
    -o output/annotated/chr22.parquet \
    --regions-dir resources/regions/ \
    --beds rmsk.bed encodeBlacklist.bed windowmaskerSdust.bed

# Re-run after adding new tracks (skips already-annotated columns)
python scripts/annotate_regions.py output/annotated/chr22.parquet \
    -o output/annotated/chr22.parquet \
    --manifest resources/region_tracks.tsv \
    --regions-dir resources/regions/ \
    --skip-existing
```

Each track adds an overlap count column (e.g., `rmsk`, `encodeBlacklist`, `windowmaskerSdust`). To add a new annotation, add a row to `resources/region_tracks.tsv`, download it, and re-run with `--skip-existing`.

### Available tracks

See `resources/region_tracks.tsv` for the full list. Key tracks:

| Track | Category | Description |
|-------|----------|-------------|
| `genomicSuperDups` | filter | Segmental duplications |
| `simpleRepeat` | filter | Simple/tandem repeats (TRF) |
| `rmsk` | filter | RepeatMasker (LINE, SINE, LTR, DNA transposons) |
| `windowmaskerSdust` | filter | Low-complexity regions |
| `gap` | filter | Assembly gaps (centromeres, telomeres) |
| `encodeBlacklist` | filter | ENCODE unified blacklist (~910 artifact regions) |
| `encodeCcreCombined` | annotate | ENCODE cis-regulatory elements (promoter, enhancer, CTCF) |
| `cpgIslandExt` | annotate | CpG islands |
| `phastConsElements100way` | annotate | Conserved elements (100 vertebrates) |
| `gwasCatalog` | annotate | NHGRI-EBI GWAS Catalog hits |

## Configuration

Edit `nextflow.config` to customize:
- Container paths
- VEP cache/plugin locations
- Python environment path
- SLURM queue settings
- Resource allocations per process

## License

MIT
