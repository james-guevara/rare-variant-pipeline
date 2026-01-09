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

## Configuration

Edit `nextflow.config` to customize:
- Container paths
- VEP cache/plugin locations
- Python environment path
- SLURM queue settings
- Resource allocations per process

## License

MIT
