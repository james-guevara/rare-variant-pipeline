# Rare Variant Annotation Pipeline

A comprehensive bioinformatics pipeline for annotating and analyzing rare genetic variants in whole genome sequencing (WGS) data from family-based studies.

## Overview

This pipeline processes joint-called VCF files through multiple stages to identify, annotate, and characterize rare variants with potential clinical or functional significance. It integrates annotations from VEP (Variant Effect Predictor), constraint metrics, conservation scores, and family-based genotype analysis.

## Pipeline Workflow

### 1. Variant Annotation (VEP)

**Script:** `RUN_ANNOTATION_PIPELINE.sb.sh`

Runs Ensembl VEP on chromosome-specific VCF files with comprehensive annotation plugins:

- **VEP Core Annotations**: Consequence predictions, gene symbols, transcript IDs, protein changes
- **LOFTEE**: Loss-of-function variant predictions with quality filters
- **SpliceAI**: Deep learning-based splice variant predictions
- **MaxEntScan**: Maximum entropy splice site scoring
- **dbNSFP v5.3a**: Extensive database including:
  - Pathogenicity scores (CADD, PrimateAI, ClinPred, AlphaMissense, MPC)
  - Population frequencies (gnomAD v4.1, 1000G, AllofUs, ALFA)
  - Conservation scores (GERP++, phyloP, phastCons)
  - ClinVar annotations

**Usage:**
```bash
# Submit as SLURM array job (chromosomes 1-24, where 23=X, 24=Y)
sbatch --array=1-24 RUN_ANNOTATION_PIPELINE.sb.sh
```

**Inputs:**
- VCF files: `vcfs/chr{CHR}_jointcall_VQSR_combined.vcf.gz`
- Reference genome: GRCh38
- VEP cache and plugins (via symlinks)

**Outputs:**
- Annotated VCF: `out_vep/chr{CHR}.vep.vcf.gz`

### 2. Post-Processing

**Script:** `RUN_ANNOTATION_PIPELINE.postprocess.sh`

Extracts VEP annotations from VCF to TSV format and filters for missense+ variants.

**Usage:**
```bash
./RUN_ANNOTATION_PIPELINE.postprocess.sh
```

**Outputs:**
- Per-chromosome TSV: `out/chr{CHR}.variants.tsv.gz`

### 3. Variant Reformatting and Enrichment

**Script:** `reformat_variants.py`

Processes multi-allelic variants and adds additional constraint and genomic annotations.

**Features:**
- Resolves multi-allelic sites (splits ALT alleles)
- Extracts per-allele INFO fields (AC, AF, AQ, MLEAC, MLEAF)
- Adds gnomAD v4.1 gene and transcript constraint metrics (pLI, LOEUF)
- Adds GeneBayes gene-level constraint
- Adds 1kb genomic constraint scores (z-scores, observed/expected ratios)
- Flags variants overlapping segmental duplications
- Flags variants in simple repeats
- Creates merged BED file of variant intervals

**Usage:**
```bash
python reformat_variants.py input.tsv.gz output.tsv.gz
```

**Inputs:**
- VEP-annotated TSV from post-processing
- Constraint databases (via `resources/` directory)

**Outputs:**
- Enriched variant TSV
- Merged variant intervals BED file

### 4. Family-Based Genotype Extraction

**Script:** `family_query.py`

Extracts genotypes for entire families when any family member carries a variant in specified genomic regions.

**Key Features:**
- Efficient family-based querying (outputs all family members when any member is a carrier)
- BED file-based region filtering
- Extracts genotype quality metrics (GT, GQ, DP, AD)
- Handles multi-allelic variants

**Usage:**
```bash
python family_query.py \
  --vcf vcfs/chr22_jointcall_VQSR_combined.vcf.gz \
  --ped pedigree.ped \
  --region regions.bed \
  --out output.genotypes.tsv
```

**Input Formats:**

*PED file* (tab-delimited, columns: FAM, IID, ...):
```
FAM001  SAMPLE001  ...
FAM001  SAMPLE002  ...
FAM002  SAMPLE003  ...
```

*BED file* (tab-delimited, 0-based coordinates):
```
chr1    1000    2000
chr1    5000    6000
```

**Outputs:**
- TSV with columns: `#CHROM, POS0, END, POS, REF, ALT, SAMPLE, GT, GQ, DP, AD, FAMILY`

### 5. Genotype Resolution for Families

**Script:** `resolve_family_genotypes.py`

Explodes family genotype data by active ALT alleles and adds carrier/compound heterozygote flags.

**Features:**
- Identifies active alleles per family
- Creates carrier flags per allele
- Detects compound heterozygotes (sample-level and family-level)
- Resolves ALT alleles and allelic depths (AD_ref, AD_alt)

**Usage:**
```bash
python resolve_family_genotypes.py input.tsv output.tsv
```

### 6. Merge Family and Variant Data

**Script:** `merge_family_variants.py`

Joins family genotype data with variant annotations.

**Usage:**
```bash
python merge_family_variants.py \
  --family family_genotypes.tsv \
  --variants annotated_variants.tsv \
  --out merged_output.tsv
```

### 7. Concatenate Chromosomes

**Script:** `concat_family_merged.py`

Concatenates per-chromosome results into genome-wide tables.

**Usage:**
```bash
python concat_family_merged.py \
  --pattern "out/chr*.family_merged.tsv.gz" \
  --out genome_wide_results.tsv.gz
```

## Dependencies

### Software Requirements

- **Python 3.9+** with packages:
  - `polars` (dataframe operations)
  - `polars-bio` (genomic interval operations)
  - `cyvcf2` (VCF parsing)
  - `numpy`

- **Singularity/Apptainer** (for containerized tools)

- **Containers:**
  - `ensembl-vep_115.2--pl5321h2a3209d_1.with_samtools`
  - `bcftools:1.22--h3a4d415_1`

- **System tools:**
  - `bgzip` / `tabix`
  - `bcftools`

### Reference Data

Required resources (typically symlinked):

- `VEP_CACHE/`: VEP annotation cache (GRCh38)
- `VEP_PLUGINS/`: VEP plugins directory
- `resources/`: Constraint databases and annotation files
  - `gnomAD/constraint/gnomad.v4.1.constraint_metrics.tsv`
  - `gnomAD/Genomic_constraint/constraint_z_genome_1kb.qc.download.txt.gz`
  - `GeneBayes/output/Supplementary_Table_1.tsv`
  - `dbNSFP/dbNSFP5.3a_grch38.gz`
  - `LOFTEE/` (human ancestor FA, GERP BigWig, conservation SQL)
  - `SpliceAI/` (SNV and indel score databases)
  - `MaxEntScan/fordownload/`
  - `repeats/` (segmental duplications, simple repeats)

## Directory Structure

```
rare_variant_pipeline/
├── RUN_ANNOTATION_PIPELINE.sb.sh          # Main VEP annotation (SLURM)
├── RUN_ANNOTATION_PIPELINE.postprocess.sh # Extract VEP to TSV
├── family_query.py                        # Family-based genotype extraction
├── reformat_variants.py                   # Multi-allelic resolution + enrichment
├── resolve_family_genotypes.py            # Explode by allele + carrier flags
├── merge_family_variants.py               # Join family + variant data
├── concat_family_merged.py                # Concatenate chromosomes
├── resolve_genotypes.py                   # (Alternative single-sample resolver)
├── containers/                            # Singularity containers
├── logs_*/                                # SLURM job logs
├── out/                                   # Final outputs
├── out_vep/                               # VEP intermediate outputs
├── resources/                             # Annotation databases (symlink)
├── vcfs/                                  # Input VCF files (symlink)
└── README.md                              # This file
```

## Quick Start

1. **Setup environment:**
   ```bash
   # Install Python dependencies
   pip install polars polars-bio cyvcf2 numpy

   # Set up symlinks to resources
   ln -s /path/to/resources resources
   ln -s /path/to/vcfs vcfs
   ln -s /path/to/VEP_CACHE VEP_CACHE
   ln -s /path/to/VEP_PLUGINS VEP_PLUGINS
   ```

2. **Run VEP annotation:**
   ```bash
   sbatch --array=1-24 RUN_ANNOTATION_PIPELINE.sb.sh
   ```

3. **Extract annotations to TSV:**
   ```bash
   ./RUN_ANNOTATION_PIPELINE.postprocess.sh
   ```

4. **Reformat and enrich variants:**
   ```bash
   for chr in {1..22} X Y; do
     python reformat_variants.py \
       out/chr${chr}.variants.tsv.gz \
       out/chr${chr}.variants.enriched.tsv.gz
   done
   ```

5. **Extract family genotypes** (if doing family analysis):
   ```bash
   # Create BED file with regions of interest first
   python family_query.py \
     --vcf vcfs/chr22_jointcall_VQSR_combined.vcf.gz \
     --ped samples.ped \
     --region regions.bed \
     --out out/chr22.family_genotypes.tsv
   ```

6. **Resolve and merge** (as needed for your analysis)

## Output Format

### Final Variant Table Columns

Core variant fields:
- `#CHROM`, `POS0`, `END`, `POS`, `REF`, `ALT`, `ALT_specific`
- `QUAL`, `INFO`

VEP annotations:
- `Consequence`, `SYMBOL`, `Gene`, `Feature`, `BIOTYPE`
- `MANE_SELECT`, `CANONICAL`
- `Protein_position`, `Amino_acids`, `Codons`
- `CADD_phred`, `CADD_raw`
- `gnomAD4.1_joint_AF`, `gnomAD4.1_joint_POPMAX_AF`
- `ClinVar` fields, pathogenicity scores, etc.

Constraint metrics:
- `gnomad_lof_pLI_canonical`, `gnomad_lof_oe_ci_upper_canonical`
- `genebayes_post_mean`, `genebayes_obs_lof`, `genebayes_exp_lof`
- `constraint_z_max`, `constraint_oe_min`

Overlap flags:
- `segDups`, `simpleRepeats`

## Cleanup Recommendations

Before committing to GitHub, remove temporary files:

```bash
# Remove core dumps (~27GB total)
rm core.*

# Remove test/temporary files
rm tmp.tsv duplicated_rows.tsv test.variants.tsv

# Clean up logs (optional - can keep for debugging)
# rm -rf logs_*/*.out logs_*/*.err
```

## Citation and References

If using this pipeline, please cite:

- **VEP**: McLaren et al. (2016) Genome Biology
- **LOFTEE**: Karczewski et al. (2020) Nature
- **dbNSFP**: Liu et al. (2020) Nucleic Acids Research
- **SpliceAI**: Jaganathan et al. (2019) Cell
- **gnomAD**: Chen et al. (2024) Nature
- **Polars**: https://pola.rs

## License

[Specify your license here]

## Contact

[Your contact information]
