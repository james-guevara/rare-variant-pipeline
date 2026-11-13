# Gene Sets for Rare Variant Burden Analysis

Downloaded: 2026-01-15

## Quick Reference

| File | Genes | Use For |
|------|-------|---------|
| `sfari_hc.txt` | 315 | ASD high-confidence genes (score 1-2) |
| `sfari_all.txt` | 1,267 | All ASD candidate genes |
| `ddg2p_confident.txt` | 1,706 | NDD genes (definitive + strong) |
| `ddg2p_all.txt` | 2,530 | All developmental disorder genes |
| `schema_sig.txt` | 10 | SCZ exome-wide significant genes |
| `schema_all.txt` | 18 | SCZ significant + FDR<5% genes |

## Sources

### SFARI Gene (Autism)
- **URL:** https://gene.sfari.org/database/human-gene/
- **Files:** `sfari_genes.csv` (raw), `sfari_hc.txt`, `sfari_all.txt`
- **Scoring:**
  - Score 1: High confidence (≥3 de novo LGD, FDR<0.1)
  - Score 2: Strong candidate (2 de novo LGD or replicated GWAS)
  - Score 3: Suggestive (1 de novo LGD)
  - S suffix: Syndromic

### DDG2P / Gene2Phenotype (Developmental Disorders)
- **URL:** https://www.ebi.ac.uk/gene2phenotype/
- **FTP:** http://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads/
- **Files:** `DDG2P.csv` (raw, 2025-12-28 release), `ddg2p_confident.txt`, `ddg2p_all.txt`
- **Confidence levels:**
  - Definitive: Highest confidence
  - Strong: High confidence
  - Moderate: Moderate evidence
  - Limited: Limited evidence

### SCHEMA (Schizophrenia)
- **Paper:** Singh et al. Nature 2022 (doi:10.1038/s41586-022-04556-w)
- **Browser:** https://schema.broadinstitute.org/
- **Files:** `schema_scz_genes.tsv` (with metadata), `schema_sig.txt`, `schema_all.txt`
- **Tiers:**
  - `exome_wide_sig`: 10 genes at P < 2.5e-6
  - `fdr_05`: 8 additional genes at FDR < 5%

## File Formats

### Simple gene lists (*.txt)
One gene symbol per line, ready for filtering:
```
CHD8
SCN2A
SYNGAP1
...
```

### Raw data files
- `sfari_genes.csv`: status, gene-symbol, gene-name, ensembl-id, chromosome, genetic-category, gene-score, syndromic, ...
- `DDG2P.csv`: g2p id, gene symbol, gene mim, hgnc id, disease name, allelic requirement, confidence, variant consequence, ...
- `schema_scz_genes.tsv`: gene_symbol, tier, p_value, OR, description

## Usage Example (Python/Polars)

```python
import polars as pl

# Load gene set
with open("resources/gene_sets/sfari_hc.txt") as f:
    sfari_genes = set(line.strip() for line in f)

# Filter variants to gene set
df = pl.scan_csv("variants.tsv", separator="\t")
df_filtered = df.filter(pl.col("SYMBOL").is_in(sfari_genes))
```

## Citation

If using these gene sets, cite:
- SFARI: https://gene.sfari.org/about-sfari-gene/
- DDG2P: Thormann et al. 2019 (doi:10.1038/s41467-019-10016-3)
- SCHEMA: Singh et al. 2022 (doi:10.1038/s41586-022-04556-w)

### Fu et al. ASD/DD/NDD Genes
- **Source:** Fu_ASD_DD_NDD_Genes.xlsx (Fu et al. TADA analysis)
- **Files:** `fu_*.txt`
- **Thresholds based on TADA FDR (inferred from data, not explicitly stated in source):**

| File | Genes | FDR Threshold |
|------|-------|---------------|
| `fu_asd72.txt` | 72 | FDR < 0.001 (exome-wide sig) |
| `fu_asd185.txt` | 185 | FDR < 0.05 |
| `fu_dd309.txt` | 309 | FDR < 0.001 (exome-wide sig) |
| `fu_dd477.txt` | 477 | FDR < 0.05 |
| `fu_ndd373.txt` | 373 | FDR < 0.001 (exome-wide sig) |
| `fu_ndd664.txt` | 664 | FDR < 0.05 |
| `fu_satterstrom102.txt` | 102 | Satterstrom et al. 2020 |
| `fu_scz10.txt` | 9 | Top SCZ (SCHEMA) |
| `fu_scz244.txt` | 234 | Extended SCZ |
| `fu_asd_dd_ndd_combined.txt` | 675 | Union: ASD185 + DD477 + NDD664 |

**Raw data file** also contains TADA p-values, FDR, and log10 Bayes factors for different variant classes (PTV, missense, CNV).
