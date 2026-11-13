# Legacy gene-set source validation

The exact files in this directory were copied from the historical Expanse pipeline resource
directory on 2026-08-29. The ignored raw bundle also contains its original `sfari_genes.csv`,
`DDG2P_2025-12-28.csv.gz`, and `Fu_ASD_DD_NDD_Genes.xlsx`, plus newly recovered primary
publication supplements.

Validation results:

- All nine Fu workbook flags (`ASD72`, `ASD185`, `DD309`, `DD477`, `NDD373`, `NDD664`,
  `Satterstrom102`, `SCZ10`, and `SCZ244`) reproduce their corresponding legacy text files
  exactly.
- The primary Fu et al. Supplementary Table 11 has the same 27 columns and 18,159 gene rows
  as the compiled historical workbook. The legacy flags above reproduce exactly from it.
- The primary Satterstrom et al. `102_ASD` sheet reproduces `fu_satterstrom102.txt` exactly
  after excluding the two footer/comment rows without Entrez identifiers.
- The historical SFARI CSV reproduces `sfari_all.txt` (1,267), `sfari_s1.txt` (244), and
  `sfari_s12.txt` (952) exactly.
- The DDG2P CSV reproduces `ddg2p_all.txt` (2,530) exactly. Its simple current-style
  `definitive OR strong` filter yields 1,972 unique genes, while the historical
  `ddg2p_confident.txt` contains 1,706 and is a strict subset. Therefore that legacy file is
  retained verbatim; its additional historical exclusion rule is not yet documented and must
  not be silently reconstructed as a plain confidence filter.
