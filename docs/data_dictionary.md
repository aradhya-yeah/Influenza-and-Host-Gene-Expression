# Data Dictionary: Host Gene Expression in Influenza B Infection

## 1. Raw Data Files

### `data/raw/GSE197143_series_matrix.txt`
- **Format:** SOFT series matrix text file
- **Source:** NCBI GEO FTP repository
- **Contents:** Sample metadata headers, experimental protocols, contributor affiliations, and platform annotations.

### `data/raw/GSE197143_norm_counts.txt`
- **Format:** Tab-separated text file
- **Source:** NCBI GEO Supplementary File
- **Rows:** 19,239 gene entries
- **Columns:**
  - `ENSEMBL_gene_ID`: Ensembl gene identifier (e.g., ENSG00000000003)
  - `gene_symbol`: Standard HUGO gene symbol (e.g., DDX3X, BATF2)
  - `entrezgene_id`: NCBI Entrez Gene ID
  - `uniprotswissprot`: UniProt accession
  - `description`: Full descriptive protein name
  - `refseq_mrna`: RefSeq mRNA transcript identifier
  - `transcript_length`: Effective transcript length in base pairs
  - `Hulu_fluB_1` to `Hulu_fluB_15`: Quantile-normalized, batch-corrected log2 expression values for 15 biological samples.

---

## 2. Processed Datasets

### `data/processed/sample_metadata.csv`
- **Format:** Comma-separated values
- **Rows:** 15 sample records
- **Columns:**
  - `sample_title`: Library identifier (`Hulu_fluB_1` to `Hulu_fluB_15`)
  - `geo_accession`: NCBI GEO sample accession (`GSM5910102` to `GSM5910116`)
  - `tissue`: Anatomical origin (`lung`)
  - `infection_status`: Experimental infection status (`influenza B` or `mock`)
  - `treatment`: Specific condition (`Mock_48hpi`, `FluB_337_24hpi`, `FluB_337_48hpi`, `FluB_Bay_24hpi`, `FluB_Bay_48hpi`)
  - `group`: Simplified cohort classification (`Mock` or `Infected`)
  - `condition`: Treatment label

### `data/processed/cleaned_expression_matrix.csv`
- **Format:** Comma-separated values
- **Rows:** 19,239 filtered gene records
- **Columns:** Gene identifiers and normalized expression measurements across all 15 explant libraries.

### `data/processed/differential_expression_results.csv`
- **Format:** Comma-separated values
- **Rows:** 19,239 genes evaluated
- **Columns:**
  - `gene_symbol`: Target gene symbol
  - `ensembl_id`: Ensembl gene identifier
  - `description`: Descriptive gene function
  - `mock_mean`: Arithmetic mean of normalized log2 expression across mock controls ($n = 3$)
  - `infected_mean`: Arithmetic mean of normalized log2 expression across infected explants ($n = 12$)
  - `log2FoldChange`: Effect size ($\bar{X}_{\text{infected}} - \bar{X}_{\text{mock}}$)
  - `pvalue`: Raw p-value from Welch's two-sample t-test
  - `padj`: False Discovery Rate (FDR) adjusted p-value via Benjamini-Hochberg method
  - `significance`: Categorical call (`Upregulated`, `Downregulated`, or `Not Significant` based on $\text{padj} < 0.05$ and $|\text{log2FC}| \geq 1.0$)
