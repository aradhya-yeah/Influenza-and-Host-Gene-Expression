# Methodology: Transcriptomic Analysis of Host Gene Expression in Influenza B Infection

## 1. Dataset Selection and Verification

Transcriptomic data was acquired from the NCBI Gene Expression Omnibus (GEO) under accession **GSE197143** (Brunotte et al., 2022).
- **Experimental Platform:** Illumina NovaSeq 6000 (GPL24676)
- **Biological System:** Human lung tissue explants
- **Infection Model:** Influenza B virus clinical isolates (B/Brisbane/60/2008 and B/Bavaria/78/2006)
- **Timepoints:** 24 and 48 hours post-infection (hpi)
- **Controls:** Mock-infected explants harvested at 48 hpi
- **Sample Count:** 15 total biological samples (3 Mock controls, 12 Infected explants)

Raw files acquired from the NCBI FTP repository:
1. `GSE197143_series_matrix.txt.gz`: Sample metadata and experimental characteristics.
2. `GSE197143_btchCORR_norm_Brunotte_Hulu_FluB_2022_120122.txt.gz`: Batch-corrected, quantile-normalized log2-transformed expression matrix.

---

## 2. Data Cleaning and Preprocessing Protocol

1. **Metadata Parsing:** Extracted sample titles (`Hulu_fluB_1` to `Hulu_fluB_15`), GEO accessions (`GSM5910102` to `GSM5910116`), tissue source (`lung`), and experimental conditions (`Mock_48hpi`, `FluB_337_24hpi`, `FluB_337_48hpi`, `FluB_Bay_24hpi`, `FluB_Bay_48hpi`).
2. **Missing Value and Identifier Filtering:** Rows without a valid HUGO gene symbol or containing missing expression entries were excluded. Duplicate entries for the same gene symbol were resolved by retaining the first occurrence. The final expression dataset comprises 19,239 distinct genes across all 15 samples.

---

## 3. Exploratory Data Analysis

1. **Distribution Inspection:** Boxplots and kernel density estimations were computed across all samples to verify equal median alignment and absence of technical dropouts.
2. **Principal Component Analysis (PCA):** 
   - The expression matrix was transposed such that samples served as rows and genes served as columns.
   - Mean centering was performed per gene.
   - The covariance matrix across samples was calculated, and eigenvalues/eigenvectors were extracted using singular value decomposition.
   - The proportion of variance explained by PC1 and PC2 was determined to assess sample clustering.

---

## 4. Differential Gene Expression Analysis

Because the dataset provides batch-corrected and quantile-normalized log2 expression values, linear statistical modeling was applied:
1. **Group Comparison:** Samples were partitioned into two cohorts:
   - Mock Controls ($n = 3$)
   - Influenza B Infected ($n = 12$)
2. **Log2 Fold Change Calculation:**
   $$\text{log2FC} = \bar{X}_{\text{infected}} - \bar{X}_{\text{mock}}$$
   where $\bar{X}$ denotes the arithmetic mean of log2 normalized counts for a given gene.
3. **Statistical Hypothesis Testing:**
   Welch's two-sample t-test (two-sided, unequal variance) was computed for each gene:
   $$t = \frac{\bar{X}_{\text{infected}} - \bar{X}_{\text{mock}}}{\sqrt{\frac{s_{\text{infected}}^2}{n_{\text{infected}}} + \frac{s_{\text{mock}}^2}{n_{\text{mock}}}}}$$
4. **Multiple Hypothesis Correction:**
   The Benjamini-Hochberg procedure was applied to control the False Discovery Rate (FDR):
   $$\text{padj}_{(i)} = \min \left( 1, \frac{P_{(i)} \cdot m}{i} \right)$$
   enforcing monotonicity from the highest to lowest p-value.
5. **Thresholding Criteria:**
   A gene was classified as significantly differentially expressed if:
   - $\text{padj} < 0.05$
   - $|\text{log2FC}| \geq 1.0$ (representing a twofold or greater change in mRNA abundance).

---

## 5. Hierarchical Clustering and Pathway Categorization

For top differentially expressed genes, expression values were standardized using sample-wise Z-score transformation:
$$Z = \frac{x - \mu}{\sigma}$$
Hierarchical agglomerative clustering using average linkage and correlation distance was applied to order genes by co-expression similarity. Functional annotation was cross-referenced with Gene Ontology (GO) biological processes, focusing on type I interferon signaling, antiviral effector function, and cellular RNA metabolism.
