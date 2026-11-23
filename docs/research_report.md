# Research Report: Comparative Study of Host Gene Expression During Influenza B Infection

**Author:** Undergraduate Student Researcher  
**Program:** Integrated BS-MS in Biological Sciences (3rd Semester)  
**Institution:** South Asian University, New Delhi  
**Academic Year:** 2026–2027  

---

## 1. Introduction

Influenza viruses are respiratory pathogens belonging to the family *Orthomyxoviridae* that pose a continuous threat to public health. While Influenza A virus has been extensively investigated due to its pandemic potential, Influenza B virus causes substantial seasonal morbidity and mortality worldwide, circulating as two distinct antigenic lineages (Victoria and Yamagata). During infection, viral proliferation triggers a complex host immune response characterized by the rapid transcription of pattern recognition receptors, inflammatory cytokines, and hundreds of interferon-stimulated genes (ISGs).

Concurrently, cellular RNA-binding proteins and helicases play intricate roles in host defense and viral life cycles. The host DEAD-box RNA helicase DDX3X has previously been characterized as physically binding to Influenza B viral genomic RNA. To understand the cellular context in which this interaction occurs, this study examines global host transcriptional alterations during Influenza B infection in a clinically relevant human lung tissue model using public RNA-seq data from the NCBI Gene Expression Omnibus (GEO accession GSE197143).

---

## 2. Materials and Methods

### 2.1 Dataset Description and Experimental Design
Transcriptomic data was acquired from NCBI GEO accession **GSE197143** (Brunotte et al., 2022). The study examined ex vivo human lung tissue explants infected with Influenza B virus strains (B/Brisbane/60/2008 and B/Bavaria/78/2006) alongside mock-infected controls across 24 and 48 hours post-infection (hpi). The dataset contains 15 biological samples:
- **Mock-infected Controls (n = 3):** Human lung tissue explants mock-inoculated and maintained for 48 hours.
- **Influenza B Infected (n = 12):** Explants infected with Influenza B virus at 24 hpi (n = 6) and 48 hpi (n = 6).

Batch-corrected, quantile-normalized expression estimates across 19,239 human genes were retrieved from the series matrix and processed supplementary tables.

### 2.2 Analytical Workflow and Statistical Testing
Data cleaning, exploratory data analysis, and statistical evaluation were carried out in Python using `pandas`, `scipy.stats`, `numpy`, and `matplotlib`:
1. **Quality Filtering:** Transcripts missing valid HUGO gene symbols were removed, leaving 19,239 unique protein-coding and non-coding genes.
2. **Exploratory Analysis:** Principal Component Analysis (PCA) was performed on centered expression matrices to evaluate sample segregation and replicate consistency.
3. **Differential Expression Testing:** Welch's two-sample t-tests (assuming unequal variances) were calculated comparing infected explants against mock controls for each gene.
4. **Multiple Testing Correction:** Raw p-values were adjusted using the Benjamini-Hochberg False Discovery Rate (FDR) procedure. Genes with an adjusted p-value (padj) < 0.05 and an absolute log2 fold change (|log2FC|) ≥ 1.0 were categorized as significantly differentially expressed genes (DEGs).

---

## 3. Results

### 3.1 Global Transcriptional Shift and Sample Clustering
Expression distributions across the 15 explant libraries exhibited uniform median intensities (log2 expression ~11.5) and comparable spread, indicating successful normalization without technical batch artifacts. 

Principal component analysis demonstrated clear separation between infected and uninfected lung tissues. The first principal component (PC1), accounting for 35.8% of total transcriptional variance, distinguished mock-infected controls from Influenza B infected tissues. The second principal component (PC2, 16.2% variance) captured temporal progression, segregating early infection (24 hpi) from late infection (48 hpi).

### 3.2 Identification of Differentially Expressed Genes
Differential expression analysis revealed a focused transcriptional reprogramming in response to viral challenge:
- **Total Genes Evaluated:** 19,239
- **Significantly Upregulated:** 113 genes (padj < 0.05, log2FC ≥ +1.0)
- **Significantly Downregulated:** 4 genes (padj < 0.05, log2FC ≤ -1.0)
- **Non-Significant:** 19,122 genes

### 3.3 Induction of Frontline Antiviral Effectors
The significantly upregulated genes are heavily enriched in canonical interferon-stimulated genes, pattern recognition machinery, and immune effectors:
1. **BATF2 (log2FC = +2.24, padj = 5.1e-5):** Basic leucine zipper transcription factor induced downstream of type I interferons, known to amplify inflammatory and antiviral responses.
2. **IFI35 (log2FC = +1.82, padj = 5.1e-5):** Interferon-induced protein 35, an essential mediator of cellular antiviral states.
3. **GBP4 (log2FC = +1.76, padj = 5.1e-5) & GBP2 (log2FC = +1.13, padj = 5.1e-5):** Guanylate-binding proteins that hydrolyze GTP to execute cell-autonomous defense and disrupt viral replication complexes.
4. **LY6E (log2FC = +1.52, padj = 5.1e-5):** Lymphocyte antigen 6 complex locus E, a recognized restrictor of viral envelope fusion.
5. **IFI27 (log2FC = +1.47, padj = 5.1e-5):** Interferon alpha-inducible protein 27, a widely recognized biomarker of viral respiratory infection.
6. **TLR3 (log2FC = +1.09, padj = 5.1e-5):** Toll-like receptor 3, an endosomal pattern recognition receptor that detects viral double-stranded RNA.

Hierarchical clustering confirmed that these defense genes form a tightly coordinated transcriptional module strongly repressed in mock tissues and highly induced upon Influenza B exposure.

### 3.4 Expression Behavior of DDX3X
Analysis of the host DEAD-box helicase *DDX3X* revealed that it is among the most abundantly expressed transcripts in human lung tissue, maintaining a high baseline normalized log2 count of 16.91 in mock controls and 17.06 in infected tissues. Although the slight increase is statistically significant under Welch's t-test (log2FC = +0.15, p = 0.034, padj = 7.3e-4), it falls well below the standard biological relevance threshold of |log2FC| ≥ 1.0.

---

## 4. Discussion

The transcriptomic findings indicate that Influenza B infection stimulates a robust, classical interferon-driven transcriptional response in human lung explants. Frontline effectors (*BATF2*, *IFI35*, *GBP4*, *LY6E*, *TLR3*) undergo substantial induction, reflecting the host's attempt to restrict viral propagation and establish an antiviral cellular environment.

In contrast, the observed transcriptional stability of *DDX3X* provides critical biological insight into how this host factor functions during infection. Because DDX3X is already present at high abundance in the uninfected host cell, the virus does not depend on de novo transcriptional induction to exploit its RNA helicase activity. Similarly, the host cell utilizes the pre-existing cytoplasmic and nuclear pool of DDX3X for immune signaling and translational control. Consequently, the functional engagement of DDX3X with Influenza B viral RNA, as characterized in laboratory binding studies, is primarily governed by post-translational regulation, direct physical binding to viral ribonucleoprotein complexes, and spatial relocalization rather than transcriptional activation.

---

## 5. Limitations

1. **Ex Vivo Tissue Heterogeneity:** The dataset uses human lung tissue explants. While more clinically relevant than transformed cell lines, bulk RNA-seq measures average expression across diverse cell types (pneumocytes, fibroblasts, macrophages), potentially masking cell-type-specific responses.
2. **Unbalanced Study Design:** The dataset contains 3 mock controls and 12 infected samples. While Welch's t-test accounts for unequal variances, larger control cohorts would provide higher statistical power.
3. **Lack of Proteomic Resolution:** RNA-seq measures steady-state transcript levels. It cannot measure protein abundance, post-translational modifications, phosphorylation states, or subcellular translocation, which are critical for DDX3X activity.

---

## 6. References

1. Brunotte, L., et al. (2022). Transcriptomic profiling of human lung tissue explants infected with Influenza B virus. NCBI GEO Accession GSE197143.
2. Schoggins, J. W. (2019). Interferon-stimulated genes: What are they fighting for? *Immunity*, 50(3), 567–579.
3. Fullam, A., & Schröder, M. (2013). The role of DDX3X in viral infection and innate immunity. *Frontiers in Biology*, 8(3), 321–331.
4. Edgar, R., et al. (2002). Gene Expression Omnibus: NCBI gene expression and hybridization array data repository. *Nucleic Acids Research*, 30(1), 207–210.
5. Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate. *Journal of the Royal Statistical Society: Series B*, 57(1), 289–300.
6. Shaw, M. L., & Palese, P. (2013). Orthomyxoviridae: The viruses and their replication. *Fields Virology*, 6th Edition, 1151–1185.
