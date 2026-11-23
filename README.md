# Comparative Study of Host Gene Expression During Influenza B Infection

**Author:** Undergraduate Student Researcher  
**Institution:** South Asian University, New Delhi  
**Program:** Integrated BS-MS in Biological Sciences (3rd Semester)  
**Research Area:** Molecular Virology, Host-Pathogen Interactions, RNA-seq Transcriptomics  

---

## Abstract

Influenza viruses perturb host cellular transcriptional networks to promote viral replication and evade intrinsic immune barriers. Following previous internship work investigating the physical interaction between the DEAD-box helicase DDX3X and Influenza B viral RNA, this study analyzes transcriptome-wide host gene expression changes in human lung tissue explants challenged with Influenza B virus using public RNA-seq data from NCBI GEO (accession GSE197143). Analysis across 15 biological samples (12 Influenza B infected and 3 mock-infected controls across 24 and 48 hours post-infection) reveals substantial transcriptional remodeling. Differential expression testing identified 113 significantly upregulated genes and 4 downregulated genes (Benjamini-Hochberg FDR < 0.05, |log2 fold change| ≥ 1.0). Upregulated transcripts are predominantly interferon-stimulated genes (ISGs) and innate pattern recognition factors, including *BATF2* (log2FC = +2.24), *IFI35* (log2FC = +1.82), *GBP4* (log2FC = +1.76), *LY6E* (log2FC = +1.52), *IFI27* (log2FC = +1.47), and *TLR3* (log2FC = +1.09). In contrast, *DDX3X* exhibits high constitutive expression in human lung explants (normalized log2 expression ≈ 17.0) with minimal transcriptional alteration (log2FC = +0.15, padj = 7.3e-4), demonstrating that its role in Influenza B infection is primarily mediated through post-transcriptional mechanisms, spatial recruitment, and direct ribonucleoprotein interactions rather than large-scale de novo transcriptional induction.

---

## Background and Research Questions

During influenza infection, host cells initiate an innate antiviral response driven by pattern recognition receptors (such as RIG-I and TLR3) that sense viral nucleic acids. This triggers signaling cascades leading to the production of type I and III interferons and the induction of hundreds of antiviral effectors. RNA helicases such as DDX3X participate in multiple stages of this process, functioning in translation regulation and viral RNA recognition.

### Research Questions
1. How does Influenza B infection alter global host gene expression in physiologically relevant human lung explants?
2. Which functional pathways and specific antiviral effector genes undergo the greatest transcriptional induction?
3. What is the baseline expression and transcriptional behavior of the host RNA helicase *DDX3X* during Influenza B infection?

---

## Key Results and Figures

### 1. Sample Expression Profiles and Quality Assessment
Boxplots of quantile-normalized log2 counts show consistent distributions across all 15 explant libraries without aberrant technical batch outliers. Group density profiles reveal a slight global shift in the transcriptome of infected tissues.

![Figure 1: Expression Distribution](figures/fig1_sample_expression_distribution.png)

### 2. Sample Segregation via Principal Component Analysis
PCA demonstrates distinct separation between mock-treated lung tissues and Influenza B infected explants along PC1 (accounting for 35.8% of total variance), while PC2 (16.2% variance) resolves progression between early (24 hpi) and late (48 hpi) infection timepoints.

![Figure 2: PCA Plot](figures/fig2_pca_plot.png)

### 3. Transcriptome-Wide Differential Expression
A volcano plot identifies a focused antiviral signature of 113 significantly induced transcripts, dominated by classical interferon effectors and viral RNA sensing components.

![Figure 3: Volcano Plot](figures/fig3_volcano_plot.png)

### 4. Hierarchical Clustering of Antiviral Programs
Bidirectional hierarchical clustering groups the top differentially expressed genes into cohesive modules representing interferon-stimulated genes (ISGs), guanylate-binding proteins, and cell death regulators.

![Figure 4: Heatmap](figures/fig4_heatmap_degs.png)

### 5. Quantitative Profiling of Key Defense Effectors and DDX3X
Comparison of normalized expression confirms significant induction of frontline defenses (*BATF2*, *IFI35*, *LY6E*, *GBP4*, *TLR3*, *DTX3L*). In contrast, *DDX3X* maintains a stable, highly abundant baseline across all samples, showing only subtle modulation (log2FC = +0.15).

![Figure 5: Antiviral Effectors vs DDX3X](figures/fig5_antiviral_genes_comparison.png)

![Figure 6: Functional Overview](figures/fig6_biological_functions_summary.png)

---

## Directory Structure

```
02_Influenza_Host_Gene_Expression/
├── README.md
├── requirements.txt
├── main.py
├── data/
│   ├── raw/
│   │   ├── GSE197143_series_matrix.txt
│   │   └── GSE197143_norm_counts.txt
│   └── processed/
│       ├── sample_metadata.csv
│       ├── cleaned_expression_matrix.csv
│       └── differential_expression_results.csv
├── figures/
│   ├── fig1_sample_expression_distribution.png
│   ├── fig2_pca_plot.png
│   ├── fig3_volcano_plot.png
│   ├── fig4_heatmap_degs.png
│   ├── fig5_antiviral_genes_comparison.png
│   └── fig6_biological_functions_summary.png
├── scripts/
│   ├── 01_download_geo.py
│   ├── 02_clean_and_inspect.py
│   ├── 03_differential_expression.py
│   └── 04_generate_figures.py
├── notebooks/
│   └── exploratory_influenza_rnaseq.ipynb
├── docs/
│   ├── research_report.md
│   ├── methodology.md
│   ├── data_dictionary.md
│   └── reproducibility_guide.md
└── references/
    ├── references.md
    └── report.tex
```

---

## How to Run

1. **Setup environment:**
   ```bash
   cd BIOLOGY_RESEARCH_PROJECTS/02_Influenza_Host_Gene_Expression
   pip install -r requirements.txt
   ```

2. **Run end-to-end pipeline:**
   ```bash
   python main.py
   ```
   This downloads GSE197143 raw data, performs quality filtering, differential expression testing with Benjamini-Hochberg FDR correction, and outputs all figures.

3. **Compile LaTeX Report (optional):**
   ```bash
   cd references
   pdflatex report.tex
   ```

---

## Limitations

- The dataset utilizes human lung explant tissue which contains diverse cell types (epithelial, endothelial, immune); bulk RNA-seq measures the average expression across all cell populations rather than cell-type-specific responses.
- The control group contains 3 mock biological replicates while the infected group contains 12 samples across two strains and timepoints, creating an unbalanced design.
- Transcriptional profiling measures steady-state mRNA abundance and does not reflect post-translational modifications, protein degradation rates, or the physical subcellular translocation of proteins such as DDX3X into viral replication factories.

---

## References

1. Brunotte, L., et al. (2022). Transcriptomic profiling of human lung tissue explants infected with Influenza B virus. NCBI GEO Accession GSE197143.
2. Schoggins, J. W. (2019). Interferon-stimulated genes: What are they fighting for? *Immunity*, 50(3), 567–579.
3. Fullam, A., & Schröder, M. (2013). The role of DDX3X in viral infection and innate immunity. *Frontiers in Biology*, 8(3), 321–331.
4. Edgar, R., et al. (2002). Gene Expression Omnibus: NCBI gene expression and hybridization array data repository. *Nucleic Acids Research*, 30(1), 207–210.
5. Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*, 57(1), 289–300.
