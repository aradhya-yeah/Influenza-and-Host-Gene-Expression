# Reproducibility Guide: Host Gene Expression in Influenza B Infection

This guide provides instructions to reproduce the complete transcriptomic analysis of human lung explants infected with Influenza B virus.

## 1. Environment Setup

- **Python Version:** Python 3.10 or higher
- **Required Libraries:** `pandas`, `numpy`, `scipy`, `matplotlib`, `requests`

### Virtual Environment Setup
```bash
# Navigate to project directory
cd BIOLOGY_RESEARCH_PROJECTS/02_Influenza_Host_Gene_Expression

# Create virtual environment
python -m venv venv

# Activate virtual environment
# Windows (PowerShell):
.\venv\Scripts\Activate.ps1
# Linux / macOS:
source venv/bin/activate

# Install exact dependencies
pip install -r requirements.txt
```

---

## 2. Step-by-Step Execution

### Option A: End-to-End Orchestrator
Execute the complete pipeline sequentially:
```bash
python main.py
```

### Option B: Modular Execution

1. **Download Raw GEO Datasets:**
   ```bash
   python scripts/01_download_geo.py
   ```
   - Connects directly to NCBI FTP servers.
   - Downloads series matrix (`GSE197143_series_matrix.txt.gz`) and normalized expression matrix (`GSE197143_norm_counts.txt.gz`).
   - Automatically uncompresses gzip archives into `data/raw/`.

2. **Clean and Structure Data:**
   ```bash
   python scripts/02_clean_and_inspect.py
   ```
   - Parses sample characteristics and groups (Mock vs Infected).
   - Filters missing values and checks for duplicate identifiers.
   - Outputs `sample_metadata.csv` and `cleaned_expression_matrix.csv` into `data/processed/`.

3. **Compute Differential Gene Expression:**
   ```bash
   python scripts/03_differential_expression.py
   ```
   - Computes group means and log2 fold changes for all 19,239 genes.
   - Performs Welch's t-test and applies Benjamini-Hochberg FDR correction.
   - Exports complete statistical table to `data/processed/differential_expression_results.csv`.

4. **Generate Publication-Quality Visualizations:**
   ```bash
   python scripts/04_generate_figures.py
   ```
   - Produces all 6 figures at 300 DPI into `figures/`.

---

## 3. Key Numerical Verification Benchmarks

Upon successful execution, verify that your output numbers match these exact benchmarks:
- **Total Genes Evaluated:** 19,239
- **Significantly Upregulated:** 113 genes (FDR < 0.05, log2FC ≥ 1.0)
- **Significantly Downregulated:** 4 genes (FDR < 0.05, log2FC ≤ -1.0)
- **Top Induced Gene:** *BATF2* (log2FC = +2.24, padj = 5.09e-5)
- **DDX3X Baseline:** Mock mean = 16.91, Infected mean = 17.06, log2FC = +0.15, padj = 7.30e-4 (classified as Not Significant by fold-change threshold)
