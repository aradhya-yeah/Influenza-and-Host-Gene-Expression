from pathlib import Path
import pandas as pd
import numpy as np
from scipy import stats

project_root = Path(__file__).resolve().parent.parent
processed_dir = project_root / "data" / "processed"

expr_df = pd.read_csv(processed_dir / "cleaned_expression_matrix.csv")
meta_df = pd.read_csv(processed_dir / "sample_metadata.csv")

sample_cols = [c for c in expr_df.columns if c.startswith("Hulu_fluB_")]
mock_samples = meta_df[meta_df["group"] == "Mock"]["sample_title"].tolist()
infected_samples = meta_df[meta_df["group"] == "Infected"]["sample_title"].tolist()

expr_df = expr_df.dropna(subset=["gene_symbol"])
expr_df = expr_df.drop_duplicates(subset=["gene_symbol"], keep="first")

mock_vals = expr_df[mock_samples].values
infected_vals = expr_df[infected_samples].values

mock_means = np.mean(mock_vals, axis=1)
infected_means = np.mean(infected_vals, axis=1)

log2_fc = infected_means - mock_means

ttest_res = stats.ttest_ind(infected_vals, mock_vals, axis=1, equal_var=False)
p_values = ttest_res.pvalue

p_values = np.nan_to_num(p_values, nan=1.0)

n_genes = len(p_values)
order = np.argsort(p_values)
ranks = np.empty(n_genes, int)
ranks[order] = np.arange(1, n_genes + 1)
p_adj = np.minimum(1.0, p_values * n_genes / ranks)
p_adj = np.minimum.accumulate(p_adj[::-1])[::-1]

dge_df = pd.DataFrame({
    "gene_symbol": expr_df["gene_symbol"].values,
    "ensembl_id": expr_df["ENSEMBL_gene_ID"].values,
    "description": expr_df["description"].values,
    "mock_mean": mock_means,
    "infected_mean": infected_means,
    "log2FoldChange": log2_fc,
    "pvalue": p_values,
    "padj": p_adj
})

dge_df["significance"] = "Not Significant"
dge_df.loc[(dge_df["padj"] < 0.05) & (dge_df["log2FoldChange"] >= 1.0), "significance"] = "Upregulated"
dge_df.loc[(dge_df["padj"] < 0.05) & (dge_df["log2FoldChange"] <= -1.0), "significance"] = "Downregulated"

dge_df = dge_df.sort_values(by="padj", ascending=True)
dge_df.to_csv(processed_dir / "differential_expression_results.csv", index=False)
