from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from scipy.cluster.hierarchy import linkage, dendrogram

project_root = Path(__file__).resolve().parent.parent
processed_dir = project_root / "data" / "processed"
figures_dir = project_root / "figures"
figures_dir.mkdir(exist_ok=True)

expr_df = pd.read_csv(processed_dir / "cleaned_expression_matrix.csv")
meta_df = pd.read_csv(processed_dir / "sample_metadata.csv")
dge_df = pd.read_csv(processed_dir / "differential_expression_results.csv")

sample_cols = [c for c in expr_df.columns if c.startswith("Hulu_fluB_")]
mock_cols = meta_df[meta_df["group"] == "Mock"]["sample_title"].tolist()
inf_cols = meta_df[meta_df["group"] == "Infected"]["sample_title"].tolist()

plt.rcParams["font.sans-serif"] = "DejaVu Sans"
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["figure.dpi"] = 300

# 1. Figure 1: Sample Expression Distribution (Boxplots & Density)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5), gridspec_kw={"width_ratios": [1.4, 1.0]})

plot_data = [expr_df[col].dropna().values for col in sample_cols]
colors_list = ["#94A3B8" if col in mock_cols else "#0284C7" for col in sample_cols]

bp = ax1.boxplot(plot_data, patch_artist=True, showfliers=False, widths=0.6)
for patch, color in zip(bp["boxes"], colors_list):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
    patch.set_edgecolor("#1E293B")
for median in bp["medians"]:
    median.set(color="#0F172A", linewidth=1.2)

ax1.set_xticklabels([c.replace("Hulu_fluB_", "S") for c in sample_cols], rotation=45, fontsize=8)
ax1.set_ylabel("Normalized Log2 Expression", fontsize=9, fontweight="bold")
ax1.set_title("A. Expression Quantile Profiles Across 15 Human Lung Explant Samples", fontsize=9.5, fontweight="bold", color="#0F2C59")
ax1.grid(True, linestyle=":", alpha=0.5)

legend_elements = [
    Patch(facecolor="#94A3B8", edgecolor="#1E293B", alpha=0.7, label="Mock (Uninfected Controls, n=3)"),
    Patch(facecolor="#0284C7", edgecolor="#1E293B", alpha=0.7, label="Influenza B Infected (n=12)")
]
ax1.legend(handles=legend_elements, loc="lower left", fontsize=7.5)

# Density estimation across groups
mock_all = expr_df[mock_cols].values.flatten()
inf_all = expr_df[inf_cols].values.flatten()

ax2.hist(mock_all, bins=50, density=True, alpha=0.45, color="#64748B", label="Mock Controls")
ax2.hist(inf_all, bins=50, density=True, alpha=0.45, color="#0284C7", label="Influenza B")
ax2.set_xlabel("Normalized Expression Value", fontsize=9, fontweight="bold")
ax2.set_ylabel("Density", fontsize=9, fontweight="bold")
ax2.set_title("B. Overall Group Density Distribution", fontsize=9.5, fontweight="bold", color="#0F2C59")
ax2.grid(True, linestyle=":", alpha=0.5)
ax2.legend(loc="upper right", fontsize=8)

plt.tight_layout()
plt.savefig(figures_dir / "fig1_sample_expression_distribution.png", dpi=300)
plt.close()

# 2. Figure 2: Principal Component Analysis (PCA)
expr_matrix = expr_df[sample_cols].values.T
expr_centered = expr_matrix - np.mean(expr_matrix, axis=0)
cov_matrix = np.cov(expr_centered, rowvar=False)
eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
sorted_idx = np.argsort(eigenvalues)[::-1]
eigenvalues = eigenvalues[sorted_idx]
eigenvectors = eigenvectors[:, sorted_idx]

pca_scores = np.dot(expr_centered, eigenvectors[:, :2])
var_exp1 = (eigenvalues[0] / np.sum(eigenvalues)) * 100
var_exp2 = (eigenvalues[1] / np.sum(eigenvalues)) * 100

fig, ax = plt.subplots(figsize=(7.5, 6))
group_color_map = {"Mock": "#64748B", "Infected": "#0284C7"}

for idx, sample in enumerate(sample_cols):
    row_meta = meta_df[meta_df["sample_title"] == sample].iloc[0]
    grp = row_meta["group"]
    trt = row_meta["treatment"]
    col = group_color_map[grp]
    marker = "o" if "Mock" in grp else ("^" if "24hpi" in trt else "s")
    ax.scatter(pca_scores[idx, 0], pca_scores[idx, 1], color=col, s=80, marker=marker,
               edgecolor="#0F172A", linewidth=1.0, alpha=0.85, zorder=4)
    ax.annotate(sample.replace("Hulu_fluB_", "S"), (pca_scores[idx, 0] + 1.2, pca_scores[idx, 1] + 1.0),
                fontsize=7.5, color="#334155")

ax.axhline(0, color="#CBD5E1", linestyle="--", linewidth=0.8)
ax.axvline(0, color="#CBD5E1", linestyle="--", linewidth=0.8)
ax.set_xlabel(f"Principal Component 1 ({var_exp1:.1f}% Variance)", fontsize=9, fontweight="bold")
ax.set_ylabel(f"Principal Component 2 ({var_exp2:.1f}% Variance)", fontsize=9, fontweight="bold")
ax.set_title("Figure 2: Principal Component Analysis of Human Lung Explants\nSegregation of Uninfected Controls and Influenza B Infection Timepoints",
             fontsize=10, fontweight="bold", color="#0F2C59")
ax.grid(True, linestyle=":", alpha=0.5)

legend_pca = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#64748B', markeredgecolor='#0F172A', markersize=9, label='Mock (Control)'),
    plt.Line2D([0], [0], marker='^', color='w', markerfacecolor='#0284C7', markeredgecolor='#0F172A', markersize=9, label='Influenza B (24 hpi)'),
    plt.Line2D([0], [0], marker='s', color='w', markerfacecolor='#0284C7', markeredgecolor='#0F172A', markersize=9, label='Influenza B (48 hpi)')
]
ax.legend(handles=legend_pca, loc="upper right", fontsize=8, framealpha=0.9)

plt.tight_layout()
plt.savefig(figures_dir / "fig2_pca_plot.png", dpi=300)
plt.close()

# 3. Figure 3: Volcano Plot
fig, ax = plt.subplots(figsize=(8.5, 6))

non_sig = dge_df[dge_df["significance"] == "Not Significant"]
up_sig = dge_df[dge_df["significance"] == "Upregulated"]
down_sig = dge_df[dge_df["significance"] == "Downregulated"]

log10_padj_all = -np.log10(np.clip(dge_df["padj"].values, 1e-15, 1.0))

ax.scatter(non_sig["log2FoldChange"], -np.log10(np.clip(non_sig["padj"], 1e-15, 1.0)),
           color="#94A3B8", s=12, alpha=0.45, label=f"Not Significant ({len(non_sig):,})")
ax.scatter(up_sig["log2FoldChange"], -np.log10(np.clip(up_sig["padj"], 1e-15, 1.0)),
           color="#DC2626", s=25, alpha=0.85, label=f"Upregulated (n = {len(up_sig)})")
ax.scatter(down_sig["log2FoldChange"], -np.log10(np.clip(down_sig["padj"], 1e-15, 1.0)),
           color="#2563EB", s=35, alpha=0.85, label=f"Downregulated (n = {len(down_sig)})")

ax.axhline(-np.log10(0.05), color="#0F172A", linestyle="--", linewidth=0.8, alpha=0.7)
ax.axvline(1.0, color="#0F172A", linestyle=":", linewidth=0.8, alpha=0.7)
ax.axvline(-1.0, color="#0F172A", linestyle=":", linewidth=0.8, alpha=0.7)

# Annotate prominent antiviral genes and DDX3X
top_label_genes = ["BATF2", "IFI35", "LY6E", "GBP4", "IFI27", "TLR3", "DTX3L", "GBP2", "DDX3X"]
for gene in top_label_genes:
    g_row = dge_df[dge_df["gene_symbol"] == gene]
    if not g_row.empty:
        fc = g_row["log2FoldChange"].values[0]
        pval = -np.log10(max(g_row["padj"].values[0], 1e-15))
        is_ddx = (gene == "DDX3X")
        ax.scatter(fc, pval, color="#F59E0B" if is_ddx else "#0F172A", s=80 if is_ddx else 40,
                   edgecolor="black", zorder=6)
        label_text = f"DDX3X (FC={fc:.2f}, padj={g_row['padj'].values[0]:.4f})" if is_ddx else gene
        ax.annotate(label_text, (fc, pval), xytext=(fc + (0.3 if fc > 0 else -0.8), pval + 0.3),
                    fontsize=8, fontweight="bold", color="#B45309" if is_ddx else "#0F172A",
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="#FEF3C7" if is_ddx else "white", alpha=0.85, edgecolor="#CBD5E1"),
                    arrowprops=dict(arrowstyle="->", color="#475569", lw=0.7))

ax.set_xlabel("Log2 Fold Change (Influenza B vs Mock)", fontsize=9, fontweight="bold")
ax.set_ylabel("-Log10 Adjusted P-Value", fontsize=9, fontweight="bold")
ax.set_title("Figure 3: Transcriptome-wide Differential Expression in Human Lung Tissue\nIdentification of Key Antiviral Interferon-Stimulated Genes (FDR < 0.05, |Log2FC| ≥ 1.0)",
             fontsize=10, fontweight="bold", color="#0F2C59")
ax.grid(True, linestyle=":", alpha=0.5)
ax.legend(loc="upper left", fontsize=8)

plt.tight_layout()
plt.savefig(figures_dir / "fig3_volcano_plot.png", dpi=300)
plt.close()

# 4. Figure 4: Hierarchical Clustering Heatmap of Selected DEGs
top_genes = dge_df[dge_df["significance"] != "Not Significant"].head(30)["gene_symbol"].tolist()
if "DDX3X" not in top_genes:
    top_genes.append("DDX3X")

heat_df = expr_df[expr_df["gene_symbol"].isin(top_genes)].set_index("gene_symbol")[sample_cols]
# Z-score standardization across samples for each gene
heat_z = heat_df.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

# Hierarchical clustering of genes
linkage_matrix = linkage(heat_z.values, method="average", metric="correlation")
dendro = dendrogram(linkage_matrix, labels=heat_z.index, no_plot=True)
ordered_genes = [heat_z.index[i] for i in dendro["leaves"]]
heat_z = heat_z.loc[ordered_genes]

fig, ax = plt.subplots(figsize=(9, 8))
im = ax.imshow(heat_z.values, cmap="RdBu_r", aspect="auto", vmin=-2.5, vmax=2.5)

ax.set_xticks(range(len(sample_cols)))
ax.set_xticklabels([f"{c.replace('Hulu_fluB_', 'S')} ({'Mock' if c in mock_cols else 'Inf'})" for c in sample_cols],
                   rotation=90, fontsize=7.5)
ax.set_yticks(range(len(ordered_genes)))
ax.set_yticklabels(ordered_genes, fontsize=7.5)

# Highlight DDX3X label in red
for tick_label in ax.get_yticklabels():
    if tick_label.get_text() == "DDX3X":
        tick_label.set_color("#DC2626")
        tick_label.set_fontweight("bold")

cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
cbar.set_label("Relative Expression (Z-score)", fontsize=8, fontweight="bold")

ax.set_title("Figure 4: Hierarchical Clustering Heatmap of Top Differentially Expressed Genes\nDistinct Transcriptional Program of Antiviral Host Factors (Mock vs Influenza B)",
             fontsize=10, fontweight="bold", color="#0F2C59")

plt.tight_layout()
plt.savefig(figures_dir / "fig4_heatmap_degs.png", dpi=300)
plt.close()

# 5. Figure 5: Selected Antiviral and RNA Helicase Gene Comparison
selected_defense_genes = ["BATF2", "IFI35", "LY6E", "GBP4", "IFI27", "TLR3", "DTX3L", "DDX3X"]
fig, axes = plt.subplots(2, 4, figsize=(11, 5.5))
axes = axes.flatten()

for idx, gene in enumerate(selected_defense_genes):
    ax = axes[idx]
    g_data = expr_df[expr_df["gene_symbol"] == gene]
    if not g_data.empty:
        m_vals = g_data[mock_cols].values.flatten()
        i_vals = g_data[inf_cols].values.flatten()
        
        bp = ax.boxplot([m_vals, i_vals], patch_artist=True, widths=0.5, showfliers=False)
        bp["boxes"][0].set_facecolor("#94A3B8")
        bp["boxes"][1].set_facecolor("#0284C7" if gene != "DDX3X" else "#F59E0B")
        for box in bp["boxes"]:
            box.set_alpha(0.75)
            box.set_edgecolor("#1E293B")
        
        # Overlay individual sample points
        ax.scatter([1]*len(m_vals), m_vals, color="#334155", s=18, alpha=0.7, zorder=4)
        ax.scatter([2]*len(i_vals), i_vals, color="#0F172A", s=18, alpha=0.7, zorder=4)
        
        g_row = dge_df[dge_df["gene_symbol"] == gene].iloc[0]
        padj_str = f"padj = {g_row['padj']:.2e}" if g_row['padj'] < 0.001 else f"padj = {g_row['padj']:.3f}"
        fc_str = f"Log2FC = {g_row['log2FoldChange']:+.2f}"
        
        ax.set_title(f"{gene}\n({fc_str}, {padj_str})", fontsize=8, fontweight="bold",
                     color="#B45309" if gene == "DDX3X" else "#0F2C59")
        ax.set_xticklabels(["Mock", "Flu B"], fontsize=7.5)
        ax.set_ylabel("Normalized Expr", fontsize=7.5)
        ax.grid(True, linestyle=":", alpha=0.5)

plt.suptitle("Figure 5: Expression Level Comparison of Key Antiviral Defense Genes and DDX3X\nNormalized RNA-seq Counts Between Mock Controls (n=3) and Influenza B Explants (n=12)",
             fontsize=10, fontweight="bold", color="#0F2C59")
plt.tight_layout()
plt.savefig(figures_dir / "fig5_antiviral_genes_comparison.png", dpi=300)
plt.close()

# 6. Figure 6: Summary of Relevant Biological Functions
fig, ax = plt.subplots(figsize=(10, 5))
ax.axis("off")

categories = [
    ("Interferon-Stimulated Effector Defense", ["BATF2", "IFI35", "LY6E", "IFI27"],
     "Direct inhibition of viral replication, entry blockade, and immune potentiation.", "#0284C7"),
    ("Guanylate-Binding Protein GTPases", ["GBP2", "GBP4"],
     "Host cell-autonomous immunity, disruption of viral replication complexes.", "#059669"),
    ("Pattern Recognition & Ubiquitin Signaling", ["TLR3", "DTX3L", "TRANK1"],
     "Viral dsRNA sensing, downstream type I interferon induction, and epigenetic defense.", "#7C3AED"),
    ("Programmed Cell Death / Necroptosis", ["MLKL"],
     "Execution of necroptotic cell death to restrict viral propagation.", "#DC2626"),
    ("DEAD-Box RNA Helicase Core (DDX3X)", ["DDX3X"],
     "High constitutive expression (Log2 ≈ 17.0); minimal transcriptional induction (Log2FC = +0.15). Functions predominantly via post-translational activation, cellular localization, and direct RNA interaction rather than dramatic mRNA transcriptional upregulation.", "#D97706")
]

y_pos = 0.90
for title, genes, desc, col in categories:
    rect = Rectangle((0.02, y_pos - 0.14), 0.96, 0.15, facecolor="#F8FAFC", edgecolor=col, linewidth=1.5)
    ax.add_patch(rect)
    ax.text(0.04, y_pos - 0.03, title, fontsize=9.5, fontweight="bold", color=col)
    ax.text(0.40, y_pos - 0.03, f"Key Genes: {', '.join(genes)}", fontsize=8.5, fontweight="bold", color="#1E293B")
    ax.text(0.04, y_pos - 0.09, desc, fontsize=8, color="#475569", wrap=True)
    y_pos -= 0.18

ax.text(0.5, 0.98, "Figure 6: Functional Overview of Host Response Pathways Modulated During Influenza B Infection",
        ha="center", fontsize=10.5, fontweight="bold", color="#0F2C59")

plt.tight_layout()
plt.savefig(figures_dir / "fig6_biological_functions_summary.png", dpi=300)
plt.close()
