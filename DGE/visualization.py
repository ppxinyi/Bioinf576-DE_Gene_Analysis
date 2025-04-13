"""
Module: Visualization
Author: Xinyi Deng
Description: This module contains visualization functions for gene expression data.
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
import os
import numpy as np



def plot_heatmap(
    expression_df: pd.DataFrame,
    top_genes: list,
    metadata: pd.DataFrame = None,
    group_col: str = "integration",
    title: str = "DEG Heatmap"
):
    """
    Plot heatmap of top differentially expressed genes and auto-save to PDF.

    Args:
        expression_df (pd.DataFrame): Z-score normalized expression matrix.
        top_genes (list): List of top gene names to include in the heatmap.
        metadata (pd.DataFrame): Metadata DataFrame indexed by sample.
        group_col (str): Column in metadata to group samples (e.g., 'integration', 'fusion').
        title (str): Title for the heatmap.
    """
    # Subset expression data
    data = expression_df.loc[top_genes]

    # Set up column color bar
    col_colors = None
    if metadata is not None and group_col in metadata.columns:
        metadata = metadata.loc[data.columns]
        group_values = metadata[group_col].unique()
        palette = sns.color_palette("Set2", len(group_values))
        color_map = dict(zip(group_values, palette))
        col_colors = metadata[group_col].map(color_map)

    # Plot clustermap
    g = sns.clustermap(
        data,
        cmap="vlag",
        col_cluster=False,
        row_cluster=True,
        col_colors=col_colors,
        figsize=(12, max(6, len(top_genes) * 0.3)),
        cbar_kws={"label": "Z-score"}
    )

    plt.suptitle(title, y=1.05, fontsize=14)

    # Add legend
    if col_colors is not None:
        for label, color in color_map.items():
            g.ax_col_dendrogram.bar(0, 0, color=color, label=label, linewidth=0)
        g.ax_col_dendrogram.legend(loc="center", ncol=len(color_map), bbox_to_anchor=(0.5, 1.1), frameon=False)

    # Auto-save figure to PDF with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"heatmap_{group_col}_{timestamp}.pdf"
    g.savefig(filename, format="pdf")
    print(f"✅ Heatmap automatically saved to: {filename}")

    plt.show()
    
def plot_volcano(deg_df, log2fc_thresh=1, pval_thresh=0.05, title="Volcano Plot"):
    deg_df = deg_df.copy()
    deg_df["-log10(pval)"] = -np.log10(deg_df["adj_pval"])
    deg_df["significant"] = (deg_df["adj_pval"] < pval_thresh) & (abs(deg_df["log2FC"]) > log2fc_thresh)

    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=deg_df, x="log2FC", y="-log10(pval)", hue="significant",
                    palette={True: "red", False: "grey"}, edgecolor=None, s=25)
    plt.axvline(x=log2fc_thresh, color="black", linestyle="--", linewidth=1)
    plt.axvline(x=-log2fc_thresh, color="black", linestyle="--", linewidth=1)
    plt.axhline(y=-np.log10(pval_thresh), color="black", linestyle="--", linewidth=1)
    plt.title(title)
    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10 Adjusted P-value")
    plt.tight_layout()
    filename = f"volcano_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf"
    plt.savefig(filename)
    print(f"✅ Volcano plot saved as {filename}")
    plt.show()


def plot_pca(expression_df, sample_info, group_col="integration", title="PCA Plot"):
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(expression_df.T)
    pca_df = pd.DataFrame(pcs, columns=["PC1", "PC2"], index=expression_df.columns)
    pca_df[group_col] = sample_info.set_index("Sample").loc[pca_df.index, group_col]

    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue=group_col, palette="Set2", s=80)
    plt.title(title)
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.tight_layout()
    filename = f"pca_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf"
    plt.savefig(filename)
    print(f"✅ PCA plot saved as {filename}")
    plt.show()


def plot_gene_boxplot(expression_df, gene_name, sample_info, group_col="integration", title=None):
    df = pd.DataFrame({
        "expression": expression_df.loc[gene_name],
        group_col: sample_info.set_index("Sample").loc[expression_df.columns, group_col]
    }).reset_index(drop=True)

    plt.figure(figsize=(6, 4))
    sns.boxplot(data=df, x=group_col, y="expression", palette="Set2")
    sns.stripplot(data=df, x=group_col, y="expression", color="black", size=4, jitter=0.15)
    plt.title(title if title else f"Expression of {gene_name}")
    plt.tight_layout()
    filename = f"boxplot_{gene_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf"
    plt.savefig(filename)
    print(f"✅ Boxplot saved as {filename}")
    plt.show()
