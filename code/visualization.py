"""
Module: Visualization
Author: Xinyi Deng
Description: This module contains visualization functions for gene expression data.
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

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
