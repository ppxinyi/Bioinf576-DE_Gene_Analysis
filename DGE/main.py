"""
Script: main.py
Author: Xinyi Deng
Description: Full RNA-seq DEG pipeline with command-line support using argparse.
"""

import argparse
import pandas as pd
import numpy as np
from datetime import datetime
from .data_processing  import (
    load_data,
    normalize_counts,
    log_transform,
    compute_z_scores,
    filter_low_variance_genes,
)
from .analysis import differential_expression, suggest_test_method
from .visualization import (
    plot_heatmap,
    plot_volcano,
    plot_pca,
    plot_gene_boxplot,
)

def main():
    # === Command-line argument parser ===
    parser = argparse.ArgumentParser(description="Run full RNA-seq DEG pipeline with visualizations.")
    parser.add_argument("--expression", required=True, help="Path to expression matrix CSV file.")
    parser.add_argument("--sample_info", required=True, help="Path to sample metadata CSV file.")
    parser.add_argument("--group_col", default="integration", help="Column in metadata to group by (e.g., integration, fusion).")
    parser.add_argument("--data_type", choices=["raw", "normalized"], default="raw", help="Specify whether expression data is raw counts or already normalized.")
    parser.add_argument("--method", default=None, help="Statistical test to use. Options: ttest, wilcoxon, anova, kruskal. Leave blank to auto-select.")

    args = parser.parse_args()
    
    # === Load expression matrix and sample metadata ===
    expression_df, sample_info = load_data(args.expression, args.sample_info)
    print(f"✅ Loaded {expression_df.shape[0]} genes and {expression_df.shape[1]} samples.")
    
    # === Preprocessing: normalization & transformation ===
    if args.data_type == "raw":
        print("🛠️  Normalizing raw counts and applying log2 transformation...")
        # Optional: use gene lengths if needed for TPM/FPKM
        # gene_lengths = pd.read_csv("gene_lengths.csv", index_col=0).squeeze()
        # normalized = normalize_counts(expression_df, gene_lengths, method="FPKM")
        normalized = expression_df
        log_expr = log_transform(normalized)
    else:
        print("🔁 Using pre-normalized expression matrix (assumed to be log2-transformed).")
        log_expr = expression_df
    
    # === Z-score normalization and low-variance gene filtering ===
    z_expr = compute_z_scores(log_expr)
    z_expr = filter_low_variance_genes(z_expr)
    print(f" Retained {z_expr.shape[0]} genes after variance filtering.")
    
    # === Differential expression analysis ===
    group_labels = sample_info[args.group_col]
    if args.method is None:
        method = suggest_test_method(group_labels)
    else:
        method = args.method
    deg_df, full_df = differential_expression(log_expr, group_labels, method=method)
    
    # === Save DEG result table ===
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    deg_filename = f"DEG_{args.group_col}_{timestamp}.csv"
    deg_df.reset_index().to_csv(deg_filename, index=False)
    print(f"📄 DEG results saved to: {deg_filename}")
    full_filename = f"DEG_full_{args.group_col}_{timestamp}.csv"
    full_df.reset_index().to_csv(full_filename, index=False)
    print(f"🧾 Full DEG result saved to: {full_filename}")
    # === Visualization ===
    if deg_df.shape[0] > 0:
        top_genes = deg_df.index.tolist()
    
        print("Generating heatmap...")
        plot_heatmap(expression_df, top_genes, metadata=sample_info, group_col=args.group_col)
    
        print("Generating volcano plot...")
        plot_volcano(deg_df, title=f"Volcano Plot - {args.group_col}")
    
        print("Generating PCA plot...")
        plot_pca(z_expr, sample_info, group_col=args.group_col)
    
        print(f"Generating boxplot for top DEG: {deg_df.index[0]}")
        for gene in top_genes:
            print(f"Plotting boxplot for: {gene}")
            plot_gene_boxplot(expression_df, gene_name=gene, sample_info=sample_info, group_col=args.group_col)
    else:
        print("⚠️ No significant DEGs found. Using top 20 genes by lowest adjusted p-value instead.")
        top_genes = full_df.sort_values("adj_pval").head(20).index.tolist()
        print("Generating heatmap...")
        plot_heatmap(expression_df, top_genes, metadata=sample_info, group_col=args.group_col)
    
        print("Generating volcano plot...")
        plot_volcano(full_df, title=f"Volcano Plot - {args.group_col}")
    
        print("Generating PCA plot...")
        plot_pca(z_expr, sample_info, group_col=args.group_col)
        
        print(f"Generating boxplot for top 5 Gene")
        print("⚠️ No genes available for boxplot. Using top 5 genes by lowest adjusted p-value instead.")
        top5_genes = full_df.sort_values("adj_pval").head(5).index.tolist()
        for gene in top5_genes:
            print(f"Plotting boxplot for: {gene}")
            plot_gene_boxplot(expression_df, gene_name=gene, sample_info=sample_info, group_col=args.group_col)

if __name__ == "__main__":
    main()
