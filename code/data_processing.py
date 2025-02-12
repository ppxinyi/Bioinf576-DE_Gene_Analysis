"""
Module: Data Processing
Author: Xinyi Deng
Description: This module contains functions for processing RNA-Seq data.
"""
import numpy as np
import pandas as pd
def normalize_counts(counts_df: pd.DataFrame, method="raw") -> pd.DataFrame:
    """
    Normalize raw count data using the specified method.
    
    Args:
        counts_df (pd.DataFrame): DataFrame with raw count data (genes as rows, samples as columns).
        method (str): Normalization method, either 'TPM' or 'FPKM'.
    
    Returns:
        pd.DataFrame: Normalized expression values.
    """
    counts_df = pd.read_csv(file_path)


# Filter out low-expression genes
    counts_df = counts_df[counts_df.sum(axis=1) >= min_expression]

    if method == "raw":
        # Compute per-sample scaling factors
        rpk = counts_df.div(gene_lengths, axis=0) * 1e3
        scaling_factors = rpk.sum(axis=0) / 1e6
        normalized_df = rpk.div(scaling_factors, axis=1)
    elif method == "FPKM":
        rpk = counts_df.div(gene_lengths, axis=0) * 1e3
        total_counts = counts_df.sum(axis=0)
        normalized_df = rpk.div(total_counts, axis=1) * 1e6
    else:
        raise ValueError("Unsupported normalization method. Choose 'TPM' or 'FPKM'.")

    return normalized_df
