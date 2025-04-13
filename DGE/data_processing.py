"""
Module: Data Processing and preprocessing
Author: Xinyi Deng
Description: This module contains functions for processing RNA-Seq data.
"""

import numpy as np
import pandas as pd

def normalize_counts(counts_df: pd.DataFrame, gene_lengths: pd.Series, method="raw", min_expression=10) -> pd.DataFrame:
    """
    Normalize raw count data using the specified method.
    
    Args:
        counts_df (pd.DataFrame): DataFrame with raw count data (genes as rows, samples as columns).
        gene_lengths (pd.Series): Series with gene lengths (index: gene names, values: lengths).
        method (str): Normalization method, either 'TPM' or 'FPKM'.
        min_expression (int): Minimum expression threshold to filter low-expression genes.
    
    Returns:
        pd.DataFrame: Normalized expression values.
    """
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
        raise ValueError("Unsupported normalization method. Choose 'raw' or 'FPKM'.")

    return normalized_df

def compute_z_scores(expression_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute Z-scores for gene expression data.

    Args:
        expression_df (pd.DataFrame): DataFrame with expression values (genes as rows, samples as columns).
    
    Returns:
        pd.DataFrame: DataFrame with Z-score normalized values.
    """
    return (expression_df - expression_df.mean(axis=1, keepdims=True)) / expression_df.std(axis=1, keepdims=True)

def filter_low_variance_genes(expression_df: pd.DataFrame, threshold: float = 0.1) -> pd.DataFrame:
    """
    Remove genes with low variance across samples.

    Args:
        expression_df (pd.DataFrame): Normalized gene expression data.
        threshold (float): Minimum variance required to retain a gene.
    
    Returns:
        pd.DataFrame: Filtered expression matrix with only high-variance genes.
    """
    variances = expression_df.var(axis=1)
    return expression_df.loc[variances > threshold]

def log_transform(expression_df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply log2 transformation to gene expression data to reduce skewness.

    Args:
        expression_df (pd.DataFrame): Normalized gene expression data.
    
    Returns:
        pd.DataFrame: Log2 transformed expression data.
    """
    return np.log2(expression_df + 1)
def load_data(expression_path: str, sample_info_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load expression matrix and sample metadata from CSV files.

    Args:
        expression_path (str): Path to the expression matrix CSV.
        sample_info_path (str): Path to the sample info CSV.

    Returns:
        tuple: (expression_df, sample_info_df) both as pandas DataFrames.
    """
    expression_df = pd.read_csv(expression_path, index_col=0)
    sample_info_df = pd.read_csv(sample_info_path)
    sample_info_df = sample_info_df.set_index("Sample")
    # Align columns
    expression_df = expression_df.loc[:, sample_info_df.index]
    return expression_df, sample_info_df
