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
    """
    return np.log2(expression_df + 1)
  def detect_outliers(expression_df: pd.DataFrame, threshold: float = 3.0) -> pd.DataFrame:
        """
        Identify outlier samples based on Z-score.
        """
     z_scores = (expression_df - expression_df.mean()) / expression_df.std()
     return (z_scores.abs() > threshold)

