"""
Module: Data Processing
Author: Xinyi Deng
Description: This module contains functions for processing RNA-Seq data.
"""
import numpy as np
import pandas as pd
def normalize_counts(counts_df: pd.DataFrame, method="TPM") -> pd.DataFrame:
    """
    Normalize raw count data using the specified method.
    
    Args:
        counts_df (pd.DataFrame): DataFrame with raw count data (genes as rows, samples as columns).
        method (str): Normalization method, either 'TPM' or 'FPKM'.
    
    Returns:
        pd.DataFrame: Normalized expression values.
    """
