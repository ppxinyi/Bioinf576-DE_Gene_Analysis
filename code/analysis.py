"""
Module: Run Deseq2 
Author: Xinyi Deng
"""
import rpy2.robjects as ro

def run_deseq2(counts_df, col_data, design_formula="~ condition"):
    """
    Perform differential expression analysis using DESeq2 in R.

    Args:
        counts_df (pd.DataFrame): DataFrame with raw count data (genes as rows, samples as columns).
        col_data (pd.DataFrame): DataFrame with sample information, including condition information.
        design_formula (str): The design formula for the DESeq2 analysis (default: "~ condition").

    Returns:
        pd.DataFrame: DataFrame with DESeq2 results (log fold change, p-values, etc.).
    """
    # Convert the DataFrames to R objects
    counts_r = ro.r['as.data.frame'](counts_df)
    col_data_r = ro.r['as.data.frame'](col_data)

    # Import DESeq2 package
    ro.r('library(DESeq2)')

    # Create DESeqDataSet
    ro.r('dds <- DESeqDataSetFromMatrix(countData=counts_r, colData=col_data_r, design=' + design_formula + ')')

    # Run DESeq2 analysis
    ro.r('dds <- DESeq(dds)')

    # Get the results
    ro.r('res <- results(dds)')

    # Return the results as a DataFrame
    res_df = ro.r('as.data.frame(res)')
    return pd.DataFrame(res_df)
