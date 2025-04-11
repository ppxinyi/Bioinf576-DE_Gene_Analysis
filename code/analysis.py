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

def differential_expression(expression_df: pd.DataFrame, group_labels: pd.Series, method="ttest", log2fc_thresh=1, pval_thresh=0.05):
    """
    Perform differential expression analysis using log2 fold change and statistical test.

    Args:
        expression_df (pd.DataFrame): Log-transformed expression matrix (genes x samples).
        group_labels (pd.Series): Series indicating group membership for each sample.
        method (str): 'ttest' or 'wilcoxon'.
        log2fc_thresh (float): Minimum absolute log2 fold change threshold.
        pval_thresh (float): Maximum p-value threshold.
    
    Returns:
        pd.DataFrame: Differential expression results with log2FC and p-values.
    """
    group1 = expression_df.loc[:, group_labels == group_labels.unique()[0]]
    group2 = expression_df.loc[:, group_labels == group_labels.unique()[1]]

    results = []
    for gene in expression_df.index:
        g1_vals = group1.loc[gene]
        g2_vals = group2.loc[gene]

        log2fc = np.log2(g2_vals.mean() + 1) - np.log2(g1_vals.mean() + 1)

        if method == "ttest":
            stat, pval = ttest_ind(g1_vals, g2_vals, equal_var=False)
        elif method == "wilcoxon":
            stat, pval = ranksums(g1_vals, g2_vals)
        else:
            raise ValueError("Unsupported method. Use 'ttest' or 'wilcoxon'.")

        results.append((gene, log2fc, pval))

    res_df = pd.DataFrame(results, columns=["gene", "log2FC", "pval"]).set_index("gene")
    res_df["adj_pval"] = res_df["pval"] * len(res_df)  # Bonferroni correction

    # Filter significant genes
    sig_df = res_df[(res_df["adj_pval"] < pval_thresh) & (abs(res_df["log2FC"]) > log2fc_thresh)].sort_values("adj_pval")
    
    return sig_df
