"""
Module: Analyse DE Genes
Author: Xinyi Deng
"""


from scipy.stats import ttest_ind, ranksums, f_oneway, kruskal
import pandas as pd
import numpy as np


def suggest_test_method(group_labels: pd.Series, verbose=True) -> str:
    """
    Suggest a statistical test method based on the number of groups and sample size.

    Args:
        group_labels (pd.Series): Series indicating group membership.
        verbose (bool): If True, print recommendation reason.

    Returns:
        str: Suggested method ('ttest', 'wilcoxon', 'anova', or 'kruskal').
    """
    unique_groups = group_labels.unique()
    num_groups = len(unique_groups)
    group_sizes = group_labels.value_counts()

    if num_groups == 2:
        n1, n2 = group_sizes.iloc[0], group_sizes.iloc[1]
        if n1 >= 20 and n2 >= 20:
            method = "ttest"
            reason = "Two groups with sufficient sample size (>=20) per group → use t-test."
        else:
            method = "wilcoxon"
            reason = "Two groups with small sample sizes → use Wilcoxon rank-sum test."

    elif num_groups > 2:
        if all(n >= 20 for n in group_sizes):
            method = "anova"
            reason = "More than two groups with decent sample sizes → use ANOVA."
        else:
            method = "kruskal"
            reason = "More than two groups with small samples → use Kruskal-Wallis test."
    else:
        raise ValueError("Invalid number of groups for comparison.")

    if verbose:
        print(f"✅ Suggested method: {method.upper()}")
        print(f" Reason: {reason}")

    return method
    
def differential_expression(expression_df: pd.DataFrame, group_labels: pd.Series, method="ttest", log2fc_thresh=1, pval_thresh=0.05):
    """
    Perform differential expression analysis using appropriate statistical test.

    Args:
        expression_df (pd.DataFrame): Log-transformed expression matrix (genes x samples).
        group_labels (pd.Series): Series indicating group membership for each sample.
        method (str): 'ttest', 'wilcoxon' (for 2 groups) or 'anova', 'kruskal' (for >2 groups).
        log2fc_thresh (float): Absolute log2FC threshold for filtering (only used for 2-group).
        pval_thresh (float): Adjusted p-value threshold.

    Returns:
        pd.DataFrame: Differential expression results.
    """
    group_labels = group_labels.loc[expression_df.columns]  # align index
    unique_groups = group_labels.unique()
    num_groups = len(unique_groups)

    results = []

    for gene in expression_df.index:
        values = [expression_df.loc[gene, group_labels == grp].values.flatten() for grp in unique_groups]

        if num_groups == 2:
            g1, g2 = values
    
            # remove NaN
            g1 = g1[~np.isnan(g1)]
            g2 = g2[~np.isnan(g2)]
            if len(g1) == 0 or len(g2) == 0:
                continue
            log2fc = np.log2(g2.mean() + 1) - np.log2(g1.mean() + 1)

            if method == "ttest":
                stat, pval = ttest_ind(g1, g2, equal_var=False)
            elif method == "wilcoxon":
                stat, pval = ranksums(g1, g2)
            else:
                raise ValueError("Unsupported 2-group test. Use 'ttest' or 'wilcoxon'.")
            results.append((gene, log2fc, pval))

        else:
            values = [v[~np.isnan(v)] for v in values]

            # skip any empty group
            if any(len(v) == 0 for v in values):
                continue

            # skip if all values are the same (kruskal/f_oneway will crash)
            all_values = np.concatenate(values)
            if np.all(all_values == all_values[0]):
                continue
            if method == "anova":
                stat, pval = f_oneway(*values)
            elif method == "kruskal":
                stat, pval = kruskal(*values)
            else:
                raise ValueError("Unsupported multi-group test. Use 'anova' or 'kruskal'.")
            results.append((gene, np.nan, pval))  # no log2FC for multi-group

    res_df = pd.DataFrame(results, columns=["gene", "log2FC", "pval"]).set_index("gene")
    res_df["adj_pval"] = res_df["pval"] * len(res_df)

    if num_groups == 2:
        sig_df = res_df[(res_df["adj_pval"] < pval_thresh) & (abs(res_df["log2FC"]) > log2fc_thresh)]
    else:
        sig_df = res_df[res_df["adj_pval"] < pval_thresh]

    return sig_df.sort_values("adj_pval")

