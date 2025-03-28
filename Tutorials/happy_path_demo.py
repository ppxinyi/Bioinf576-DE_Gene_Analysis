import pandas as pd
import numpy as np
from data_processing import normalize_counts, compute_z_scores, filter_low_variance_genes, log_transform
from analysis import run_deseq2
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

# Enable R <-> pandas conversion
pandas2ri.activate()


#load data
counts_file = "./data/example_counts.csv"
counts_df = pd.read_csv(counts_file, index_col=0)

# Normalize raw counts (TPM or FPKM, here using 'raw' method = TPM)
normalized_df = normalize_counts(counts_df, gene_lengths, method="raw")

# Filter low variance genes
filtered_df = filter_low_variance_genes(normalized_df, threshold=0.1)

# Log-transform expression matrix
log_df = log_transform(filtered_df)

results_df = run_deseq2(counts_df=counts_df, col_data=col_data, design_formula="~ condition")
results_df.head()

results_df.to_csv("deseq2_results.csv")
log_df.to_csv("log_transformed_expression.csv")
