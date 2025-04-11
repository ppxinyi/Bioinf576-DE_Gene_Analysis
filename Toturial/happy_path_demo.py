# Happy Path RNA-seq DEG Analysis Tutorial

### 🧪 Goal of the Tutorial
In this tutorial, we demonstrate a complete RNA-seq data analysis pipeline to identify differentially expressed genes (DEGs) across sample groups. We aim to answer the biological question:

> *Are there significantly differentially expressed genes between tumor samples with and without gene fusion events?*

---

### 📊 Dataset Description
We use two input files:
- `expression_matrix.csv`: Gene expression matrix where rows are genes and columns are sample names. This is either raw count data or normalized expression values.
- `sample_info.csv`: Metadata table containing sample labels and grouping variables (e.g., `fusion`, `integration`, `type`). It must contain a `Sample` column that matches the expression matrix columns.

---

### 🔧 Step 1: Environment Setup

```python
import pandas as pd
from data_processing import load_data, log_transform, compute_z_scores, filter_low_variance_genes
from analysis import differential_expression
from visual import plot_heatmap, plot_volcano, plot_pca, plot_gene_boxplot
```

---

### 📥 Step 2: Load Input Data

```python
expression_df, sample_info = load_data("expression_matrix.csv", "sample_info.csv")
expression_df.shape, sample_info.shape
```

---

### 🔬 Step 3: Data Preprocessing

```python
log_expr = log_transform(expression_df)
z_expr = compute_z_scores(log_expr)
z_expr = filter_low_variance_genes(z_expr)
```

---

### 📐 Step 4: Differential Expression Analysis

```python
group_col = "fusion"  # or "integration", "type"
group_labels = sample_info[group_col]
method = suggest_test_method(group_labels)
deg_df = differential_expression(log_expr, group_labels, method=method)
deg_df.head()
```

---

### 🖼️ Step 5: Visualization

#### 🔥 Top DEG Heatmap
```python
top_genes = deg_df.head(20).index.tolist()
plot_heatmap(z_expr, top_genes, metadata=sample_info, group_col=group_col)
```

#### 🌋 Volcano Plot
```python
plot_volcano(deg_df, title=f"Volcano Plot - {group_col}")
```

#### 🧭 PCA Plot
```python
plot_pca(z_expr, sample_info, group_col=group_col)
```

#### 📦 Boxplot of Top Gene
```python
plot_gene_boxplot(log_expr, gene_name=deg_df.index[0], sample_info=sample_info, group_col=group_col)
```

---

### 📊 Step 6: Result Summary
- Total number of DEGs: `len(deg_df)`
- Top genes include: `deg_df.head(5).index.tolist()`
- Result saved as: `DEG_fusion_YYYYMMDD.csv`

---

### 🧬 Biological Interpretation
Our results show that there are several genes significantly associated with fusion status. For example, `TP53`, `MMP9`, and `CXCL10` are upregulated in fusion-positive samples.

This pipeline enables the identification of potential biomarkers and pathways involved in tumor subtypes, supporting future biological discovery and clinical research.

---

### ✅ Conclusion
This notebook demonstrates how to:
- Run differential expression testing
- Visualize and interpret the output
