## 1. Data Upload & Quality Control (Week 1)
- **1.1 Data Upload**
  - Support CSV, TSV, and Excel formats.
  - need raw count matrix
- **1.2 Quality Control (QC)**
  - Compute basic statistics (e.g., missing values,distribution, outliers).
  - Generate QC visualizations (box plots, histograms).

## 2. Normalization & Batch Correction (Week 2)
- **2.1 Data Preprocessing**
  - Convert raw counts to TPM/FPKM (if necessary).
  - Perform log transformation.
- **2.2 Normalization Methods**
  - Apply DESeq2's variance stabilizing transformation (VST) or regularized log transformation (rlog).
- **2.3 Batch Effect Correction**
  - Integrate Combat-Seq for batch effect correction.

## 3. PCA & Differential Expression Analysis (Week 3)
- **3.1 Principal Component Analysis (PCA)**
  - Compute and visualize PCA for quality assessment.
- **3.2 Differential Expression Analysis**
  - choose suitable statistic method to calculate the expression counts between different conditions
  - Generate results as a table with adjusted p-values and log fold changes.
  - Using p-value combined with log2 fold changes to select significant genes

## 4. Interactive Visualization (Week 4)
- **4.1 Gene Expression Heatmap**
  - Create an interactive heatmap using Plotly/Seaborn.
- **4.2 Volcano Plot**
  - Visualize significantly differentially expressed genes.


## 5. Testing & Documentation (Week 5)
- **5.1 Testing**
  - Conduct unit testing and user testing.
- **6.2 Documentation**
  - Provide user guides and API documentation.

## 6. Final Review & Launch (Week 6)
- **6.1 Review & Optimization**
  - Optimize performance and fix bugs.
- **6.2 Launch & User Feedback Collection**
  - Release the tool and gather feedback for improvements.
