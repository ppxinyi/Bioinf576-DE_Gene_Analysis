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
  - Apply some other normalized method(like 
- **2.3 Batch Effect Correction**
  - Integrate Combat-Seq for batch effect correction

## 3. PCA & Differential Expression Analysis (Week 3)
- **3.1 Principal Component Analysis (PCA)**
  - Compute and visualize PCA for quality assessment.
  - if smaples doen't clustered by conditions, find some other function to fix batch effect.
    
## 4. Differential Expression Analysis (Week 5)
- **4.1 Differential Expression Analysis**
  - choose suitable statistic method to calculate the expression counts between different conditions
  - Generate results as a table with adjusted p-values and log fold changes.
  - Using p-value combined with log2 fold changes to select significant genes

#      -------Above finished by 03/28/2025--------------------  #

## 5. Interactive Visualization (Week 6)
- **5.1 Gene Expression Heatmap**
  - Create an interactive heatmap using Plotly/Seaborn.
- **5.2 Volcano Plot**
  - Visualize significantly differentially expressed genes.

## 6. Gene notation (Week 6)
- **6.1 Gene ID Conversion**
  - Convert gene IDs using databases like org.Hs.eg.db.
  - Ensure consistency of gene identifiers.
  - Annotate genes with GO terms and pathway information.

## 7. Testing & Documentation (Week 7)
- **7.1 Testing**
  - Conduct unit testing and user testing.
- **7.2 Documentation**
  - Provide user guides and API documentation.
    
## 8. Compare with GSEA Analysis result (Week 7)
- **8.1 GSEA Analysis Comparison**
  -  based on the DEG table, add functional enrichment analysis function.
  - Compare differentially expressed genes with enriched pathways from GSEA.

## 9. Final Review & Launch (Week 8)
- **9.1 Review & Optimization**
  - Optimize performance and fix bugs.
- **9.2 Launch & User Feedback Collection**
  - Release the tool and gather feedback for improvements.
    


