## 1 Data Requirements  

### 1.1 Input Data Format  
The tool will accept RNA-Seq **count matrices** in the following formats:  
- **CSV/TSV**: Gene expression raw counts with genes in rows and samples in columns.  
- **TXT (HTSeq Counts)**: Text-based format where each row contains a gene ID and raw read count.  

#### **Example CSV Format:**  
| Gene ID | Sample_1 | Sample_2 | Sample_3 | ... |  
|---------|----------|----------|----------|-----|  
| BRCA1   | 1034     | 2043     | 3098     | ... |  
| TP53    | 589      | 679      | 820      | ... |  

### 1.2 Preprocessing Requirements  
Before differential expression analysis, the following preprocessing steps are required:  
1. **Remove genes with low counts** (e.g., genes with counts <10 in all samples).  
2. **Normalize counts** using DESeq2’s variance stabilizing transformation (VST).  
3. **Perform quality control (QC)** to filter out poor-quality samples.  

### 1.3 Expected Dataset Size  
- **Number of genes**: ~20,000  
- **Number of samples**: 50–500 (depending on dataset)  
- **File size**: Typically 10MB–500MB  

### 1.4 Data Sources  
| Dataset | Description | URL |  
|---------|-------------|-----|  
| **GSE102073** | Head and Neck Cancer RNA-Seq (GEO) | [GSE102073](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102073) |  
| **TCGA-HNSC** | Tumor vs. Normal RNA-Seq Data (TCGA) | [TCGA-HNSC](https://portal.gdc.cancer.gov/projects/TCGA-HNSC) |  
