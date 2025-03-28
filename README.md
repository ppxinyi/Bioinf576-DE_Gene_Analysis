# Bioinf576-Project

Objective: Provide a one-stop, intelligent RNA-Seq differential expression analysis tool for non-programming users.

-Input: Raw gene expression data from Bulk RNA Seq.

-Output: Significant different gene between different conditions calculation report and the gene expression heatmap. 

Project Lifecycle
1. Initiation
•	Objective: Create an intelligent RNA-Seq analysis tool for non-programmers.
•	Stakeholders: Researchers, bioinformaticians, developers.
•	Technology Stack: Python (backend), R (DESeq2/edgeR), Streamlit/Flask (frontend), Plotly/Seaborn (visualization).
2. Planning
•	Milestones:
1.	Data upload and QC.
2.	Normalization and batch correction.
3.	PCA and differential analysis.
4.	Interactive visualizations (heatmaps, volcano plots).
1 weeks per milestone.

## Toturil 
You can view the example demo of how to use this tool with sample data in this Jupyter notebook:

📎 [Happy Path Demo Notebook](./tutorials/happy_path_demo.ipynb)

This notebook showcases the entire pipeline — from loading raw counts and metadata to normalization, filtering, and running DESeq2 to obtain differentially expressed genes.
