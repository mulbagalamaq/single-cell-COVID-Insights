# scRNA-seq Analysis of COVID-19 PBMC Samples

This repository contains my analysis of single-cell RNA sequencing (scRNA-seq) data from Peripheral Blood Mononuclear Cells (PBMCs) derived from three COVID-19 patients and three healthy control individuals. Using **Scanpy**, a powerful tool for large-scale scRNA-seq analysis, I explored the cellular and molecular landscape to better understand the immune response to COVID-19.

## Project Overview

The goal of this project was to identify differential gene expression patterns, cluster cell types, and visualize key immune features. Each sample was subsampled to 1,500 cells to ensure consistency and facilitate meaningful comparisons between COVID-19 and healthy PBMC samples.

---

## Data Summary

- **Samples**: 6 PBMC samples (3 from COVID-19 patients, 3 from healthy controls)
- **Data Format**: Raw data in `.h5` format
- **Subsampling**: Each sample was subsampled to 1,500 cells to balance the analysis.
- **Focus**: Immune-related gene expression patterns and pathway enrichment.

---

## Analysis Workflow

### 1. Data Loading
- Imported data using `scanpy.read_10x_h5`.
- Ensured variable names were unique across all samples.

### 2. Quality Control (QC)
- Filtered cells based on mitochondrial, ribosomal, and hemoglobin gene expression.
- Removed low-quality cells and potential artifacts.

### 3. Normalization and Scaling
- Normalized expression values to a depth of 10,000 reads per cell.
- Log-transformed the data for downstream analysis.

### 4. Doublet Detection
- Identified potential doublets using **Scrublet** and flagged them for filtering.

### 5. Highly Variable Gene Selection
- Identified genes with the most variation across cells to focus the analysis.

### 6. Dimensionality Reduction
- Performed PCA for dimensionality reduction.
- Used t-SNE and UMAP for visualizing clusters.

### 7. Batch Correction
- Applied **Scanorama** to correct for batch effects and integrate data from all samples.

### 8. Clustering
- Clustered cells using **Leiden** and **Louvain** algorithms.
- Annotated clusters based on known marker genes.

### 9. Differential Expression Analysis
- Compared gene expression between COVID-19 and control samples to identify key marker genes.

### 10. Gene Set Enrichment Analysis (GSEA)
- Used **gseapy** to identify enriched pathways and biological processes in COVID-19 samples.

---

## Results Summary

### Clustering of PBMCs
- **UMAP and t-SNE Visualizations**: Showed distinct immune cell populations in COVID-19 and control samples.
- **Cell Type Identification**: Annotated clusters as T cells, B cells, NK cells, monocytes, and other immune cell types.

### Differential Gene Expression
- **Upregulated Genes**: Immune activation and inflammatory response genes were highly expressed in COVID-19 samples.
- **Downregulated Genes**: Genes related to immune regulation and cellular homeostasis were less expressed in COVID-19 samples.

### Enriched Pathways
- **Immune Response Pathways**: Interferon signaling, T-cell activation, and antigen presentation pathways were enriched.
- **Inflammation**: Increased activity in TNF and IL-6 signaling pathways highlighted the inflammatory state in COVID-19 samples.
- **Cellular Stress**: Pathways related to oxidative stress and apoptosis were more active in COVID-19 samples.

### Visualization of Key Genes
- **Dot Plots and Violin Plots**: Highlighted gene expression patterns across different cell types.
- **Heatmaps**: Displayed transcriptional differences between COVID-19 and control samples.

---

## Repository Structure

```plaintext
├── data/
│   ├── Raw data files (.h5 format)
├── results/
│   ├── Processed files (.h5ad format)
├── notebooks/
│   ├── Jupyter notebooks detailing each analysis step
├── scripts/
│   ├── Python scripts for batch processing and analysis automation
```

## Key Dependencies

This project relies on the following Python packages:

- **Scanpy**: Comprehensive toolkit for scRNA-seq data analysis.
- **Scrublet**: For detecting and filtering doublets.
- **Seaborn** and **Matplotlib**: For creating informative visualizations.
- **Pandas** and **NumPy**: For efficient data manipulation.
- **scikit-learn**: For additional clustering methods.
- **gseapy**: For performing Gene Set Enrichment Analysis (GSEA).

### Installation

To install the required dependencies, run:

```bash
pip install scanpy scrublet seaborn matplotlib pandas numpy scikit-learn gseapy
```

## Key Findings

This analysis revealed significant differences in immune cell composition and gene expression between COVID-19 and healthy PBMC samples, providing valuable insights into COVID-19’s impact on the immune system.

- **Distinct Clusters**: Immune cell types were clearly separated between conditions.
- **Pathway Insights**: Immune response pathways were enriched in COVID-19 samples, suggesting a heightened inflammatory state.
- **Key Genes**: Identified marker genes provided insights into immune activation and regulation.
