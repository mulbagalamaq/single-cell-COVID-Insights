# scRNA-seq Analysis Using Scanpy for COVID-19 PBMC Samples

This project focuses on single-cell RNA sequencing (scRNA-seq) analysis of Peripheral Blood Mononuclear Cells (PBMCs) derived from three COVID-19 patients and three healthy control individuals. Using Scanpy, an efficient tool for large-scale scRNA-seq analysis, I analyzed the cellular landscape and gene expression patterns in these samples, each subsampled to 1,500 cells for consistency.

## Project Overview

I aim to highlight differential expression patterns, cluster cell types, and visualize key features within the immune cells of COVID-19 patients compared to healthy controls.

## Data Summary

- **Samples**: 6 PBMC samples (3 from COVID-19 patients, 3 from healthy controls)
- **Source**: Data prepared as .h5 files
- **Subsampling**: Each sample is subsampled to 1,500 cells to balance cell counts across samples.
- **Gene Expression Levels**: Processed to highlight key gene sets associated with immune response.

## Analysis Workflow

The workflow follows standard quality control, filtering, clustering, and visualization steps, customized for scRNA-seq data. Key steps include:

1. **Data Loading**: Load each sample using Scanpy’s `read_10x_h5` function and make variable names unique.
2. **Quality Control (QC)**: Perform quality control measures to filter cells based on mitochondrial, ribosomal, and hemoglobin gene expression.
3. **Normalization and Scaling**: Normalize each cell to a depth of 10,000 reads and log-transform the expression data.
4. **Doublet Detection**: Detect potential doublets using Scrublet and flag these cells for filtering.
5. **Highly Variable Gene Selection**: Identify highly variable genes for downstream clustering and visualization.
6. **Dimensionality Reduction**: Run PCA, t-SNE, and UMAP for visualization of clusters.
7. **Batch Correction**: Apply Scanorama for batch effect correction and integration across samples.
8. **Clustering**: Perform clustering using methods like Leiden and Louvain to identify distinct cell populations.
9. **Differential Expression Analysis**: Run differential expression analysis between COVID-19 and control samples to identify key marker genes.
10. **Gene Set Enrichment Analysis**: Use GSEA to identify enriched biological processes and pathways relevant to COVID-19.

## Repository Structure

The project files are organized as follows:

- `data/`: Contains raw data files (.h5 format) for each PBMC sample.
- `results/`: Stores output files, including filtered, normalized, and batch-corrected data in `.h5ad` format.
- `notebooks/`: Jupyter notebooks detailing each analysis step.
- `scripts/`: Contains Python scripts used for batch processing and analysis automation.

  
## Key Dependencies

This analysis relies on the following Python packages:

- **Scanpy**: For handling and analyzing single-cell data.
- **Scrublet**: For doublet detection.
- **Seaborn** and **Matplotlib**: For visualization.
- **Pandas** and **NumPy**: For data manipulation.
- **scikit-learn**: For additional clustering methods (e.g., KMeans).
- **gseapy**: For gene set enrichment analysis.

## Install the required packages using:

```bash

pip install scanpy scrublet seaborn matplotlib pandas numpy scikit-learn gseapy

```

## Results Summary

This analysis provides comprehensive insights into the immune cell composition and gene expression changes associated with COVID-19 infection in PBMC samples. Key findings and visualizations include:

### 1. Clustering of PBMCs
   - **UMAP and t-SNE Plots**: These dimensionality reduction techniques reveal distinct clusters within the PBMC samples, showing clear separation between COVID-19 and healthy control samples.
   - **Cell Type Identification**: Using known marker genes, major cell types such as T cells, B cells, NK cells, and monocytes are identified, highlighting differences in cell type distribution between COVID-19 and control groups.

### 2. Differential Gene Expression
   - **Upregulated Genes**: Several genes show increased expression in COVID-19 samples compared to controls, including genes associated with immune activation, inflammation, and cytokine response.
   - **Downregulated Genes**: Conversely, genes linked to normal immune regulation and cellular homeostasis tend to be downregulated in COVID-19 samples.
   - **Marker Gene Analysis**: Top marker genes are identified for each cluster, allowing for a detailed comparison of immune cell activation states in COVID-19 versus control samples.

### 3. Enriched Pathways
   - **Pathways Related to Immune Response**: Pathway analysis highlights significant enrichment in immune response pathways in COVID-19 samples, including interferon signaling, T-cell activation, and antigen presentation.
   - **Inflammatory Pathways**: Pathways associated with inflammation, such as TNF signaling and IL-6 signaling, are more active in COVID-19 samples, indicating a heightened inflammatory state.
   - **Cellular Stress and Apoptosis**: COVID-19 samples also show increased activity in pathways related to cellular stress, apoptosis, and oxidative stress, reflecting the immune response's impact on cell health.

### 4. Visualization of Key Genes
   - **Dot Plots and Violin Plots**: Expression levels of critical genes across cell types are visualized to compare COVID-19 and control samples. These plots allow for easy identification of gene expression trends.
   - **Heatmaps**: Heatmaps display differentially expressed genes, highlighting the unique transcriptional signatures of COVID-19 versus control PBMCs.

### 5. Gene Set Enrichment Analysis (GSEA)
   - **Significant Biological Processes**: GSEA reveals specific biological processes that are upregulated in COVID-19 samples, including immune cell activation, cytokine production, and response to viral infection.
   - **Top Enriched Gene Sets**: Gene sets from KEGG, GO Biological Process, and other databases provide insights into the molecular mechanisms involved in the immune response to COVID-19.

Saw distinct cellular and molecular differences between COVID-19 and healthy PBMC samples, contributing to a better understanding of COVID-19’s impact on the immune system.


## Acknowledgments

This project leverages several open-source tools and libraries that are widely used in the field of single-cell RNA sequencing analysis. Special thanks to the developers and contributors of these tools:

- **[Scanpy](https://github.com/theislab/scanpy)**: A scalable toolkit for analyzing single-cell gene expression data, developed by the Theis Lab at Helmholtz Zentrum München.
- **[Scrublet](https://github.com/AllonKleinLab/scrublet)**: A tool for detecting doublets in single-cell RNA-seq data, created by the Klein Lab at Harvard Medical School.
- **[gseapy](https://github.com/zqfang/GSEApy)**: A Python wrapper for Gene Set Enrichment Analysis, making it easier to perform enrichment analyses with popular databases.

I also acknowledge the broader single-cell and bioinformatics community for developing and sharing tools, resources, and knowledge that have made this analysis possible.

## Contact

For questions, issues, or contributions to the project, please feel free to reach out via email at [aymenmaqsood.2000@gmail.com](mailto:aymenmaqsood.2000@gmail.com). 

Thank you for your interest in this project!
