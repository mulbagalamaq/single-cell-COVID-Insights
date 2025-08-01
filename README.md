# Single-Cell RNA-seq Analysis of COVID-19 PBMC Samples

## Overview

This repository contains my analysis of single-cell RNA sequencing data from Peripheral Blood Mononuclear Cells (PBMCs) comparing COVID-19 patients to healthy controls. I used Scanpy to process the data and identify key differences in immune cell gene expression patterns.

## What I Found

### Key Results
- **5,312 high-quality cells** analyzed after quality control filtering
- **1,247 differentially expressed genes** between COVID-19 and control samples
- **23 significantly enriched pathways** related to immune response
- **85% statistical power** achieved with current sample size

### Main Biological Insights
1. **Strong interferon response** in COVID-19 patients (IFIT1, IFIT3, ISG15 upregulated)
2. **T-cell exhaustion** signatures (CD3D, CD4, CD8A downregulated)
3. **Enhanced monocyte activation** and inflammatory pathways
4. **Batch effects detected** but controlled for in analysis

## Analysis Pipeline

I followed a standard scRNA-seq workflow with some enhancements:

```
Raw Data → Quality Control → Normalization → 
Dimensionality Reduction → Clustering → 
Differential Expression → Pathway Analysis
```

### Quality Control
- Filtered cells with <200 genes or <1000 UMIs
- Removed cells with >20% mitochondrial content
- Detected and removed statistical outliers using IQR method
- Assessed batch effects using ANOVA (F=52.7, p<0.001)

### Statistical Methods
- **Differential Expression**: Wilcoxon rank-sum test with Benjamini-Hochberg correction
- **Effect Sizes**: Calculated Cohen's d for all significant genes
- **Power Analysis**: 85% power achieved with current sample size
- **Pathway Analysis**: Enrichr across 5 databases (KEGG, GO, Reactome, Hallmark)

## Visualizations

### UMAP Clustering
![UMAP by Condition](plots/umap_by_condition.png)
*Clear separation between COVID-19 and control samples*

![UMAP by Clusters](plots/umap_by_clusters.png)
*18 distinct cell populations identified by Leiden clustering*

### Quality Control
![Quality Control Metrics](plots/quality_control_metrics.png)
*QC metrics across conditions - genes per cell, UMI counts, mitochondrial content*

### Differential Expression
![Volcano Plot](plots/volcano_plot.png)
*Volcano plot showing significantly upregulated (red) and downregulated (blue) genes*

![Marker Genes Heatmap](plots/marker_genes_heatmap.png)
*Heatmap of top marker genes by condition*

### Pathway Analysis
![KEGG Pathways](plots/kegg_pathways.png)
*Top enriched KEGG pathways - interferon signaling most significant*

![Hallmark Pathways](plots/hallmark_pathways.png)
*Enriched Hallmark gene sets showing immune response patterns*

### Summary Statistics
![Summary Statistics](plots/summary_statistics.png)
*Overview of QC metrics, DE results, pathway enrichment, and performance*

## Technical Details

### Data
- **6 PBMC samples**: 3 COVID-19, 3 healthy controls
- **10x Genomics Chromium** sequencing
- **9,000 cells initially**, 5,312 after QC
- **17,993 genes** detected after filtering

### Software Stack
```
Python 3.10
Scanpy 1.9.0
Pandas 2.2.0
NumPy 1.26.0
Seaborn 0.11.0
GSEApy 0.12.0
Statsmodels 0.14.0
```

### Performance
- **Processing time**: ~1 minute
- **Memory usage**: 16 GB RAM
- **Statistical power**: 85%

## Key Findings

### Upregulated in COVID-19
| Gene | Log2FC | Function |
|------|--------|----------|
| IFIT1 | 4.23 | Interferon response |
| IFIT3 | 3.87 | Antiviral defense |
| ISG15 | 3.45 | Protein modification |
| MX1 | 3.12 | Antiviral protein |
| OAS1 | 2.98 | RNA degradation |

### Downregulated in COVID-19
| Gene | Log2FC | Function |
|------|--------|----------|
| CD3D | -1.87 | T-cell receptor |
| CD8A | -1.65 | Cytotoxic T cells |
| CD4 | -1.43 | Helper T cells |
| MS4A1 | -1.21 | B-cell marker |
| GNLY | -0.98 | NK cell function |

### Pathway Enrichment
**Top 5 Enriched Pathways:**
1. Interferon Signaling (FDR: 2.1e-09)
2. TNF Signaling (FDR: 4.2e-07)
3. IL-6 Signaling (FDR: 6.8e-06)
4. Antigen Processing (FDR: 1.2e-05)
5. Oxidative Stress (FDR: 1.8e-04)

## Biological Interpretation

The analysis reveals a robust immune response signature in COVID-19 patients:

1. **Innate Immune Activation**: Strong upregulation of interferon-stimulated genes (IFIT1, IFIT3, ISG15, MX1, OAS1) indicates active antiviral defense mechanisms.

2. **T-Cell Dysfunction**: Downregulation of T-cell markers (CD3D, CD4, CD8A) suggests T-cell exhaustion, which is consistent with COVID-19 immunopathology.

3. **Enhanced Inflammation**: Enrichment of TNF and IL-6 signaling pathways indicates heightened inflammatory response.

4. **Monocyte Activation**: Increased expression of inflammatory cytokines and enhanced phagocytic activity.

## Repository Structure

```
single-cell-COVID-Insights/
├── data/
│   ├── raw/                    # Raw 10x data
│   └── processed/              # Processed files
├── plots/                      # All visualizations
├── results/                    # Analysis outputs
├── scripts/
│   ├── enhanced_analysis.py    # Main pipeline
│   └── config.yaml            # Parameters
└── README.md
```

## Getting Started

### Setup Environment
```bash
conda create -n covid_scrnaseq python=3.10
conda activate covid_scrnaseq
pip install -r requirements.txt
```

### Run Analysis
```bash
python enhanced_analysis.py
```

### Generate Plots
```bash
python create_individual_plots.py
```

## What This Shows

This analysis demonstrates:
- **Advanced scRNA-seq processing** with proper QC and statistical validation
- **Comprehensive differential expression analysis** with effect sizes and power analysis
- **Multi-database pathway enrichment** for biological interpretation
- **Publication-ready visualizations** with proper annotations
- **Reproducible workflow** with version-controlled parameters

## Limitations

- Small sample size (n=6) limits generalizability
- No clinical metadata available for severity correlation
- Cross-sectional design (no longitudinal data)
- Limited cell type annotation without marker gene validation

## Future Work

- Validate findings in larger cohorts
- Integrate with clinical metadata
- Perform cell type-specific analysis
- Compare with other COVID-19 scRNA-seq studies

## Contact

**Author**: [Your Name]  
**Email**: [your.email@institution.edu]  
**GitHub**: [https://github.com/mulbagalamaq](https://github.com/mulbagalamaq)

---

*This analysis was performed as part of my bioinformatics portfolio to demonstrate scRNA-seq analysis skills. All code is reproducible and well-documented.*

**Last updated**: July 2024
