#!/usr/bin/env python3
"""
Generate comprehensive plots and metrics for COVID-19 scRNA-seq analysis
Author: Bioinformatics Analysis Pipeline
Date: June 2024
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set scanpy settings
sc.settings.set_figure_params(dpi=100, frameon=False)
sc.settings.verbosity = 3

# Set style
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

def create_quality_control_plots(adata, output_dir="plots"):
    """Generate quality control plots"""
    Path(output_dir).mkdir(exist_ok=True)
    
    # Create QC metrics plot
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Genes per cell
    sc.pl.violin(adata, 'n_genes_by_counts', groupby='sample', ax=axes[0,0])
    axes[0,0].set_title('Genes per Cell by Sample')
    
    # UMI counts per cell
    sc.pl.violin(adata, 'total_counts', groupby='sample', ax=axes[0,1])
    axes[0,1].set_title('UMI Counts per Cell by Sample')
    
    # Mitochondrial percentage
    sc.pl.violin(adata, 'pct_counts_mt', groupby='sample', ax=axes[1,0])
    axes[1,0].set_title('Mitochondrial Content by Sample')
    
    # Scatter plot of genes vs counts
    sc.pl.scatter(adata, 'total_counts', 'n_genes_by_counts', ax=axes[1,1])
    axes[1,1].set_title('Genes vs UMI Counts')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/qc_metrics.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return "Quality control plots generated"

def create_umap_plots(adata, output_dir="plots"):
    """Generate UMAP visualization plots"""
    Path(output_dir).mkdir(exist_ok=True)
    
    # UMAP by sample
    sc.pl.umap(adata, color='sample', size=20, legend_loc='on data')
    plt.title('UMAP by Sample')
    plt.savefig(f'{output_dir}/umap_by_sample.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # UMAP by condition
    sc.pl.umap(adata, color='condition', size=20, legend_loc='on data')
    plt.title('UMAP by Condition (COVID-19 vs Control)')
    plt.savefig(f'{output_dir}/umap_by_condition.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # UMAP by cell type (if available)
    if 'leiden' in adata.obs.columns:
        sc.pl.umap(adata, color='leiden', size=20, legend_loc='on data')
        plt.title('UMAP by Cell Type Clusters')
        plt.savefig(f'{output_dir}/umap_by_cluster.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    return "UMAP plots generated"

def create_heatmap_plots(adata, output_dir="plots"):
    """Generate heatmap visualizations"""
    Path(output_dir).mkdir(exist_ok=True)
    
    # Top marker genes heatmap
    if 'rank_genes_groups' in adata.uns:
        sc.pl.rank_genes_groups_heatmap(adata, n_genes=20, groupby='condition')
        plt.title('Top Marker Genes by Condition')
        plt.savefig(f'{output_dir}/marker_genes_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    return "Heatmap plots generated"

def create_violin_plots(adata, output_dir="plots"):
    """Generate violin plots for key genes"""
    Path(output_dir).mkdir(exist_ok=True)
    
    # Key COVID-19 related genes
    covid_genes = ['IFIT1', 'IFIT3', 'ISG15', 'MX1', 'OAS1', 'CD3D', 'CD4', 'CD8A']
    available_genes = [gene for gene in covid_genes if gene in adata.var_names]
    
    if available_genes:
        sc.pl.violin(adata, available_genes, groupby='condition', rotation=45)
        plt.title('Expression of Key Genes by Condition')
        plt.savefig(f'{output_dir}/key_genes_violin.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    return "Violin plots generated"

def generate_metrics_report(adata, output_dir="results"):
    """Generate comprehensive metrics report"""
    Path(output_dir).mkdir(exist_ok=True)
    
    # Basic statistics
    stats = {
        'total_cells': len(adata),
        'total_genes': len(adata.var),
        'mean_genes_per_cell': adata.obs['n_genes_by_counts'].mean(),
        'median_genes_per_cell': adata.obs['n_genes_by_counts'].median(),
        'mean_umi_per_cell': adata.obs['total_counts'].mean(),
        'median_umi_per_cell': adata.obs['total_counts'].median(),
        'mean_mito_pct': adata.obs['pct_counts_mt'].mean(),
        'median_mito_pct': adata.obs['pct_counts_mt'].median()
    }
    
    # Sample statistics
    sample_stats = adata.obs.groupby('sample').agg({
        'n_genes_by_counts': ['mean', 'median', 'std'],
        'total_counts': ['mean', 'median', 'std'],
        'pct_counts_mt': ['mean', 'median', 'std']
    }).round(2)
    
    # Condition statistics
    condition_stats = adata.obs.groupby('condition').agg({
        'n_genes_by_counts': ['mean', 'median', 'std'],
        'total_counts': ['mean', 'median', 'std'],
        'pct_counts_mt': ['mean', 'median', 'std']
    }).round(2)
    
    # Save statistics
    with open(f'{output_dir}/analysis_metrics.txt', 'w') as f:
        f.write("COVID-19 scRNA-seq Analysis Metrics Report\n")
        f.write("=" * 50 + "\n\n")
        
        f.write("Overall Statistics:\n")
        f.write("-" * 20 + "\n")
        for key, value in stats.items():
            f.write(f"{key}: {value}\n")
        
        f.write("\nSample Statistics:\n")
        f.write("-" * 20 + "\n")
        f.write(sample_stats.to_string())
        
        f.write("\n\nCondition Statistics:\n")
        f.write("-" * 20 + "\n")
        f.write(condition_stats.to_string())
    
    return stats, sample_stats, condition_stats

def create_summary_plots(adata, output_dir="plots"):
    """Create summary visualization plots"""
    Path(output_dir).mkdir(exist_ok=True)
    
    # Cell distribution by sample
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    
    # Sample distribution
    sample_counts = adata.obs['sample'].value_counts()
    axes[0].pie(sample_counts.values, labels=sample_counts.index, autopct='%1.1f%%')
    axes[0].set_title('Cell Distribution by Sample')
    
    # Condition distribution
    condition_counts = adata.obs['condition'].value_counts()
    axes[1].pie(condition_counts.values, labels=condition_counts.index, autopct='%1.1f%%')
    axes[1].set_title('Cell Distribution by Condition')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/cell_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return "Summary plots generated"

def main():
    """Main function to run the analysis"""
    print("Starting COVID-19 scRNA-seq Analysis Visualization...")
    
    # Check if processed data exists
    try:
        # Try to load processed data
        adata = sc.read_h5ad('data/processed/processed_data.h5ad')
        print("Loaded processed data successfully")
    except:
        print("Processed data not found. Please run the analysis notebook first.")
        return
    
    # Generate all plots and metrics
    print("Generating quality control plots...")
    create_quality_control_plots(adata)
    
    print("Generating UMAP plots...")
    create_umap_plots(adata)
    
    print("Generating heatmap plots...")
    create_heatmap_plots(adata)
    
    print("Generating violin plots...")
    create_violin_plots(adata)
    
    print("Generating summary plots...")
    create_summary_plots(adata)
    
    print("Generating metrics report...")
    stats, sample_stats, condition_stats = generate_metrics_report(adata)
    
    print("\nAnalysis complete! Generated files:")
    print("- plots/qc_metrics.png")
    print("- plots/umap_by_sample.png")
    print("- plots/umap_by_condition.png")
    print("- plots/cell_distribution.png")
    print("- results/analysis_metrics.txt")
    
    # Print summary statistics
    print(f"\nSummary Statistics:")
    print(f"Total cells: {stats['total_cells']:,}")
    print(f"Total genes: {stats['total_genes']:,}")
    print(f"Mean genes per cell: {stats['mean_genes_per_cell']:.1f}")
    print(f"Mean UMI per cell: {stats['mean_umi_per_cell']:.1f}")
    print(f"Mean mitochondrial content: {stats['mean_mito_pct']:.1f}%")

if __name__ == "__main__":
    main() 