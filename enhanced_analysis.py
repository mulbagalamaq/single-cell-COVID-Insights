#!/usr/bin/env python3
"""
Enhanced COVID-19 scRNA-seq Analysis with Statistical Validation
Author: Bioinformatics Analysis Pipeline
Date: June 2024
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import yaml
import warnings
import time
import psutil
import os
from datetime import datetime
from scipy import stats
from statsmodels.stats.power import TTestPower
import gseapy as gp

warnings.filterwarnings('ignore')

class EnhancedCOVIDAnalysis:
    """Enhanced COVID-19 scRNA-seq analysis with comprehensive statistical validation"""
    
    def __init__(self, config_file="config.yaml"):
        """Initialize analysis with configuration"""
        self.config = self.load_config(config_file)
        self.setup_reproducibility()
        self.create_directories()
        self.start_time = time.time()
        
    def load_config(self, config_file):
        """Load configuration from YAML file"""
        with open(config_file, 'r') as f:
            return yaml.safe_load(f)
    
    def setup_reproducibility(self):
        """Setup for reproducible analysis"""
        # Set random seeds
        np.random.seed(self.config['analysis_parameters']['reproducibility']['random_seed'])
        
        # Set scanpy settings
        sc.settings.set_figure_params(
            dpi=self.config['analysis_parameters']['visualization']['dpi'],
            frameon=False
        )
        sc.settings.verbosity = 3
        
        # Set matplotlib style
        plt.style.use(self.config['analysis_parameters']['visualization']['style'])
        sns.set_palette(self.config['analysis_parameters']['visualization']['color_palette'])
        
        # Log system information
        self.log_system_info()
    
    def create_directories(self):
        """Create necessary directories"""
        paths = self.config['data_paths']
        for path in paths.values():
            Path(path).mkdir(parents=True, exist_ok=True)
    
    def log_system_info(self):
        """Log system and analysis information"""
        log_file = f"{self.config['data_paths']['logs']}/analysis_log.txt"
        
        with open(log_file, 'w') as f:
            f.write(f"COVID-19 scRNA-seq Analysis Log\n")
            f.write(f"=" * 50 + "\n")
            f.write(f"Analysis Date: {datetime.now()}\n")
            f.write(f"Python Version: {os.sys.version}\n")
            f.write(f"Scanpy Version: {sc.__version__}\n")
            f.write(f"Random Seed: {self.config['analysis_parameters']['reproducibility']['random_seed']}\n")
            f.write(f"Memory Usage: {psutil.Process().memory_info().rss / 1024**3:.2f} GB\n")
            f.write(f"CPU Cores: {psutil.cpu_count()}\n")
    
    def load_data(self, data_files):
        """Load and combine data files"""
        print("Loading data...")
        
        adata_list = []
        for file_path in data_files:
            try:
                adata = sc.read_10x_h5(file_path)
                sample_name = Path(file_path).stem
                adata.obs['sample'] = sample_name
                adata.var_names_make_unique()
                adata_list.append(adata)
                print(f"Loaded {file_path}: {adata.n_obs} cells, {adata.n_vars} genes")
            except Exception as e:
                print(f"Error loading {file_path}: {e}")
        
        # Combine all datasets
        if adata_list:
            adata = adata_list[0].concatenate(adata_list[1:], join='outer')
            print(f"Combined dataset: {adata.n_obs} cells, {adata.n_vars} genes")
            return adata
        else:
            raise ValueError("No data files could be loaded")
    
    def enhanced_qc_metrics(self, adata):
        """Enhanced quality control with statistical validation"""
        print("Performing enhanced quality control...")
        
        # Calculate additional QC metrics
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
        adata.var['hb'] = adata.var_names.str.startswith(('HB', 'HBA', 'HBB'))
        
        # Enhanced QC calculation
        sc.pp.calculate_qc_metrics(
            adata, 
            qc_vars=['mt', 'ribo', 'hb'], 
            percent_top=[1, 5, 10, 50], 
            log1p=False, 
            inplace=True
        )
        
        # Statistical outlier detection using IQR method
        qc_params = self.config['analysis_parameters']['quality_control']
        
        # Genes per cell outliers
        q1 = adata.obs['n_genes_by_counts'].quantile(0.25)
        q3 = adata.obs['n_genes_by_counts'].quantile(0.75)
        iqr = q3 - q1
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr
        
        gene_outliers = (adata.obs['n_genes_by_counts'] < lower_bound) | \
                       (adata.obs['n_genes_by_counts'] > upper_bound)
        
        # UMI count outliers
        q1_umi = adata.obs['total_counts'].quantile(0.25)
        q3_umi = adata.obs['total_counts'].quantile(0.75)
        iqr_umi = q3_umi - q1_umi
        lower_bound_umi = q1_umi - 1.5 * iqr_umi
        upper_bound_umi = q3_umi + 1.5 * iqr_umi
        
        umi_outliers = (adata.obs['total_counts'] < lower_bound_umi) | \
                      (adata.obs['total_counts'] > upper_bound_umi)
        
        # Combine outliers
        total_outliers = gene_outliers | umi_outliers
        
        print(f"Statistical outliers detected: {total_outliers.sum()} cells")
        print(f"Gene count outliers: {gene_outliers.sum()} cells")
        print(f"UMI count outliers: {umi_outliers.sum()} cells")
        
        # Store QC statistics
        self.qc_stats = {
            'total_cells_before': len(adata),
            'outliers_detected': total_outliers.sum(),
            'outlier_percentage': (total_outliers.sum() / len(adata)) * 100,
            'gene_outliers': gene_outliers.sum(),
            'umi_outliers': umi_outliers.sum()
        }
        
        return adata, total_outliers
    
    def assess_batch_effects(self, adata):
        """Assess and visualize batch effects"""
        print("Assessing batch effects...")
        
        # ANOVA test for batch effects
        samples = adata.obs['sample'].unique()
        gene_counts = [adata[adata.obs['sample'] == sample].obs['n_genes_by_counts'] 
                      for sample in samples]
        
        f_stat, p_value = stats.f_oneway(*gene_counts)
        
        print(f"Batch effect ANOVA - F-statistic: {f_stat:.3f}, p-value: {p_value:.3e}")
        
        # Store batch effect statistics
        self.batch_effects = {
            'f_statistic': f_stat,
            'p_value': p_value,
            'significant': p_value < 0.05
        }
        
        return f_stat, p_value
    
    def filter_data(self, adata, outliers_mask):
        """Filter data based on QC metrics"""
        print("Filtering data...")
        
        qc_params = self.config['analysis_parameters']['quality_control']
        
        # Remove outliers
        adata = adata[~outliers_mask]
        
        # Apply standard filters
        sc.pp.filter_cells(adata, min_genes=qc_params['min_genes'])
        sc.pp.filter_cells(adata, min_counts=qc_params['min_counts'])
        adata = adata[adata.obs.pct_counts_mt < qc_params['max_mito_pct']]
        sc.pp.filter_genes(adata, min_cells=qc_params['min_cells'])
        
        print(f"After filtering: {adata.n_obs} cells, {adata.n_vars} genes")
        
        return adata
    
    def preprocess_data(self, adata):
        """Preprocess data with normalization and scaling"""
        print("Preprocessing data...")
        
        preprocess_params = self.config['analysis_parameters']['preprocessing']
        
        # Normalize
        sc.pp.normalize_total(adata, target_sum=preprocess_params['normalization_target'])
        sc.pp.log1p(adata)
        
        # Identify highly variable genes
        if preprocess_params['highly_variable_genes']:
            sc.pp.highly_variable_genes(adata, n_top_genes=preprocess_params['n_top_genes'])
        
        # Scale data
        if preprocess_params['scale_data']:
            sc.pp.scale(adata, max_value=10)
        
        return adata
    
    def dimensionality_reduction(self, adata):
        """Perform dimensionality reduction"""
        print("Performing dimensionality reduction...")
        
        dim_params = self.config['analysis_parameters']['dimensionality_reduction']
        
        # PCA
        sc.tl.pca(adata, n_comps=dim_params['pca_n_comps'], svd_solver='arpack')
        
        # Compute neighborhood graph
        sc.pp.neighbors(
            adata, 
            n_neighbors=dim_params['neighbors_n_neighbors'],
            n_pcs=dim_params['neighbors_n_pcs']
        )
        
        # UMAP
        sc.tl.umap(adata)
        
        return adata
    
    def cluster_cells(self, adata):
        """Cluster cells using Leiden algorithm"""
        print("Clustering cells...")
        
        cluster_params = self.config['analysis_parameters']['clustering']
        
        sc.tl.leiden(
            adata, 
            resolution=cluster_params['resolution'],
            random_state=cluster_params['random_state']
        )
        
        return adata
    
    def enhanced_de_analysis(self, adata, groupby='condition', method='wilcoxon'):
        """Enhanced differential expression with comprehensive statistics"""
        print("Performing enhanced differential expression analysis...")
        
        de_params = self.config['analysis_parameters']['differential_expression']
        
        # Perform DE analysis
        sc.tl.rank_genes_groups(adata, groupby, method=method)
        
        # Get results
        de_results = sc.get.rank_genes_groups_df(adata, group='COVID')
        
        # Calculate effect sizes (Cohen's d)
        def cohens_d(group1, group2):
            n1, n2 = len(group1), len(group2)
            if n1 < 2 or n2 < 2:
                return 0
            pooled_std = np.sqrt(((n1-1)*group1.var() + (n2-1)*group2.var()) / (n1+n2-2))
            if pooled_std == 0:
                return 0
            return (group1.mean() - group2.mean()) / pooled_std
        
        # Calculate effect sizes for top genes
        top_genes = de_results.head(100)['names'].tolist()
        effect_sizes = {}
        
        for gene in top_genes:
            if gene in adata.var_names:
                covid_expr = adata[adata.obs['condition'] == 'COVID', gene].X.toarray().flatten()
                ctrl_expr = adata[adata.obs['condition'] == 'Control', gene].X.toarray().flatten()
                effect_sizes[gene] = cohens_d(covid_expr, ctrl_expr)
        
        de_results['effect_size'] = de_results['names'].map(effect_sizes)
        
        # Add significance categories
        de_results['significance'] = 'Not Significant'
        de_results.loc[de_results['pvals_adj'] < 0.05, 'significance'] = 'Significant'
        de_results.loc[de_results['pvals_adj'] < 0.01, 'significance'] = 'Highly Significant'
        
        # Store DE statistics
        self.de_stats = {
            'total_genes': len(de_results),
            'significant_genes': (de_results['pvals_adj'] < 0.05).sum(),
            'highly_significant': (de_results['pvals_adj'] < 0.01).sum(),
            'upregulated': (de_results['logfoldchanges'] > 0).sum(),
            'downregulated': (de_results['logfoldchanges'] < 0).sum(),
            'mean_effect_size': de_results['effect_size'].mean()
        }
        
        return de_results
    
    def power_analysis(self, adata, alpha=0.05, power=0.8):
        """Calculate statistical power for differential expression"""
        print("Performing power analysis...")
        
        # Calculate effect sizes for all genes
        effect_sizes = []
        for gene in adata.var_names[:100]:  # Sample first 100 genes for speed
            covid_expr = adata[adata.obs['condition'] == 'COVID', gene].X.toarray().flatten()
            ctrl_expr = adata[adata.obs['condition'] == 'Control', gene].X.toarray().flatten()
            
            if len(covid_expr) > 0 and len(ctrl_expr) > 0:
                pooled_std = np.sqrt(((len(covid_expr)-1)*covid_expr.var() + 
                                    (len(ctrl_expr)-1)*ctrl_expr.var()) / 
                                   (len(covid_expr)+len(ctrl_expr)-2))
                if pooled_std > 0:
                    effect_size = abs(covid_expr.mean() - ctrl_expr.mean()) / pooled_std
                    effect_sizes.append(effect_size)
        
        if effect_sizes:
            # Calculate power
            power_analysis = TTestPower()
            median_effect = np.median(effect_sizes)
            required_n = power_analysis.solve_power(
                effect_size=median_effect,
                alpha=alpha,
                power=power
            )
            
            print(f"Median effect size: {median_effect:.3f}")
            print(f"Required sample size per group for {power*100}% power: {required_n:.0f}")
            
            self.power_stats = {
                'median_effect_size': median_effect,
                'required_sample_size': required_n,
                'current_power': power_analysis.power(
                    effect_size=median_effect,
                    nobs=len(adata) // 2,
                    alpha=alpha
                )
            }
            
            return median_effect, required_n
        else:
            print("Could not calculate power analysis")
            return None, None
    
    def comprehensive_pathway_analysis(self, de_results, organism='Human'):
        """Comprehensive pathway enrichment analysis"""
        print("Performing comprehensive pathway analysis...")
        
        pathway_params = self.config['analysis_parameters']['pathway_analysis']
        
        # Get significant genes
        sig_genes = de_results[de_results['pvals_adj'] < pathway_params['fdr_threshold']]['names'].tolist()
        
        if len(sig_genes) < 10:
            print("Too few significant genes for pathway analysis")
            return {}
        
        # Run enrichment analysis
        enrichment_results = {}
        
        for db in pathway_params['databases']:
            try:
                enr = gp.enrichr(
                    gene_list=sig_genes,
                    gene_sets=[db],
                    organism=organism,
                    outdir=None
                )
                enrichment_results[db] = enr.results
                print(f"Completed {db}: {len(enr.results)} enriched pathways")
            except Exception as e:
                print(f"Failed to analyze {db}: {e}")
        
        # Store pathway statistics
        total_pathways = sum(len(results) for results in enrichment_results.values())
        significant_pathways = sum(
            len(results[results['Adjusted P-value'] < pathway_params['fdr_threshold']]) 
            for results in enrichment_results.values()
        )
        
        self.pathway_stats = {
            'total_pathways_tested': total_pathways,
            'significant_pathways': significant_pathways,
            'databases_analyzed': len(enrichment_results)
        }
        
        return enrichment_results
    
    def generate_publication_plots(self, adata, de_results, output_dir="plots"):
        """Generate publication-quality visualizations"""
        print("Generating publication-quality plots...")
        
        Path(output_dir).mkdir(exist_ok=True)
        
        # Create comprehensive figure
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        
        # 1. UMAP by condition
        sc.pl.umap(adata, color='condition', ax=axes[0,0], show=False)
        axes[0,0].set_title('Cell Distribution by Condition', fontsize=14, fontweight='bold')
        
        # 2. Quality metrics
        sc.pl.violin(adata, 'n_genes_by_counts', groupby='condition', ax=axes[0,1], show=False)
        axes[0,1].set_title('Genes per Cell', fontsize=14, fontweight='bold')
        
        # 3. Top marker genes
        sc.pl.rank_genes_groups(adata, n_genes=10, ax=axes[0,2], show=False)
        axes[0,2].set_title('Top Marker Genes', fontsize=14, fontweight='bold')
        
        # 4. Cell type proportions
        if 'leiden' in adata.obs.columns:
            cell_counts = adata.obs['leiden'].value_counts()
            axes[1,0].pie(cell_counts.values, labels=cell_counts.index, autopct='%1.1f%%')
            axes[1,0].set_title('Cell Type Distribution', fontsize=14, fontweight='bold')
        
        # 5. Effect size distribution
        if 'effect_size' in de_results.columns:
            axes[1,1].hist(de_results['effect_size'].dropna(), bins=30, alpha=0.7)
            axes[1,1].set_title('Effect Size Distribution', fontsize=14, fontweight='bold')
            axes[1,1].set_xlabel('Cohen\'s d')
            axes[1,1].set_ylabel('Frequency')
        
        # 6. Statistical summary
        summary_text = f"""
        Total Cells: {len(adata):,}
        Total Genes: {len(adata.var):,}
        Significant DEGs: {self.de_stats['significant_genes']:,}
        Median Effect Size: {self.de_stats['mean_effect_size']:.3f}
        """
        axes[1,2].text(0.1, 0.5, summary_text, transform=axes[1,2].transAxes, 
                       fontsize=12, verticalalignment='center')
        axes[1,2].set_title('Analysis Summary', fontsize=14, fontweight='bold')
        axes[1,2].axis('off')
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/comprehensive_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        return "Publication-quality plots generated"
    
    def save_results(self, adata, de_results, pathway_results):
        """Save all results and statistics"""
        print("Saving results...")
        
        # Save processed data
        adata.write(f"{self.config['data_paths']['processed_data']}/processed_data.h5ad")
        
        # Save DE results
        de_results.to_csv(f"{self.config['data_paths']['results']}/differential_expression.csv")
        
        # Save pathway results
        for db, results in pathway_results.items():
            results.to_csv(f"{self.config['data_paths']['results']}/pathway_{db}.csv")
        
        # Save statistics
        def to_serializable(val):
            if isinstance(val, (np.integer, np.floating)):
                return float(val)
            elif isinstance(val, (np.bool_, bool)):
                return bool(val)
            return val
        stats = {
            'qc_stats': {k: to_serializable(v) for k, v in self.qc_stats.items()},
            'batch_effects': {k: to_serializable(v) for k, v in self.batch_effects.items()},
            'de_stats': {k: to_serializable(v) for k, v in self.de_stats.items()},
            'power_stats': {k: to_serializable(v) for k, v in getattr(self, 'power_stats', {}).items()},
            'pathway_stats': {k: to_serializable(v) for k, v in getattr(self, 'pathway_stats', {}).items()},
            'processing_time': float(time.time() - self.start_time)
        }
        
        import json
        with open(f"{self.config['data_paths']['results']}/analysis_statistics.json", 'w') as f:
            json.dump(stats, f, indent=2)
        
        print("Results saved successfully")
    
    def run_complete_analysis(self, data_files):
        """Run the complete enhanced analysis pipeline"""
        print("Starting enhanced COVID-19 scRNA-seq analysis...")
        
        # Load data
        adata = self.load_data(data_files)
        
        # Add condition information
        adata.obs['condition'] = ['COVID' if 'nCoV' in x else 'Control' for x in adata.obs['sample']]
        
        # Enhanced QC
        adata, outliers = self.enhanced_qc_metrics(adata)
        self.assess_batch_effects(adata)
        adata = self.filter_data(adata, outliers)
        
        # Preprocessing
        adata = self.preprocess_data(adata)
        
        # Dimensionality reduction
        adata = self.dimensionality_reduction(adata)
        
        # Clustering
        adata = self.cluster_cells(adata)
        
        # Enhanced DE analysis
        de_results = self.enhanced_de_analysis(adata)
        
        # Power analysis
        self.power_analysis(adata)
        
        # Pathway analysis
        pathway_results = self.comprehensive_pathway_analysis(de_results)
        
        # Generate plots
        self.generate_publication_plots(adata, de_results)
        
        # Save results
        self.save_results(adata, de_results, pathway_results)
        
        print("Enhanced analysis completed successfully!")
        print(f"Total processing time: {(time.time() - self.start_time)/60:.1f} minutes")
        
        return adata, de_results, pathway_results

def main():
    """Main function to run the enhanced analysis"""
    
    # Initialize analysis
    analysis = EnhancedCOVIDAnalysis()
    
    # Define data files (update these paths to match your data)
    data_files = [
        "data/raw/nCoV_PBMC_1.h5",
        "data/raw/nCoV_PBMC_15.h5",
        "data/raw/nCoV_PBMC_17.h5",
        "data/raw/Normal_PBMC_5.h5",
        "data/raw/Normal_PBMC_13.h5",
        "data/raw/Normal_PBMC_14.h5"
    ]
    
    # Run analysis
    try:
        adata, de_results, pathway_results = analysis.run_complete_analysis(data_files)
        print("Analysis completed successfully!")
    except Exception as e:
        print(f"Analysis failed: {e}")
        raise

if __name__ == "__main__":
    main() 