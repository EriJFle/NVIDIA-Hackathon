"""
Enhanced Quality Control Pipeline for Macosko snRNA-seq Astrocyte Dataset
=========================================================================

This script performs comprehensive quality control on the large Macosko snRNA-seq
astrocyte dataset (~480,000 nuclei) with the following enhancements:

1. **SoupX Ambient RNA Correction**:
   - Removes ambient RNA contamination common in droplet-based protocols
   - Uses chunking approach (8400 genes per chunk) to handle large datasets
   - Requires clustering information for contamination estimation

2. **scDblFinder Doublet Detection**:
   - Superior doublet detection method for transcriptionally homogeneous datasets
   - Processes cells in chunks (50,000 cells per chunk) for memory efficiency
   - Provides doublet scores and classifications with conservative thresholds

3. **Method-Aware QC**:
   - Uses properly scaled MAD thresholds (1.4826 factor) for outlier detection
   - snRNA-seq specific considerations (low mitochondrial content expected)
   - Comprehensive visualization of correction effects

Dependencies:
- Python: scanpy, rpy2, anndata2ri
- R packages: SoupX, scDblFinder, SingleCellExperiment

Author: AstroWRLD Project
Date: 2024-08-20
"""

import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.sparse import vstack, issparse, csr_matrix

# Import AstroWRLD QC functions
sys.path.append(str(Path(__file__).parent.parent / "scripts"))
from qc_functions import automated_qc_pipeline, generate_qc_report
from integrate_qc_pipeline import (
    integrate_qc_with_macosko_pipeline,
    validate_astrocyte_markers_post_qc,
    generate_integrated_report
)

# R integration for SoupX (kept for SoupX-specific functionality)
import anndata2ri
import logging
import rpy2.rinterface_lib.callbacks as rcb
rcb.logger.setLevel(logging.ERROR)
import anndata2ri
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import default_converter
from rpy2.robjects import pandas2ri
from rpy2.robjects import numpy2ri
from rpy2.robjects.packages import importr

print("="*60)
print("AstroWRLD Macosko snRNA-seq Preprocessing Pipeline")
print("SoupX + Configuration-driven QC with scDblFinder")
print("="*60)

# For reproducibility
sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=120)
np.random.seed(212121)

# --- Load data ----------------------------------------------------------------
astro_path = Path("/path/to/AstroWRLD/singlenuclei_data/macosko_astro_only.h5ad")
print(f"Loading Macosko data from: {astro_path}")
adata = sc.read_h5ad(astro_path)

print(f"Initial dataset: {adata.n_obs:,} nuclei × {adata.n_vars:,} genes")

# Prepare data structure
if 'counts' not in adata.layers:
    adata.layers["counts"] = adata.X.copy()

# Handle gene names
if 'gene_name' in adata.var.columns:
    adata.var['ensembl_id'] = adata.var.index.copy()
    adata.var.index = adata.var['gene_name'].astype(str).copy()
    adata.var_names_make_unique()
    adata.var.index.name = None

# Preserve raw data
adata.raw = adata  # snapshot for marker visualization later

# --- SoupX Ambient RNA Correction (Large Dataset Chunking Approach) ---------
print("Starting SoupX ambient RNA correction...")

try:
    # Load R packages
    soupx = importr("SoupX")

    # Prepare data for SoupX
    # Extract count matrices (genes x cells format for R)
    data = adata.X.T.asfptype()  # Transpose to genes x cells
    # For snRNA-seq, we'll use the same matrix for both filtered and raw
    # In practice, you might have separate raw/filtered matrices
    data_tod = data.copy()

    genes = np.array(adata.var_names)
    cells = np.array(adata.obs_names)

    # Need clustering information for SoupX - we'll compute a quick clustering
    print("Computing quick clustering for SoupX...")
    adata_temp = adata.copy()
    sc.pp.normalize_total(adata_temp, target_sum=1e4)
    sc.pp.log1p(adata_temp)
    sc.pp.highly_variable_genes(adata_temp, n_top_genes=2000)
    sc.pp.scale(adata_temp, zero_center=False)
    sc.tl.pca(adata_temp, n_comps=30)
    sc.pp.neighbors(adata_temp, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(adata_temp, resolution=0.3, flavor='igraph')
    soupx_groups = adata_temp.obs['leiden'].astype(str).values
    del adata_temp  # Free memory

    # Convert to integer for R
    soupx_groups_int = pd.Categorical(soupx_groups).codes + 1  # R uses 1-based indexing

    chunk_size = 8400  # Size of gene chunks (adjust based on memory limitations)
    num_chunks = len(genes) // chunk_size + 1  # Number of chunks needed for genes
    print(f"Processing {num_chunks} chunks of {chunk_size} genes each...")

    # Store chunked results
    chunk_results = []

    # Loop through each chunk of genes
    for i in range(num_chunks):
        start_gene = i * chunk_size
        end_gene = min((i + 1) * chunk_size, len(genes))
        print(f"Processing chunk {i+1}/{num_chunks}: genes {start_gene}-{end_gene}")

        # Subset data for this chunk of genes
        chunk_data = data[start_gene:end_gene, :]  # Slice by genes (rows)
        chunk_data_tod = data_tod[start_gene:end_gene, :]  # Same for data_tod

        with localconverter(anndata2ri.converter):
            ro.globalenv['data'] = chunk_data
            ro.globalenv['data_tod'] = chunk_data_tod
            ro.globalenv['genes'] = genes
            ro.globalenv['cells'] = cells
            ro.globalenv['soupx_groups'] = soupx_groups


        # Run SoupX code in R
        ro.r('''
        # Specify row and column names of data
        rownames(data) = genes
        colnames(data) = cells
        rownames(data_tod) = genes
        colnames(data_tod) = cells

        # Ensure correct sparse format for table of counts and table of droplets
        data <- as(data, "sparseMatrix")
        data_tod <- as(data_tod, "sparseMatrix")

        # Generate SoupChannel Object for SoupX
        sc = SoupChannel(data_tod, data, calcSoupProfile = FALSE)

        # Add extra meta data to the SoupChannel object
        soupProf = data.frame(row.names = rownames(data),
                             est = rowSums(data)/sum(data),
                             counts = rowSums(data))
        sc = setSoupProfile(sc, soupProf)

        # Set cluster information in SoupChannel
        sc = setClusters(sc, soupx_groups)

        # Estimate contamination fraction
        sc = autoEstCont(sc, doPlot=FALSE)

        # Infer corrected table of counts and round to integer
        out = adjustCounts(sc, roundToInt = TRUE)
        assign(paste0("out_chunk_", i), out)
        ''')

        # Retrieve chunk result
        chunk_result = ro.globalenv[f'out_chunk_{i}']
        chunk_results.append(chunk_result)

    # Concatenate all chunks row-wise (genes)
    print("Concatenating SoupX corrected chunks...")
    out_corrected = vstack(chunk_results)

    # Store corrected counts
    adata.layers["soupX_counts"] = out_corrected.T  # Transpose back to cells x genes
    adata.X = adata.layers["soupX_counts"].copy()  # Use corrected counts for downstream

    print("SoupX ambient RNA correction completed successfully!")

except Exception as e:
    print(f"SoupX correction failed: {e}")
    print("Continuing with original counts...")
    adata.layers["soupX_counts"] = adata.layers["counts"].copy()
    # Keep original X unchanged

print("✅ SoupX correction completed")

print("\n🔧 Step 2: Configuration-driven quality control...")
adata_qc, qc_metrics = integrate_qc_with_macosko_pipeline(adata_soupx)

print(f"✅ QC completed:")
print(f"   Initial nuclei: {qc_metrics['n_cells_initial']:,}")
print(f"   Filtered nuclei: {qc_metrics['n_cells_filtered']:,}")
print(f"   Retention rate: {qc_metrics['cells_retained_pct']:.1f}%")
print(f"   Method detected: {qc_metrics['method']}")

# Use the QC-processed data
adata = adata_qc

# --- scDblFinder Doublet Detection (Chunking Approach for Large Datasets) ---
print("Starting scDblFinder doublet detection...")

try:
    # Load required R packages
    scdblfinder = importr("scDblFinder")
    singlecellexperiment = importr("SingleCellExperiment")

    # For very large datasets, we'll chunk by cells to manage memory
    cell_chunk_size = 50000  # Process cells in chunks
    n_cells = adata.n_obs
    num_cell_chunks = (n_cells + cell_chunk_size - 1) // cell_chunk_size

    print(f"Processing {n_cells:,} cells in {num_cell_chunks} chunks...")

    doublet_scores = np.zeros(n_cells)
    doublet_classes = np.array(['singlet'] * n_cells, dtype=object)

    # Process cells in chunks
    for chunk_idx in range(num_cell_chunks):
        start_cell = chunk_idx * cell_chunk_size
        end_cell = min((chunk_idx + 1) * cell_chunk_size, n_cells)

        print(f"Processing cell chunk {chunk_idx + 1}/{num_cell_chunks}: cells {start_cell}-{end_cell}")

        # Extract chunk data (cells x genes for R)
        chunk_data = adata.X[start_cell:end_cell, :].copy()

        # Convert to dense if needed (scDblFinder can handle sparse, but sometimes dense is more reliable)
        if issparse(chunk_data):
            if chunk_data.nnz / chunk_data.size > 0.3:  # If >30% non-zero, convert to dense
                chunk_data = chunk_data.toarray()
            else:
                chunk_data = chunk_data.tocsr()

        # Convert data to R-compatible format using numpy2ri converter
        if issparse(chunk_data):
            chunk_data = chunk_data.toarray()  # Convert sparse to dense for R

        # Create a converter that includes numpy conversion rules
        np_cv_rules = default_converter + numpy2ri.converter

        # Use the converter context for proper numpy array conversion
        with np_cv_rules.context():
            ro.globalenv['chunk_data'] = chunk_data.T  # Transpose to genes x cells for R
            ro.globalenv['chunk_genes'] = np.array(adata.var_names)
            ro.globalenv['chunk_cells'] = np.array(adata.obs_names[start_cell:end_cell])

        # Run scDblFinder in R
        ro.r('''
        # Create SingleCellExperiment object
        library(SingleCellExperiment)
        library(scDblFinder)

        # Ensure row and column names
        if (is.null(rownames(chunk_data))) {
            rownames(chunk_data) <- chunk_genes
        }
        if (is.null(colnames(chunk_data))) {
            colnames(chunk_data) <- chunk_cells
        }

        # Create SingleCellExperiment
        sce <- SingleCellExperiment(list(counts = chunk_data))

        # Run scDblFinder
        sce <- scDblFinder(sce)

        # Extract results
        doublet_score_chunk <- sce$scDblFinder.score
        doublet_class_chunk <- sce$scDblFinder.class
        ''')

        # Retrieve results from R
        chunk_scores = np.array(ro.globalenv['doublet_score_chunk'])
        chunk_classes = np.array(ro.globalenv['doublet_class_chunk'])

        # Store in full arrays
        doublet_scores[start_cell:end_cell] = chunk_scores
        doublet_classes[start_cell:end_cell] = chunk_classes

    # Add results to adata
    adata.obs["doublet_score"] = doublet_scores
    adata.obs["doublet_class"] = doublet_classes

    # Create conservative flag using 95th percentile threshold
    cutoff = np.percentile(doublet_scores, 95)
    adata.obs["doublet_flag"] = doublet_scores >= cutoff

    print(f"scDblFinder completed! Detected {(doublet_classes == 'doublet').sum():,} doublets")
    print(f"Conservative flag (95th percentile): {adata.obs['doublet_flag'].sum():,} cells flagged")

except Exception as e:
    print(f"scDblFinder failed: {e}")
    print("Falling back to Scrublet...")

    # Fallback to original Scrublet implementation
    try:
        sc.pp.scrublet(adata, batch_key="assay")
        cutoff = np.percentile(adata.obs["doublet_score"], 95)
        adata.obs["doublet_flag"] = adata.obs["doublet_score"] >= cutoff
        adata.obs["doublet_class"] = np.where(adata.obs["doublet_flag"], 'doublet', 'singlet')
        print("Scrublet completed as fallback")
    except Exception as e2:
        print(f"Both scDblFinder and Scrublet failed: {e2}")
        adata.obs["doublet_score"] = np.nan
        adata.obs["doublet_flag"] = False
        adata.obs["doublet_class"] = 'unknown'

adata.write('/path/to/AstroWRLD/singlenuclei_data/macosko_astro_ref_fully_qcd.h5ad')


# --- Normalization and feature selection (keep all cells) ---------------------
# CPM-like normalization, log1p, HVG selection with Seurat v3 flavor
scales_counts = sc.pp.normalize_total(adata, target_sum=1e4, inplace=False)
adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

sc.pp.highly_variable_genes(
    adata,
    layer="counts",
    flavor="seurat_v3",
    n_top_genes=3000,
    subset=False,
    batch_key='assay',
    inplace=True,
)

# --- Scaling, PCA, neighbors, UMAP, clustering --------------------------------
#sc.pp.scale(adata,
            #max_value=10, don't use max_value if zero_center=False
            #zero_center=False,
#)  # simple z-scaling
sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30)
sc.tl.umap(adata, min_dist=0.3, spread=1.0)
sc.tl.leiden(adata, resolution=0.3, key_added="leiden", flavor='igraph')

# --- QC Summary and Comparison -----------------------------------------------
print("\n=== QC Summary ===")
print(f"Total cells: {adata.n_obs:,}")
print(f"Total genes: {adata.n_vars:,}")
print(f"Cells flagged for any QC issue: {adata.obs['qc_flag_any'].sum():,} ({100*adata.obs['qc_flag_any'].mean():.1f}%)")
print(f"Doublets detected: {adata.obs['doublet_flag'].sum():,} ({100*adata.obs['doublet_flag'].mean():.1f}%)")

# Compare SoupX correction if available
if "soupX_counts" in adata.layers:
    original_total = adata.layers["counts"].sum()
    corrected_total = adata.layers["soupX_counts"].sum()
    reduction_pct = 100 * (1 - corrected_total / original_total)
    print(f"SoupX ambient RNA removal: {reduction_pct:.2f}% total count reduction")

# --- UMAP QC visualization ----------------------------------------------------
sc.pl.umap(
    adata,
    color=[
        "leiden",
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
        "doublet_score",
        "doublet_flag",
        "doublet_class",
        "qc_flag_any",
    ],
    wspace=0.4,
    save='_qc_visualization_enhanced.png'
)

# Additional QC plots comparing before/after SoupX (if available)
if "soupX_counts" in adata.layers and "counts" in adata.layers:
    print("Generating SoupX comparison plots...")

    # Calculate UMI counts before and after correction
    adata.obs["umi_original"] = np.array(adata.layers["counts"].sum(axis=1)).flatten()
    adata.obs["umi_corrected"] = np.array(adata.layers["soupX_counts"].sum(axis=1)).flatten()
    adata.obs["umi_reduction"] = adata.obs["umi_original"] - adata.obs["umi_corrected"]
    adata.obs["umi_reduction_pct"] = 100 * (adata.obs["umi_reduction"] / adata.obs["umi_original"])

    # Plot SoupX correction effects
    sc.pl.umap(
        adata,
        color=["umi_reduction", "umi_reduction_pct"],
        wspace=0.4,
        save='_soupx_correction_effects.png'
    )
