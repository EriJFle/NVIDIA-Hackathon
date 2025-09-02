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
    integrate_qc_with_allen_pipeline,
    validate_astrocyte_markers_post_qc,
    generate_integrated_report
)

# R integration for scDblFinder
import anndata2ri
import logging
import rpy2.rinterface_lib.callbacks as rcb
rcb.logger.setLevel(logging.ERROR)
anndata2ri.activate()
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

# Activate automatic pandas <-> R dataframe conversion
pandas2ri.activate()

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=120)
np.random.seed(212121)

print("="*60)
print("AstroWRLD Allen Brain Atlas scRNA-seq Preprocessing Pipeline")
print("Configuration-driven QC with scDblFinder")
print("="*60)

astro_path = Path("/path/to/abc_cache/allen_wmb10x_astro_ref.h5ad")
print(f"Loading Allen data from: {astro_path}")
adata = sc.read_h5ad(astro_path)

print(f"Initial dataset: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# Prepare data structure
if 'counts' not in adata.layers:
    adata.layers["counts"] = adata.X.copy()

# Handle gene names
if 'gene_symbol' in adata.var.columns:
    adata.var['ensembl_id'] = adata.var.index.copy()
    adata.var.index = adata.var['gene_symbol'].astype(str).copy()
    adata.var_names_make_unique()
    adata.var.index.name = None

# Preserve raw data
adata.raw = adata  # snapshot for marker visualization later

print("\n🔧 Step 1: Configuration-driven quality control...")
adata_qc, qc_metrics = integrate_qc_with_allen_pipeline(adata)

print(f"✅ QC completed:")
print(f"   Initial cells: {qc_metrics['n_cells_initial']:,}")
print(f"   Filtered cells: {qc_metrics['n_cells_filtered']:,}")
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

    # For Allen dataset (~300k cells), use moderate chunking
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

        # Pass data to R
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

adata.write('/path/to/AstroWRLD/allen_wmb_astro_ref_fully_qcd.h5ad')


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

# Generate integrated report\nprint(\"\\n📊 Generating comprehensive QC report...\")\nreport_path = \"results/allen_integrated_qc_report.txt\"\nos.makedirs(\"results\", exist_ok=True)\ngenerate_integrated_report(adata, qc_metrics, marker_validation, report_path)\nprint(f\"📋 Detailed report saved: {report_path}\")\n\n# --- UMAP QC visualization ----------------------------------------------------
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
    save='_allen_qc_visualization_enhanced.png'
)


# --- Differential expression per cluster --------------------------------------
# Wilcoxon rank-sum on raw-backed data structure
# (Scanpy uses the current .X, which is log1p-normalized; adata.raw preserves counts snapshot)
sc.tl.rank_genes_groups(
    adata,
    groupby="leiden",
    method="wilcoxon",
    use_raw=True,
    pts=True,
    corr_method="benjamini-hochberg",
    tie_correct=True,
)

# Inspect top markers per cluster (on screen)
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)
