import matplotlib.pyplot as plt
import scanpy as sc

adata = sc.read_h5ad("gl261.h5ad")

rsc.get.anndata_to_GPU(adata)

rsc.pp.normalize_total(adata, target_sum=1e4)
rsc.pp.log1p(adata)
rsc.pp.highly_variable_genes(adata, n_top_genes=2000)

# Run NMF on normalized, log-transformed data
nmf_gene_programs_gpu(
    adata,
    layer="X",
    n_programs=50
)

rsc.pp.filter_genes(adata, min_cells=3)
rsc.pp.pca(adata)
rsc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
rsc.tl.umap(adata, min_dist=0.3)

for i in range(50):
    key = f"Program_{i}"
    adata.obs[key] = adata.obsm["NMF_cell_programs"][:, i]
    sc.pl.umap(adata, color=f'Program_{i}', cmap='viridis', save=f"_{key}_nmf.png", vmin=0, vmax='p99')