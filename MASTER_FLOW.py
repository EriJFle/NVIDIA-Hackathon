import os
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import median_abs_deviation

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=150)
sc.settings.figdir = "real_run/figures"

np.random.seed(212121)

# ====================================================
#               scVI Run - Biolean w/ 3230 hvgs
# ====================================================

import scvi
import torch
torch.set_float32_matmul_precision('high')
import rapids_singlecell as rsc

scvi.settings.seed = 212121

adata = sc.read_h5ad("real_run/concat_with_raw.h5ad")

df = pd.read_csv("real_run/shared_hvgs_only.csv")

hvgs = df['gene_symbol']

# Can keep just hvgs in adata, adata.raw has all genes
adata = adata[:, hvgs].copy()

adata.write("real_run/scVI_ready_raw_full_genes.h5ad")

# # 1) Bio-leaning (usually improves NMI/ARI/BRAS without tanking iLISI/kBET)

scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="modality",
    categorical_covariate_keys=["dataset"],
)

model = scvi.model.SCVI(
    adata,
    n_layers=2,
    n_hidden=256,
    n_latent=56,                # slightly lower to sharpen structure
    dropout_rate=0.1,
    gene_likelihood="nb",
    dispersion="gene",          # less batch-specific variance absorption
    encode_covariates=False,    # keep nuisance out of the encoder
    use_layer_norm="both",
    use_batch_norm="none",
    use_observed_lib_size=True,
)

plan_kwargs = dict(
    n_steps_kl_warmup=800,      # milder warmup than 2000
    reduce_lr_on_plateau=True,
    lr_factor=0.7,
    lr_patience=6,
    lr_scheduler_metric="elbo_validation",
    weight_decay=1e-6,
)

model.train(
    max_epochs=80,
    early_stopping=True,
    early_stopping_monitor="elbo_validation",
    early_stopping_patience=18,
    check_val_every_n_epoch=1,
    plan_kwargs=plan_kwargs,
    gradient_clip_val=1.0,
    accelerator="gpu",
    batch_size=2048,
    precision="16-mixed",
)

adata.obsm["X_scVI_biolean"] = model.get_latent_representation()
model.save("real_run/hvgs3230_biolean", overwrite=True)

rsc.pp.neighbors(
    adata,
    use_rep="X_scVI_biolean",
)

rsc.tl.umap(
    adata,
    min_dist=0.3,
)

sc.pl.umap(
    adata,
    color=['region_broad', 'dataset'],
    save="_biolean_umap.png",
)

adata.write("real_run/thru_scvi_umap.h5ad")


# ====================================================
#      scVI Run - Mixing-optimized config (suggested)
# ====================================================

# Reload HVG-filtered data to avoid scVI registry conflicts
adata_mix = sc.read_h5ad("real_run/thru_scvi_umap.h5ad")

# Setup: use dataset as batch; modality as covariate
scvi.model.SCVI.setup_anndata(
    adata_mix,
    layer='counts',
    batch_key='dataset',
    categorical_covariate_keys=['modality'],
)

model_mix = scvi.model.SCVI(
    adata_mix,
    n_layers=2,
    n_hidden=256,
    n_latent=64,
    dropout_rate=0.1,
    gene_likelihood='nb',
    dispersion='gene-batch',
    encode_covariates=True,
    use_layer_norm='both',
    use_batch_norm='encoder',
    use_observed_lib_size=True,
)

plan_kwargs_mix = dict(
    n_steps_kl_warmup=400,
    reduce_lr_on_plateau=True,
    lr_factor=0.7,
    lr_patience=6,
    lr_scheduler_metric='elbo_validation',
    weight_decay=1e-6,
)

model_mix.train(
    max_epochs=120,
    early_stopping=True,
    early_stopping_monitor='elbo_validation',
    early_stopping_patience=18,
    check_val_every_n_epoch=1,
    plan_kwargs=plan_kwargs_mix,
    accelerator='gpu',
    batch_size=2048,
    precision='16-mixed',
)

# Store latent under a distinct key
adata_mix.obsm['X_scVI_mixing'] = model_mix.get_latent_representation()
model_mix.save('real_run/hvgs3230_mix_optimized', overwrite=True)

# Neighbors/UMAP tuned for mixing
rsc.pp.neighbors(adata_mix, use_rep='X_scVI_mixing', metric='cosine', n_neighbors=30)
rsc.tl.umap(
    adata_mix,
    min_dist=0.6,
    key_added="X_umap_mix"
)

sc.pl.umap(
    adata_mix,
    color=['dataset', 'region_broad', 'qc_flag_any', 'doublet_flag'],
    ncols=2,
    save='_mixing_umap.png',
)

adata_mix.write('real_run/thru_scvi_umap_biolean_and_mixing.h5ad')

# ==============================================================================
#                           scib-metrics run
# ==============================================================================
adata = sc.read_h5ad("real_run/thru_scvi_umap_biolean_and_mixing.h5ad")

model = scvi.model.SCVI.load("real_run/hvgs3230_biolean", adata=adata)

adata.layers['scvi_norm_10k'] = model.get_normalized_expression(
    library_size=1e4,
    return_numpy=True,
)

adata.layers['scvi_norm_10k_log1p'] = np.log1p(adata.layers['scvi_norm_10k'])

from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection

import jax
print(jax.devices())
print(jax.default_backend())

rsc.pp.pca(
    adata,
    n_comps=50,
    zero_center=False,
    mask_var="highly_variable",
    layer="scvi_norm_10k_log1p"
)

rsc.pp.harmony_integrate(
    adata,
    key="dataset",
    basis="X_pca",
    adjusted_basis="harmony",
    theta=1.0,
    sigma=0.2,
    max_iter_harmony=20,
    verbose=True,
)

adata.obsm["Unintegrated"] = adata.obsm["X_pca"].copy()

adata.write("real_run/thru_pca_and_harmony_scvi_norm_10k_log1p_as_X.h5ad")

biocons = BioConservation(isolated_labels=False)

bm = Benchmarker(
    adata,
    batch_key="dataset",
    label_key="region_broad",
    embedding_obsm_keys=[ "X_scVI_biolean", "X_scVI_mixing", "harmony", "Unintegrated"],
    pre_integrated_embedding_obsm_key="X_pca",
    bio_conservation_metrics=biocons,
    batch_correction_metrics=BatchCorrection(),
    n_jobs=-1,
)

bm.benchmark()

df = bm.get_results(min_max_scale=False, clean_names=True)

print(df)

df.to_csv("real_run_2_scVI.csv")

bm.plot_results_table(min_max_scale=False, show=True, save_dir="real_run/bm_results")

# ====================================================
#               Leiden Clustering
# ====================================================
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=150)
sc.settings.figdir = "real_run/figures"

np.random.seed(212121)

adata = sc.read_h5ad("real_run/thru_scvi_umap_biolean_and_mixing.h5ad")

rsc.pp.neighbors(
    adata,
    use_rep="X_scVI_biolean",
)

rsc.tl.umap(
    adata,
    min_dist=0.3,
)

rsc.tl.leiden(
    adata,
    resolution=0.3,
    random_state=212121,
    key_added='leiden_03',
)

rsc.tl.leiden(
    adata,
    resolution=0.4,
    random_state=212121,
    key_added='leiden_04',
)

sc.pl.umap(
    adata,
    color=['leiden_03', 'region_broad',
    'leiden_04', 'dataset'],
    ncols=2,
    use_raw=False,
    layer='scvi_norm_10k_log1p',
    wspace=0.4,
    save='_2_leiden_res_overview.png'
)

rsc.tl.rank_genes_groups_logreg(
    adata,
    groupby='leiden_03',
    use_raw=False,
    n_genes=25,
    layer='scvi_norm_10k_log1p',
)

key = "rank_genes_groups"
groups = list(pd.DataFrame(adata.uns[key]["names"]).columns)

dfs = []
for g in groups:
    dfg = sc.get.rank_genes_groups_df(adata, group=str(g), key=key)
    if "group" not in dfg.columns:
        dfg.insert(0, "group", str(g))
    dfs.append(dfg)

df = pd.concat(dfs, ignore_index=True)
df.to_csv("leiden_03_rank_genes_groups_all.csv", index=False)

sc.pl.rank_genes_groups(
    adata,
    save="_leiden_03_ranked.png"
)

rsc.tl.rank_genes_groups_logreg(
    adata,
    groupby='leiden_04',
    use_raw=False,
    n_genes=25,
    layer='scvi_norm_10k_log1p',
)

key = "rank_genes_groups"
groups = list(pd.DataFrame(adata.uns[key]["names"]).columns)

dfs = []
for g in groups:
    dfg = sc.get.rank_genes_groups_df(adata, group=str(g), key=key)
    if "group" not in dfg.columns:
        dfg.insert(0, "group", str(g))
    dfs.append(dfg)

df = pd.concat(dfs, ignore_index=True)
df.to_csv("rank_genes_groups_all.csv", index=False)

sc.pl.rank_genes_groups(
    adata,
    save="_leiden_04_ranked.png"
)

rsc.tl.embedding_density(
    adata,
    basis='scVI_biolean'
) # ValueError: Cannot find the embedded representation `adata.obsm['X_scvi_biolean']`. Compute the embedding first. --> Could be an easy PR to fix (simple grammatical issue)!
adata.write("real_run/thru_rank_genes_and_leiden.h5ad")

sc.pl.rank_genes_groups_dotplot(
    adata,
    n_genes=3,
    use_raw=False,
    layer="scvi_norm_10k_log1p",
    standard_scale="var",
    cmap="Purples",
    save="_rank_genes_top_3.png"
)
