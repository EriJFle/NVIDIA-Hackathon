from pathlib import Path
import pandas as pd
import numpy as np
import anndata as ad
from tqdm import tqdm
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

# =========================
# CONFIG
# =========================
CACHE_DIR = Path("/path/to/allen_wmb/abc_cache")  # put this on fast scratch
OUTDIR    = Path("/path/to/allen_wmb/outputs")
OUTDIR.mkdir(parents=True, exist_ok=True)
OUT_H5AD  = OUTDIR / "Allen_WMB10X_v2v3_Astrocytes_raw.h5ad"

# Process BOTH datasets
EXPR_DIRS = ["WMB-10Xv2", "WMB-10Xv3"]

# =========================
# INIT CACHE
# =========================
print(f"Using cache: {CACHE_DIR}")
abc = AbcProjectCache.from_cache_dir(CACHE_DIR)
try:
    abc.load_latest_manifest()
except Exception:
    pass
print("Manifest:", abc.current_manifest)

# =========================
# LOAD MOUSE METADATA & TAXONOMY
# =========================
# Cell metadata (covers both v2 and v3)
cell = abc.get_metadata_dataframe(
    directory="WMB-10X",
    file_name="cell_metadata",
    dtype={"cell_label": str},
).set_index("cell_label")

# Taxonomy (class/subclass/supertype/cluster + neurotransmitter)
term_set   = abc.get_metadata_dataframe("WMB-taxonomy", "cluster_annotation_term_set")
membership = abc.get_metadata_dataframe("WMB-taxonomy", "cluster_to_cluster_annotation_membership")

pivot = (
    membership.groupby(["cluster_alias", "cluster_annotation_term_set_name"])["cluster_annotation_term_name"]
    .first()
    .unstack()
)

# Identify astro classes (robust: case-insensitive contains)
if "class" not in pivot.columns:
    raise RuntimeError("Taxonomy pivot missing 'class' column — check WMB-taxonomy inputs.")
astro_cluster_aliases = set(pivot.index[pivot["class"].str.contains("astro", case=False, na=False)])
print("Distinct astrocyte clusters (class level):", len(astro_cluster_aliases))

tax_cols = [c for c in ["neurotransmitter", "class", "subclass", "supertype", "cluster"] if c in pivot.columns]

# =========================
# HELPER: add safe metadata to obs
# =========================
def add_obs_columns_safe(adata_pkg: ad.AnnData, meta: pd.DataFrame, cols):
    cols = [c for c in cols if c in meta.columns]
    # avoid overwriting existing columns in adata_pkg.obs
    to_add = [c for c in cols if c not in adata_pkg.obs.columns]
    if to_add:
        adata_pkg.obs = adata_pkg.obs.join(meta[to_add], how="left")

# =========================
# SCAN, SLICE, CONCAT
# =========================
adatas = []
for EXPR_DIR in EXPR_DIRS:
    # Restrict to the dataset version (v2 or v3)
    cell_ds = cell[cell["dataset_label"] == EXPR_DIR]
    if cell_ds.empty:
        print(f"[{EXPR_DIR}] No cells in metadata; skipping.")
        continue

    # Astro cells in this dataset
    astro_cells = cell_ds[cell_ds["cluster_alias"].isin(astro_cluster_aliases)].copy()
    print(f"[{EXPR_DIR}] Astrocyte cells in metadata: {len(astro_cells):,}")
    if astro_cells.empty:
        continue

    # Join taxonomy names by cluster_alias for convenience
    tax_map = pivot.loc[list(astro_cluster_aliases), tax_cols].copy()
    astro_cells = astro_cells.join(tax_map, on="cluster_alias", how="left")

    # Discover RAW expression packages for this dataset
    files = abc.list_data_files(EXPR_DIR)  # e.g., 'WMB-10Xv2-TH/raw', '.../log2'
    raw_entries = sorted([f for f in files if f.endswith("/raw")])
    print(f"[{EXPR_DIR}] RAW packages discovered: {len(raw_entries)}")

    for entry in tqdm(raw_entries, desc=f"Packages {EXPR_DIR}"):
        pkg_label = entry.split("/")[0]  # e.g., 'WMB-10Xv2-TH'

        # Cells assigned to this package
        pkg_cells = astro_cells.index[astro_cells["feature_matrix_label"] == pkg_label]
        if len(pkg_cells) == 0:
            continue

        fpath = abc.get_data_path(directory=EXPR_DIR, file_name=entry)
        ad_backed = ad.read_h5ad(fpath, backed="r")

        # Intersect defensively with what's in the file
        pkg_cells = pd.Index(pkg_cells).intersection(ad_backed.obs_names)
        if len(pkg_cells) == 0:
            ad_backed.file.close()
            continue

        # Slice to memory
        ad_pkg = ad_backed[pkg_cells, :].to_memory()
        ad_backed.file.close()

        # Add dataset/package tags
        ad_pkg.obs["dataset_version"] = EXPR_DIR
        ad_pkg.obs["wmb_package"] = pkg_label

        # Add useful metadata safely (avoid collisions)
        candidate_cols = [
            "dataset_label", "feature_matrix_label",
            "region_of_interest_acronym", "brain_section_label",
            "donor_label", "donor_genotype", "donor_sex", "cluster_alias"
        ] + tax_cols
        add_obs_columns_safe(ad_pkg, astro_cells.loc[pkg_cells], candidate_cols)

        # QC per chunk (cheap on sparse)
        ad_pkg.obs["n_counts"] = np.asarray(ad_pkg.X.sum(axis=1)).ravel()
        ad_pkg.obs["n_genes"]  = np.asarray((ad_pkg.X > 0).sum(axis=1)).ravel()

        adatas.append(ad_pkg)

# Final concat across BOTH datasets
if not adatas:
    raise RuntimeError("No astrocyte cells matched in v2 or v3 RAW packages.")
combined = ad.concat(
    adatas,
    axis=0,
    join="outer",         # union of genes across v2 and v3
    merge="unique",
    label="batch",
    index_unique=None,    # keep original cell labels
)

combined.uns["source_datasets"] = EXPR_DIRS
combined.uns["value_type"]      = "raw"
combined.uns["note"]            = "Allen WMB-10X mouse atlas; astrocytes from v2+v3 RAW matrices"

print("Final AnnData shape:", combined.n_obs, "cells x", combined.n_vars, "genes")
print("Writing:", OUT_H5AD)
combined.write(OUT_H5AD, compression="lzf")
print("Done.")
