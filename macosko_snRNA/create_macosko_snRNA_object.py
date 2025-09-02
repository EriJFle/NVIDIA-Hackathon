"""
Download the entire anndata object from the NeMo Google Drive.

gdown --folder 'https://drive.google.com/drive/folders/16XFtWxSfX7CWYzm4XkIXowuom8La83YQ' \
      -O macosko_snrna_drive \
      -c

Then, navigate to that directory and execute the code below in bash.

#!/bin/bash
pigz -dc Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.h5ad.gz > macosko_snRNA.h5ad

"""

"""
Subset a huge .h5ad to astrocytes without loading the full matrix into RAM.

1) Open in backed='r' to read only metadata (cheap).
2) Build a boolean mask for astrocytes (optionally exclude Bergmann glia).
3) Write the slice directly to a new .h5ad (Scanpy/AnnData handles chunking).
"""

import numpy as np
import scanpy as sc
import pandas as pd
from scipy import sparse
import os

os.chdir('/path/to/singlenuclei_data')

A = sc.read_h5ad("macosko_snRNA.h5ad", backed="r")

# 1) Build masks directly from ClusterNm
cn = A.obs["ClusterNm"].astype(str)

# Match labels that with "Astro-" prefix, case doesn't matter givent he argument
is_astro = cn.str.contains(r"^astro", case=False, na=False)

# Exclude Bergmann glia if present
is_berg = cn.str.contains("bergmann", case=False, na=False)
keep = np.flatnonzero(is_astro & ~is_berg)

print(f"Astro nuclei: {keep.size:,}")

# Key: Read subset in chunks and sparsify
chunk_size = 10000  # Adjust based on available memory
n_chunks = (len(keep) + chunk_size - 1) // chunk_size

# Create empty lists to collect data
X_chunks = []
obs_subset = A.obs.iloc[keep].copy()
var_subset = A.var.copy()

print(f"Processing {n_chunks} chunks...")

for i in range(n_chunks):
    start_idx = i * chunk_size
    end_idx = min((i + 1) * chunk_size, len(keep))
    chunk_indices = keep[start_idx:end_idx]

    print(f"Chunk {i+1}/{n_chunks}: indices {start_idx}-{end_idx}")

    # Read chunk and convert to sparse
    X_chunk = A.X[chunk_indices, :]
    if not sparse.issparse(X_chunk):
        # Convert to sparse (most scRNA-seq data should be sparse)
        X_chunk = sparse.csr_matrix(X_chunk)

    X_chunks.append(X_chunk)

# Combine chunks
X_combined = sparse.vstack(X_chunks)

# Create new AnnData object
adata = sc.AnnData(X=X_combined, obs=obs_subset, var=var_subset)

# Copy other attributes if needed
if A.obsm.keys():
    for key in A.obsm.keys():
        adata.obsm[key] = A.obsm[key][keep]

# Add metadata from the methods section of Langlieb et al., 2023
adata.obs['age'] = '8w'
adata.obs['assay'] = '10x_v3_sn'
adata.obs['strain'] = 'C57BL/6J'
adata.obs['tissue_prep'] = 'from_ln_frozen'

# --- Merge Library_Metadata.tsv into adata.obs on derived_cell_libs -> library
lib = pd.read_csv('Library_Metadata.tsv', sep='\t')
lib = lib.drop_duplicates(subset=['library'])

# Ensure the join keys exist and are strings
if 'derived_cell_libs' not in adata.obs.columns:
    # Fallbacks: try an existing 'library' column or parse from obs_names like "BARCODE-LIBID"
    if 'library' in adata.obs.columns:
        adata.obs['derived_cell_libs'] = adata.obs['library'].astype(str)
    elif adata.obs_names.str.contains('-').any():
        adata.obs['derived_cell_libs'] = adata.obs_names.str.rsplit('-', n=1).str[-1]
    else:
        raise KeyError("Missing 'derived_cell_libs' in adata.obs and unable to infer it.")


adata.obs['derived_cell_libs'] = adata.obs['derived_cell_libs'].astype(str)

lib['library'] = lib['library'].astype(str)

# Prefix all lib columns except the key for clarity, then merge
lib_pref = lib.rename(columns={c: (f'{c}' if c != 'library' else c) for c in lib.columns})
adata.obs = adata.obs.merge(
    lib_pref,
    how='left',
    left_on='derived_cell_libs',
    right_on='library',
    validate='many_to_one'
)

del adata.obs['library']

# --- Merge CellType_Metadata.tsv into adata.obs on ClusterNm -> Annotation
cell_lib = pd.read_csv('CellType_Metadata.tsv', sep='\t')
cell_lib = cell_lib.drop_duplicates(subset=['Annotation'])
adata.obs['ClusterNm'] = adata.obs['ClusterNm'].astype(str)
cell_lib['Annotation'] = cell_lib['Annotation'].astype(str)

# Prefix all cell-type columns except the key, then merge
cell_pref = cell_lib[['Annotation', 'Max_TopStruct', 'Max_DeepCCF']].copy()
adata.obs = adata.obs.merge(
    cell_pref,
    how='left',
    left_on='ClusterNm',
    right_on='Annotation',
    validate='many_to_one'
)
# Disambiguate the right-side key
del adata.obs['Annotation']

# --- HDF5-safe: ensure all obs columns are scalar-friendly (strings/categoricals) ---
import numpy as _np
import pandas as _pd
from collections.abc import Sequence as _Seq

def _to_scalar_string(x):
    """Convert lists/dicts/mixed objects into a string safe for HDF5 vlen."""
    if _pd.isna(x):
        return ""
    if isinstance(x, (str, bytes)):
        return str(x)
    if isinstance(x, (int, float, bool, _np.integer, _np.floating, _np.bool_)):
        return str(x)
    # datetimes with isoformat
    if hasattr(x, "isoformat"):
        try:
            return x.isoformat()
        except Exception:
            pass
    # list/tuple/set -> pipe-joined
    if isinstance(x, _Seq) and not isinstance(x, (str, bytes)):
        return "|".join(str(y) for y in x)
    # dict -> key:value pairs
    if isinstance(x, dict):
        return "|".join(f"{k}:{v}" for k, v in x.items())
    return str(x)

# Fix the known offender first (if present)
if "removal_reason" in adata.obs.columns:
    adata.obs["removal_reason"] = adata.obs["removal_reason"].map(_to_scalar_string).astype("category")

# Coerce all object-dtype obs columns to HDF5-safe strings
for _col in adata.obs.select_dtypes(include=["object"]).columns:
    adata.obs[_col] = adata.obs[_col].map(_to_scalar_string)

# Write with lossless compression
adata.write_h5ad("/path/to/AstroWRLD/singlenuclei_data/macosko_astro_only.h5ad", compression="lzf")
