import scanpy as sc
import numpy as np

from _rank_genes_groups_wilcoxon import rank_genes_groups_wilcoxon

print(sc.__version__)
print(rsc.__version__)
print(cp.__version__)

'''
1.11.4
0.13.1
13.6.0
'''

sc.set_figure_params(dpi=300)
sc.settings.verbosity = 3
sc.settings.figdir = "nvidia-hackathon/figures"

adata = sc.read_h5ad("gpu_de_ranked.h5ad")

print(adata.shape)

'''
(737011, 3259)
'''

del adata.uns['rank_genes_groups']

adata.X = adata.layers['scvi_norm_10k_log1p'].copy()

import time

t = time.perf_counter()

rank_genes_groups_wilcoxon(
    adata,
    groupby="region_broad",
    use_raw=False,
    tie_correct=False,
)

print(f"gpu_rank_genes_groups_wilcoxon took {time.perf_counter()-t:.2f}s")

'''
gpu_rank_genes_groups_wilcoxon took 21.47s
'''


sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby="region_broad",
    standard_scale="var",
    n_genes=3,
    use_raw=False,
    layer=None,
    cmap="Purples",
    save="dotplot_gpu_wilcoxon_region_ranked.png"
)

del t

t=time.perf_counter()

sc.tl.rank_genes_groups(
    adata,
    groupby="region_broad",
    groups="all",
    use_raw=False,
    tie_correct=False,
    method="wilcoxon",
    layer=None,
    key_added="cpu_ranks_final",
)
print(f"cpu_rank_genes_groups_wilcoxon took {time.perf_counter()-t:.2f}s")

'''
cpu_rank_genes_groups_wilcoxon took 222.50s
'''

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby="region_broad",
    standard_scale="var",
    n_genes=3,
    use_raw=False,
    layer=None,
    cmap="Purples",
    key="cpu_ranks_final",
    save="dotplot_cpu_wilcoxon_region_ranked.png"
)

cpu_res = adata.uns["cpu_ranks_final"].copy()
gpu_res = adata.uns["rank_genes_groups"].copy()

cpu_scores = cpu_res["scores"]
gpu_scores = gpu_res["scores"]

for group in cpu_scores.dtype.names:
    a = cpu_scores[group]
    b = gpu_scores[group]
    np.testing.assert_array_equal(a, b)
    identical = np.array_equal(a, b)
    print(f"{group}: {'match' if identical else 'DIFFER'}")

'''
CB: match
CTXsp: match
HPF: match
HY: match
Isocortex: match
MB: match
MY: match
OLF: match
P: match
PAL: match
STR: match
TH: match
'''

# Test with Larger Gene Pool (6x Larger)
adata = adata.raw.to_adata()
adata.shape

'''
(737011, 18982)
'''

del t

t = time.perf_counter()

rank_genes_groups_wilcoxon(
    adata,
    groupby="region_broad",
    tie_correct=False,
)

print(f"gpu_rank_genes_groups_wilcoxon took {time.perf_counter()-t:.2f}s")

# LogReg Speed test
del t

t = time.perf_counter()

rsc.tl.rank_genes_groups_logreg(
    adata,
    groupby="region_broad",
    use_raw=False,
)

print(f"gpu_rank_genes_groups_region_logreg took {time.perf_counter()-t:.2f}s")

'''
gpu_rank_genes_groups_region_logreg took 82.33s
'''

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby="region_broad",
    standard_scale="var",
    n_genes=3,
    use_raw=False,
    layer=None,
    cmap="Purples",
    save="logreg_gpu_region_ranked.png"
)

# GPU wilcoxon vs logreg on leiden

del t

t = time.perf_counter()

rank_genes_groups_wilcoxon(
    adata,
    groupby="leiden",
    use_raw=False,
    tie_correct=False,
)

print(f"gpu_leiden_rank_genes_groups_wilcoxon took {time.perf_counter()-t:.2f}s")

'''
gpu_leiden_rank_genes_groups_wilcoxon took 41.2s
'''

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby="leiden",
    standard_scale="var",
    n_genes=3,
    use_raw=False,
    layer=None,
    cmap="Purples",
    save="wilcoxon_leiden_gpu_ranked.png"
)

del t

t = time.perf_counter()

rsc.tl.rank_genes_groups_logreg(
    adata,
    groupby="leiden",
    use_raw=False,
)

print(f"gpu_rank_genes_groups_leiden_logreg took {time.perf_counter()-t:.2f}s")

'''
gpu_rank_genes_groups_leiden_logreg took 68.10s
'''

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby="leiden",
    standard_scale="var",
    n_genes=3,
    use_raw=False,
    layer=None,
    cmap="Purples",
    save="logreg_leiden_gpu_ranked.png"
)
