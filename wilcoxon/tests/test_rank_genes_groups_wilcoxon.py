import numpy as np
import scanpy as sc
import pytest

# Skip entire module if GPU deps not available
pytest.importorskip("cupy")
pytest.importorskip("cudf")

from src.wilcoxon.tl._rank_genes_groups_wilcoxon import rank_genes_groups

def _load_pbmc3k_or_synth():
    try:
        return sc.datasets.pbmc3k()
    except Exception:
        rng = np.random.default_rng(0)
        X = rng.negative_binomial(2, 0.2, size=(6000, 2000)).astype(np.float32)
        ad = sc.AnnData(X)
        ad.obs["leiden"] = rng.integers(0, 8, size=ad.n_obs).astype(str)
        return ad


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_gpu_wilcoxon_parity_small():
    ad = _load_pbmc3k_or_synth()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    if "leiden" not in ad.obs:
        sc.pp.pca(ad)
        sc.pp.neighbors(ad)
        sc.tl.leiden(ad, resolution=0.5, key_added="leiden")

    sc.tl.rank_genes_groups(ad, "leiden", method="wilcoxon", tie_correct=True, n_genes=30, key_added="cpu")
    rank_genes_groups(ad, "leiden", n_genes=30)

    cpu = ad.uns["cpu"]; gpu = ad.uns["rank_genes_groups"]
    for g in gpu["params"]["groups"]:
        top_cpu = list(cpu["names"][g]); top_gpu = list(gpu["names"][g])
        overlap = len(set(top_cpu) & set(top_gpu))
        assert overlap >= 20, f"Low overlap for group {g}: {overlap}/30"

        zc = np.array(cpu["scores"][g], dtype=float)
        zg = np.array(gpu["scores"][g], dtype=float)
        rmse = np.sqrt(np.mean((zc - zg) ** 2))
        assert rmse < 0.35, f"High z RMSE for group {g}: {rmse}"

def test_gpu_wilcoxon_degenerate_group():
    ad = _load_pbmc3k_or_synth()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    ad.obs["deg"] = "A"
    ad.obs["deg"].iloc[:100] = "B"  # two groups; small A, large rest
    # should not crash
    rank_genes_groups(ad, "deg", n_genes=10)
