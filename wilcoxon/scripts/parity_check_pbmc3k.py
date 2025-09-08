import scanpy as sc
import numpy as np
from src.wilcoxon.tl._rank_genes_groups_wilcoxon import rank_genes_groups

def main():
    ad = sc.datasets.pbmc3k()  # small, easy
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    sc.tl.leiden(ad, resolution=0.5, key_added="leiden")

    # CPU baseline
    sc.tl.rank_genes_groups(ad, "leiden", method="wilcoxon", tie_correct=True, n_genes=50, key_added="cpu")

    # GPU
    rank_genes_groups(ad, "leiden", n_genes=50)

    cpu = ad.uns["cpu"]
    gpu = ad.uns["rank_genes_groups"]
    groups = gpu["params"]["groups"]
    for g in groups:
        top_cpu = list(cpu["names"][g])
        top_gpu = list(gpu["names"][g])
        overlap = len(set(top_cpu) & set(top_gpu))
        print(f"{g}: Top50 overlap={overlap}/50")
        z_cpu = np.array(cpu["scores"][g], dtype=float)
        z_gpu = np.array(gpu["scores"][g], dtype=float)
        print(f"  Δz RMSE: {np.sqrt(np.mean((z_cpu - z_gpu)**2)):.3f}")

if __name__ == "__main__":
    main()
