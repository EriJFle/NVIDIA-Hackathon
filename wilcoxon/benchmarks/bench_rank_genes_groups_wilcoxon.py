import time
import argparse
import numpy as np
import scanpy as sc

def synth_adata(n_cells, n_genes, n_groups):
    rng = np.random.default_rng(42)
    X = rng.negative_binomial(2, 0.2, size=(n_cells, n_genes)).astype(np.float32)
    ad = sc.AnnData(X)
    ad.obs["grp"] = rng.integers(0, n_groups, size=n_cells).astype(str)
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    return ad

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cells", type=int, default=100000)
    parser.add_argument("--genes", type=int, default=2000)
    parser.add_argument("--groups", type=int, default=10)
    parser.add_argument("--batch-genes", type=int, default=None)
    args = parser.parse_args()

    ad = synth_adata(args.cells, args.genes, args.groups)

    # CPU baseline (subset to reasonable size to avoid hours)
    ad_cpu = ad[:, :min(2000, args.genes)].copy()
    t0 = time.time()
    sc.tl.rank_genes_groups(ad_cpu, "grp", method="wilcoxon", tie_correct=True, n_genes=50)
    t_cpu = time.time() - t0
    print(f"CPU (Scanpy) on {ad_cpu.n_obs}x{ad_cpu.n_vars}: {t_cpu:.2f}s")

    # GPU full (local import from repo)
    from src.wilcoxon.tl._rank_genes_groups_wilcoxon import rank_genes_groups as gpu_rg
    t0 = time.time()
    gpu_rg(ad, "grp", n_genes=50, batch_genes=args.batch_genes)
    t_gpu = time.time() - t0
    print(f"GPU (wilcoxon_gpu) on {ad.n_obs}x{ad.n_vars}: {t_gpu:.2f}s")

if __name__ == "__main__":
    main()
