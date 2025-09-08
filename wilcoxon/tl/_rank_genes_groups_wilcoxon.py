from __future__ import annotations

import numpy as np
from math import ceil

def _require_gpu():
    try:
        import cupy as cp  # noqa: F401
        import cudf  # noqa: F401
    except Exception as e:
        raise NotImplementedError("GPU Wilcoxon requires CuPy + cuDF.") from e

def _init_rmm_if_available():
    # Optional: use RMM pool for stability/speed; safe if RMM not installed
    try:
        import rmm
        import cupy as cp
        rmm.reinitialize(pool_allocator=True)
        cp.cuda.set_allocator(rmm.rmm_cupy_allocator)
    except Exception:
        pass

def rank_genes_groups(
    adata,
    groupby: str,
    *,
    groups=None,
    reference: str = "rest",
    layer: str | None = None,
    n_genes: int = 100,
    tie_correct: bool = True,
    use_raw: bool | None = None,      # accepted but unused; kept for Scanpy parity
    batch_genes: int | None = None,
    rank_dtype: str = "float32",      # ranks dtype; values ranked in float64
    zero_fastpath: bool = True,
    continuity_correction: bool = False,
    corr_method: str = "benjamini-hochberg",
    copy: bool = False,
):
    """
    GPU-accelerated Wilcoxon rank-sum (1-vs-rest) with tie-corrected variance.
    Returns Scanpy-compatible results in `adata.uns["rank_genes_groups"]`.
    """
    _require_gpu()
    _init_rmm_if_available()

    import cupy as cp
    import cudf
    import scanpy as sc

    if reference != "rest":
        # v0 scope: only 1-vs-rest
        raise NotImplementedError("Only reference='rest' is supported in GPU Wilcoxon v0.")

    if groupby not in adata.obs:
        raise KeyError(f"{groupby!r} not found in adata.obs")

    # ---------- pull X as a reader that yields (cells x G_batch) cudf.DataFrame ----------
    X = adata.layers[layer] if layer is not None else adata.X
    n_cells, n_genes_total = X.shape
    varnames = np.asarray(adata.var_names)

    def _read_gene_batch(col_idx: np.ndarray) -> "cudf.DataFrame":
        # host -> gpu dense (float64) to minimize spurious ties
        if hasattr(X, "tocsr"):  # sparse anndata
            sub = X[:, col_idx].toarray()
        else:
            sub = np.asarray(X[:, col_idx])
        df = cudf.DataFrame.from_numpy(sub.astype(np.float64, copy=False))
        df.columns = adata.var_names[col_idx].to_numpy()
        # ensure no NaN: Wilcoxon expects real numbers (scRNA should be NaN-free)
        df = df.fillna(0.0)
        return df

    # ---------- groups / masks ----------
    cat = adata.obs[groupby].astype("category")
    adata.obs[groupby] = cat
    all_groups = list(cat.cat.categories)
    if groups is None or groups == "all":
        groups = all_groups
    groups = [g for g in groups if g in all_groups]
    codes_host = cat.cat.codes.values
    codes = cudf.Series(codes_host)  # GPU
    code_map = {g: i for i, g in enumerate(all_groups)}

    # ---------- batch size heuristic ----------
    def _auto_batch_genes(n_cells: int, user: int | None) -> int:
        if user is not None:
            return int(user)
        # aim ~1–1.5 GB working set per batch (double-buffer-ish)
        # bytes ≈ cells * genes * 8 (float64) * ~1.2 (overhead)
        # => genes ≈ (1.2e9 / 8) / cells
        G = max(32, min(2048, int(1.2e9 / 8 / max(n_cells, 1))))
        return G

    G = _auto_batch_genes(n_cells, batch_genes)
    rank_dtype_np = np.float32 if rank_dtype == "float32" else np.float64

    # ---------- accumulators ----------
    out = {
        g: {
            "z": cp.full(n_genes_total, cp.nan, dtype=cp.float32),
            "p": cp.full(n_genes_total, 1.0, dtype=cp.float32),
            "auc": cp.full(n_genes_total, 0.5, dtype=cp.float32),
            "logfc": cp.full(n_genes_total, 0.0, dtype=cp.float32),
        }
        for g in groups
    }

    # helpers
    from cupyx.scipy.special import ndtr

    def _bh_fdr(p_vec: "cp.ndarray") -> "cp.ndarray":
        p = p_vec.astype(cp.float64, copy=False)
        n = p.size
        order = cp.argsort(p)
        ranked = p[order]
        q = ranked * n / (cp.arange(1, n + 1))
        q = cp.minimum.accumulate(q[::-1])[::-1]
        out = cp.empty_like(q)
        out[order] = cp.minimum(q, 1.0)
        return out.astype(cp.float32, copy=False)

    # ---------- iterate gene batches ----------
    n_batches = ceil(n_genes_total / G)
    for b in range(n_batches):
        lo = b * G
        hi = min((b + 1) * G, n_genes_total)
        idx = np.arange(lo, hi, dtype=np.int64)
        V = _read_gene_batch(idx)  # cudf.DataFrame, (cells x G_b), float64
        G_b = len(V.columns)
        if G_b == 0:
            continue

        # Optional fast path: all-zero columns
        if zero_fastpath:
            any_nonzero = V.any(axis=0)  # cudf Series bool per column
            if (~any_nonzero).any():
                # drop zero-only columns from further work
                V = V.loc[:, V.columns[any_nonzero.to_pandas().values]]
                if len(V.columns) == 0:
                    continue

        # Column-wise ranks with average ties on pooled sample
        R = V.rank(method="average", axis=0, na_option="keep")  # cudf.DataFrame
        if rank_dtype_np is np.float32:
            R = R.astype(np.float32)

        # Precompute tie term T per column on pooled values:
        # T = sum_k (t_k^3 - t_k), t_k are tie block sizes among pooled (A∪rest).
        # Using value_counts over V per column.
        T_list = []
        for c in V.columns:
            vc = V[c].value_counts()         # cudf Series (value -> count)
            t = vc.values.to_cupy()          # cupy array
            T_col = (t ** 3 - t).sum()
            T_list.append(T_col)
        T = cp.stack(T_list)                  # shape (G_eff,)

        # Pooled means for rest computation
        mean_all = cp.asarray(V.mean(axis=0).to_cupy(), dtype=cp.float64)

        N = len(codes)
        for g in groups:
            code = code_map[g]
            is_A = (codes == code)
            n1 = int(is_A.sum())
            n2 = N - n1
            if n1 == 0 or n2 == 0:
                # degenerate; stats remain default
                continue
            n1n2 = n1 * n2

            # rank-sum for group A
            rank_sum = cp.asarray(R.loc[is_A].sum(axis=0).to_cupy(), dtype=cp.float64)
            U = rank_sum - (n1 * (n1 + 1) / 2.0)
            auc = (U / n1n2).astype(cp.float32)

            # tie-corrected variance
            # Var(U) = n1*n2/12 * [ (N+1) - T/(N*(N-1)) ]
            if tie_correct:
                corr = (N + 1) - (T / (N * (N - 1)))
            else:
                corr = (N + 1) * cp.ones_like(T)
            varU = (n1n2 / 12.0) * corr
            sdU = cp.sqrt(varU)

            # z and two-sided p (normal approx), optional continuity correction
            center = n1n2 / 2.0
            if continuity_correction:
                cc = 0.5 * cp.sign(U - center)
                z = (U - center - cc) / sdU
            else:
                z = (U - center) / sdU
            p = 2.0 * (1.0 - ndtr(cp.abs(z)))
            z = z.astype(cp.float32)
            p = p.astype(cp.float32)

            # log fold-change: log(mean(A)+eps) - log(mean(rest)+eps)
            eps = 1e-9
            mean_A = cp.asarray(V.loc[is_A].mean(axis=0).to_cupy(), dtype=cp.float64)
            mean_rest = (N * mean_all - n1 * mean_A) / max(n2, 1)
            logfc = (cp.log(mean_A + eps) - cp.log(mean_rest + eps)).astype(cp.float32)

            # write into global slots (note: V may be a subset if zero-fastpath pruned)
            col_names = np.array(V.columns.to_pandas())
            # map local cols -> global indices in [lo:hi)
            local_map = {name: lo + i for i, name in enumerate(varnames[lo:hi])}
            # The above assumes V kept column order; safer mapping:
            global_idx = np.array([np.where(varnames == name)[0][0] for name in col_names], dtype=np.int64)

            out[g]["z"][global_idx] = z
            out[g]["p"][global_idx] = p
            out[g]["auc"][global_idx] = auc
            out[g]["logfc"][global_idx] = logfc

        # free big frames early
        del V, R, T, mean_all

    # ---------- assemble Scanpy-compatible outputs ----------
    results = {
        "params": {
            "groupby": groupby,
            "groups": groups,
            "reference": "rest",
            "method": "wilcoxon_gpu",
            "tie_correct": tie_correct,
            "continuity_correction": continuity_correction,
            "n_genes": n_genes,
            "layer": layer,
            "rank_dtype": rank_dtype,
            "batch_genes": G,
        },
        "names": {},
        "scores": {},
        "pvals": {},
        "pvals_adj": {},
        "logfoldchanges": {},
        "auc": {},
    }

    for g in groups:
        z = out[g]["z"].get()
        p = out[g]["p"]
        q = _bh_fdr(p).get()
        order = np.argsort(-z)  # descending z
        top = order[:n_genes]
        results["names"][g] = varnames[top]
        results["scores"][g] = z[top]
        results["pvals"][g] = p.get()[top]
        results["pvals_adj"][g] = q[top]
        results["logfoldchanges"][g] = out[g]["logfc"].get()[top]
        results["auc"][g] = out[g]["auc"].get()[top]

    if copy:
        return results

    adata.uns["rank_genes_groups"] = results
    return None