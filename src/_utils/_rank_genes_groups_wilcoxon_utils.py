from __future__ import annotations

from typing import Any

import cupy as cp
import cupyx.scipy.special
import numpy as np
import pandas as pd


EPS = 1e-9


def _norm_params(
    *,
    groupby,
    groups,
    reference,
    layer,
    use_raw,
    n_top,
    method_kwargs,
):
    """Normalize Wilcoxon parameters while preserving Scanpy semantics."""
    if method_kwargs:
        method_kwargs.clear()
    if layer and use_raw is True:
        raise ValueError("Cannot specify `layer` and have `use_raw=True`.")
    if reference != "rest":
        raise NotImplementedError(
            "GPU Wilcoxon implementation currently only supports reference='rest'. "
            "For vs-specific-group comparisons, please use scanpy.tl.rank_genes_groups"
        )
    if groups in ("all", None):
        normalized_groups = None
    elif isinstance(groups, (str, int)):
        raise ValueError("Specify a sequence of groups")
    else:
        normalized_groups = [str(g) for g in groups]
    return groupby, normalized_groups, reference, layer, use_raw, n_top


def _get_X_obs_var(
    adata: Any,
    *,
    groupby: str,
    layer: str | None,
    use_raw: bool | None,
):
    """Return expression matrix, categorical labels, and feature index."""
    labels = pd.Series(adata.obs[groupby]).reset_index(drop=True)
    if not pd.api.types.is_categorical_dtype(labels):
        labels = labels.astype("category")
    if layer:
        X = adata.layers[layer]
        var_names = adata.var_names
    elif use_raw is None and hasattr(adata, "raw") and adata.raw:
        X = adata.raw.X
        var_names = adata.raw.var_names
    elif use_raw is True:
        X = adata.raw.X
        var_names = adata.raw.var_names
    else:
        X = adata.X
        var_names = adata.var_names
    return X, labels, var_names


def _group_masks(labels: pd.Series, groups: list[str] | None):
    """Materialize group order, observation masks, and group sizes."""
    categories = labels.cat.categories
    codes = labels.cat.codes.to_numpy()
    cat_strings = np.array([str(cat) for cat in categories])
    if groups is None:
        idx = np.arange(len(categories), dtype=int)
    else:
        idx = []
        for name in groups:
            match = np.where(cat_strings == name)[0]
            if match.size == 0:
                raise ValueError(f"Group {name!r} not found in observation categories.")
            idx.append(int(match[0]))
        idx = np.asarray(idx, dtype=int)
    order = categories[idx].to_numpy()
    masks = codes[None, :] == idx[:, None]
    group_sizes = masks.sum(axis=1)
    n_cells = labels.shape[0]
    for group, size in zip(order, group_sizes):
        rest = n_cells - size
        if size <= 25 or rest <= 25:
            print(
                f"Warning: Group {group} has ≤25 cells. "
                "Normal approximation may be less accurate."
            )
    return order, masks, group_sizes


def _rank_chunks(X: Any, *, block: int):
    """Yield CuPy rank arrays and raw chunks in gene blocks."""
    n_genes = X.shape[1]
    for left in range(0, n_genes, block):
        right = min(left + block, n_genes)
        chunk = X[:, left:right]
        if hasattr(chunk, "toarray"):
            chunk = chunk.toarray()
        chunk_dev = cp.asarray(chunk, dtype=cp.float64)
        ranks = cp.empty_like(chunk_dev, dtype=cp.float64)
        for j in range(chunk_dev.shape[1]):
            column = chunk_dev[:, j].ravel()
            sorter = cp.argsort(column)
            sorted_column = column[sorter]
            unique = cp.concatenate(
                (cp.array([True]), sorted_column[1:] != sorted_column[:-1])
            )
            dense = cp.empty(unique.size, dtype=cp.int64)
            dense[sorter] = cp.cumsum(unique)
            bounds = cp.concatenate((cp.flatnonzero(unique), cp.array([unique.size])))
            ranks[:, j] = 0.5 * (bounds[dense] + bounds[dense - 1] + 1.0)
        yield ranks, chunk_dev, slice(left, right)


def _wilcoxon(
    ranks: cp.ndarray,
    group_matrix: cp.ndarray,
    group_sizes: cp.ndarray,
    *,
    tie_correct: bool,
    n_cells: int,
):
    """Project ranks to group sums and compute z-scores and p-values."""

    def _tie_corr(vals: cp.ndarray) -> cp.ndarray:
        corr = cp.ones(vals.shape[1], dtype=cp.float64)
        for j in range(vals.shape[1]):
            arr = cp.sort(vals[:, j].ravel())
            idx = cp.flatnonzero(
                cp.concatenate((cp.array([True]), arr[1:] != arr[:-1], cp.array([True])))
            )
            diff = cp.diff(idx).astype(cp.float64)
            size = cp.float64(arr.size)
            if size >= 2:
                corr[j] = 1.0 - (diff**3 - diff).sum() / (size**3 - size)
        return corr

    tie = _tie_corr(ranks) if tie_correct else cp.ones(ranks.shape[1], dtype=cp.float64)
    rank_sums = group_matrix.T @ ranks
    expected = group_sizes[:, None] * (n_cells + 1) / 2.0
    rest_sizes = n_cells - group_sizes
    variance = tie[None, :] * group_sizes[:, None] * rest_sizes[:, None]
    variance *= (n_cells + 1) / 12.0
    std = cp.sqrt(variance)
    z = (rank_sums - expected) / std
    cp.nan_to_num(z, copy=False)
    abs_z = cp.abs(z)
    p = 2.0 * (1.0 - cupyx.scipy.special.ndtr(abs_z))
    return z, p


def _effects(
    X_chunk: cp.ndarray,
    group_matrix: cp.ndarray,
    group_sizes: cp.ndarray,
    *,
    expm1_func,
    eps: float = EPS,
):
    """Mirror Scanpy fold-change computation on the device for one chunk."""
    group_sums = group_matrix.T @ X_chunk
    group_means = group_sums / group_sizes[:, None]
    n_cells = X_chunk.shape[0]
    total_mean = cp.mean(X_chunk, axis=0)
    rest_sizes = n_cells - group_sizes
    rest_sum = total_mean * n_cells - group_means * group_sizes[:, None]
    mean_rest = rest_sum / rest_sizes[:, None]
    numerator = expm1_func(group_means) + eps
    denominator = expm1_func(mean_rest) + eps
    return cp.log2(numerator / denominator)


def _adjust_pvals(*, pvals: np.ndarray, method: str, total_genes: int) -> tuple[np.ndarray, np.ndarray]:
    """Benjamini-Hochberg or Bonferroni correction on host values."""
    cleaned = np.array(pvals, copy=True)
    cleaned[np.isnan(cleaned)] = 1.0
    if method == "benjamini-hochberg":
        from statsmodels.stats.multitest import multipletests

        _, adjusted, _, _ = multipletests(cleaned, alpha=0.05, method="fdr_bh")
        return cleaned, adjusted
    return cleaned, np.minimum(cleaned * total_genes, 1.0)


def _pack_result(
    var_names,
    groups,
    *,
    scores,
    logfc,
    pvals,
    gene_indices,
    n_top: int,
    method: str,
    total_genes: int,
    params: dict,
):
    """Assemble the Scanpy-compatible rank genes result."""
    ordered: dict[str, dict[str, np.ndarray]] = {}
    var_array = np.asarray(var_names)
    for group in groups:
        score_arr = np.concatenate(scores[group]) if scores[group] else np.empty(0)
        logfc_arr = np.concatenate(logfc[group]) if logfc[group] else np.empty(0)
        pval_arr = np.concatenate(pvals[group]) if pvals[group] else np.empty(0)
        idx_arr = np.concatenate(gene_indices[group]) if gene_indices[group] else np.empty(0, dtype=int)
        if pval_arr.size:
            clean, adj = _adjust_pvals(pvals=pval_arr, method=method, total_genes=total_genes)
        else:
            clean = np.array([], dtype=float)
            adj = np.array([], dtype=float)
        sort_keys = np.column_stack((score_arr, -idx_arr)) if score_arr.size else np.empty((0, 2))
        sort_idx = np.lexsort((sort_keys[:, 1], sort_keys[:, 0]))[::-1] if sort_keys.size else np.empty(0, dtype=int)
        limit = min(n_top, sort_idx.size)
        keep = sort_idx[:limit]
        ordered[group] = {
            "scores": score_arr[keep],
            "logfc": logfc_arr[keep],
            "pvals": clean[keep],
            "pvals_adj": adj[keep],
            "names": var_array[idx_arr[keep]],
        }
    groups_order = list(ordered.keys())
    scores_rec = np.rec.fromarrays(
        [ordered[g]["scores"].astype(np.float32) for g in groups_order],
        dtype=[(g, "float32") for g in groups_order],
    )
    names_rec = np.rec.fromarrays(
        [ordered[g]["names"].astype("U50") for g in groups_order],
        dtype=[(g, "U50") for g in groups_order],
    )
    logfc_rec = np.rec.fromarrays(
        [ordered[g]["logfc"].astype(np.float32) for g in groups_order],
        dtype=[(g, "float32") for g in groups_order],
    )
    pvals_rec = np.rec.fromarrays(
        [ordered[g]["pvals"].astype(np.float64) for g in groups_order],
        dtype=[(g, "float64") for g in groups_order],
    )
    adj_rec = np.rec.fromarrays(
        [ordered[g]["pvals_adj"].astype(np.float64) for g in groups_order],
        dtype=[(g, "float64") for g in groups_order],
    )
    return {
        "params": params,
        "scores": scores_rec,
        "names": names_rec,
        "logfoldchanges": logfc_rec,
        "pvals": pvals_rec,
        "pvals_adj": adj_rec,
    }
