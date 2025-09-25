from __future__ import annotations

from typing import TYPE_CHECKING, Iterable, Literal

import cupy as cp
import numpy as np

from _rank_genes_groups_wilcoxon_helpers import (
    _effects,
    _get_X_obs_var,
    _group_masks,
    _norm_params,
    _pack_result,
    _rank_chunks,
    _wilcoxon,
)

if TYPE_CHECKING:
    from anndata import AnnData


def rank_genes_groups_wilcoxon(
    adata: "AnnData",
    groupby: str,
    *,
    groups: Literal["all"] | Iterable[str] = "all",
    use_raw: bool | None = None,
    reference: str = "rest",
    n_genes: int | None = None,
    tie_correct: bool = False,
    layer: str | None = None,
    chunk_size: int | None = None,
    corr_method: str = "benjamini-hochberg",
    **kwds,
) -> None:
    """GPU-accelerated Wilcoxon rank-sum test with Scanpy parity."""
    try:
        import cupyx.scipy.special
    except ImportError as exc:
        raise ImportError(
            "CuPy is required for GPU-accelerated Wilcoxon test. "
            "Please install CuPy or use scanpy.tl.rank_genes_groups with method='wilcoxon'"
        ) from exc

    avail_corr = {"benjamini-hochberg", "bonferroni"}
    if corr_method not in avail_corr:
        raise ValueError(f"Correction method must be one of {avail_corr}.")

    (
        groupby,
        normalized_groups,
        reference,
        layer,
        use_raw,
        n_top,
    ) = _norm_params(
        groupby=groupby,
        groups=groups,
        reference=reference,
        layer=layer,
        use_raw=use_raw,
        n_top=n_genes,
        method_kwargs=kwds,
    )

    X, labels, var_names = _get_X_obs_var(
        adata,
        groupby=groupby,
        layer=layer,
        use_raw=use_raw,
    )
    groups_order, group_masks, group_sizes = _group_masks(labels, normalized_groups)

    n_cells, n_genes_total = X.shape
    n_top = n_genes_total if n_top is None else min(n_top, n_genes_total)

    def _resolve_chunk(value: int | None) -> int:
        if value is not None:
            return value
        try:
            free_mem, _ = cp.cuda.runtime.memGetInfo()
            bytes_per_gene = n_cells * 8 * 4
            max_genes = int(0.6 * free_mem / bytes_per_gene)
            return max(min(max_genes, 1000), 100)
        except Exception:
            return 500

    block = _resolve_chunk(chunk_size)

    group_matrix = cp.asarray(group_masks.T, dtype=cp.float64)
    group_sizes_dev = cp.asarray(group_sizes, dtype=cp.float64)

    base = adata.uns.get("log1p", {}).get("base")
    if base is not None:
        expm1_func = lambda x: cp.expm1(x * cp.log(base))
    else:
        expm1_func = cp.expm1

    groups_order_str = [str(g) for g in groups_order]
    scores_store = {group: [] for group in groups_order_str}
    logfc_store = {group: [] for group in groups_order_str}
    pvals_store = {group: [] for group in groups_order_str}
    gene_indices_store = {group: [] for group in groups_order_str}

    for ranks_chunk, X_chunk, gene_slice in _rank_chunks(X, block=block):
        z_chunk, p_chunk = _wilcoxon(
            ranks_chunk,
            group_matrix,
            group_sizes_dev,
            tie_correct=tie_correct,
            n_cells=n_cells,
        )
        logfc_chunk = _effects(
            X_chunk,
            group_matrix,
            group_sizes_dev,
            expm1_func=expm1_func,
        )

        z_host = z_chunk.get()
        p_host = p_chunk.get()
        logfc_host = logfc_chunk.get()
        gene_ids = np.arange(gene_slice.start, gene_slice.stop, dtype=int)
        for group_idx, group in enumerate(groups_order_str):
            scores_store[group].append(z_host[group_idx])
            logfc_store[group].append(logfc_host[group_idx])
            pvals_store[group].append(p_host[group_idx])
            gene_indices_store[group].append(gene_ids)

    params = {
        "groupby": groupby,
        "method": "wilcoxon",
        "reference": reference,
        "use_raw": use_raw,
        "tie_correct": tie_correct,
        "layer": layer,
        "corr_method": corr_method,
    }
    adata.uns["rank_genes_groups"] = _pack_result(
        var_names,
        groups_order_str,
        scores=scores_store,
        logfc=logfc_store,
        pvals=pvals_store,
        gene_indices=gene_indices_store,
        n_top=n_top,
        method=corr_method,
        total_genes=n_genes_total,
        params=params,
    )
