"""GPU-accelerated differential expression utilities (experimental).

This package currently provides a GPU Wilcoxon rank-sum implementation
compatible with Scanpy's ``rank_genes_groups`` output structure.

Usage (from repo root without packaging):

  - CLI parity check: ``python -m src.wilcoxon.scripts.parity_check_pbmc3k``
  - Benchmarks: ``python -m src.wilcoxon.benchmarks.bench_rank_genes_groups_wilcoxon``

"""

__all__ = [
    "tl",
]

