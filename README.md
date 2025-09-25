# NVIDIA Hackathon

GPU-accelerated implementation of Scanpy's rank_genes_groups Wilcoxon rank-sum test. Developed for the NVIDIA x scverse Hackathon and for integration with rapids_singlecell ecosystem. The repo holds notebooks, scripts, and source modules for benchmarking and prototyping analyses.

## Results
- 5-15x speedup of wilcoxon rank-sum test
- Perfect parity with Scanpy's rank_genes_groups(method = "wilcoxon") output
- Up to 2x faster than rapids_singlecell.tl.rank_genes_groups_logreg() for large (500k+ cell) datasets
- Automatic memory-based chunking strategy enables practically any GPU to be used. The GPU-accelerated Wilcoxon rank-sum was tested and validated on an RTX2080Ti and A100.

## Dependencies
- All notebooks were run with NVIDIA's preconfigued NVIDIA Brev resources (https://github.com/NVIDIA-AI-Blueprints/single-cell-analysis-blueprint)
- Install the required packages from `requirements.txt` to run AstroWRLD tests and analysis.


### Installation
```bash
uv pip install -r requirements.txt
```
