# AstroWRLD

A comprehensive single-cell transcriptomic reference atlas of mouse astrocytes for disease modeling and seamless, reproducible dataset integration.

---

**Note**: This is an active research project. The codebase and methods are under development and subject to change.

## Overview

This project aims to build a robust reference atlas of mouse astrocytes in health and disease. Using ~300,000 astrocytes from the Allen Brain Cell Atlas (scRNA-seq) and ~480,000 astrocytes from BrainCellData (snRNA-seq) as a reference atlas, this project will enable researchers to perform reference-guided analysis of query datasets through probabilistic cell type prediction, batch-corrected integration, and disease state characterization. Beyond annotation transfer, the atlas will support novel astrocyte subtype discovery, disease progression modeling, and cross-modal data integration.

## Objectives

1. **Build a comprehensive reference** from Allen Brain Cell Atlas astrocytes across all brain regions
2. **Enable disease dataset integration** using scVI-tools and scArches for query-to-reference mapping
3. **Characterize disease-associated astrocyte states** across multiple neurological conditions (EAE/MS, Alzheimer's disease, glioma, TBI, epilepsy)
4. **Systematically benchmark integration methods** using scib-metrics to identify optimal algorithms (scVI, scANVI, scArches variants) for astrocyte-specific integration challenges
5. **Link transcriptomic signatures to functional phenotypes** by integrating gene expression with astrocyte functional properties and disease-relevant behaviors


## Data Sources

### Reference Datasets
- **Allen Brain Cell Atlas**: ~300,000 mouse astrocytes
- 10x Genomics scRNA-seq
- Comprehensive brain coverage

- **BrainCellData**: ~480,000 mouse astrocytes
- 10x Genomics snRNA-seq (Macosko lab)
- Complementary to Allen Atlas
- Provides single-nucleus perspective for cross-modality validation

### Disease Datasets (Public)
- Multiple scRNA-seq and snRNA-seq datasets from GEO/SRA
- Disease models: EAE, AD (5xFAD, APP/PS1), glioma, TBI, epilepsy
- *Specific datasets will be documented as integrated*

### Functional Annotations (Literature-derived)
- Morphological measurements from published imaging studies
- Calcium signaling profiles from GCaMP recordings
- Metabolic states from seahorse/proteomics studies
- Regional functional specializations from circuit mapping

## Methods

### Core Technologies
- **Data Processing**: Scanpy
- **Preservation of Astrocyte Biology and Batch-Aware Integration**: scVI-tools (scVI and scANVI)
- **Reference Mapping**: scArches
- **Benchmarking**: scib-metrics

### Integration Benchmarking Strategy

A critical component of this project is systematic benchmarking to identify the optimal integration approach for astrocyte data:

- **Comprehensive metrics evaluation** using scib-metrics:
  - Biological conservation (ARI, NMI, cell type ASW, isolated label scores)
  - Batch correction (kBET, graph connectivity, iLISI)
  - Overall integration quality (aggregate scores)
- **Method comparison matrix**:
  - scVI with varying latent dimensions
  - scANVI with different semi-supervised strategies
  - scArches with multiple surgery configurations
  - Harmony, BBKNN, and Scanorama as baseline comparisons
- **Astrocyte-specific optimization**:
  - Evaluating performance on subtypes (e.g., hippocampal vs cortical)
  - Testing robustness to technical variation (sc vs sn)
  - Assessing disease state mapping accuracy

### Pipeline Architecture

## Installation

This project uses `uv` for dependency management and includes access to the Allen Brain Cell Atlas via `abc-atlas-access`.

    # Install uv if needed
    curl -LsSf https://astral.sh/uv/install.sh | sh
    # or on macOS: brew install astral-sh/uv/uv

    # Clone and set up the core env
    git clone https://github.com/EriJFle/AstroWRLD.git
    cd AstroWRLD
    uv venv
    source .venv/bin/activate

    uv pip install \
        --extra-index-url=https://pypi.nvidia.com \
        "cudf-cu12==25.8.*" "dask-cudf-cu12==25.8.*" "cuml-cu12==25.8.*" \
        "cugraph-cu12==25.8.*" "nx-cugraph-cu12==25.8.*" "cuxfilter-cu12==25.8.*" \
        "cucim-cu12==25.8.*" "pylibraft-cu12==25.8.*" "raft-dask-cu12==25.8.*" \
        "cuvs-cu12==25.8.*" "nx-cugraph-cu12==25.8.*"

    uv pip install scib-metrics

    uv pip install scanpy scvi-tools

## Usage

*Coming soon - will include example code for:*
- Loading the reference atlas
- Mapping your own astrocyte dataset
- Running standard analyses

## Contributing

This project will be a large undertaking and will require a team of passionate researchers working towards a common goal. We welcome contributions! Areas where help is particularly needed:

1. **Functional annotation curation**: Literature mining for gene-function relationships
2. **Disease dataset curation**: Identifying and preprocessing public datasets
3. **Prediction model development**: Improving function prediction algorithms
4. **Benchmarking extensions**: Testing additional integration methods
5. **Validation**: Testing predictions with new functional datasets
6. **Documentation**: Improving tutorials and examples

Please open an issue to discuss potential contributions or submit a pull request.

## Key Design Decisions

- **Astrocyte-specific QC**: Accounting for unique astrocyte characteristics (e.g., high Aqp4, variable GFAP)
- **sc/sn compatibility**: Using layer normalization in models for cross-modality integration
- **Disease focus**: Feature selection includes established disease-associated gene programs
- **Hierarchical annotation**: Leveraging Allen Brain Atlas annotation depth for flexible mapping
- **Rigorous benchmarking**: Using scib-metrics to ensure optimal integration method selection rather than arbitrary choice
- **Function-first approach**: Prioritizing biological interpretability through systematic gene-function mapping

## Current Status

**Under Development**

- [x] Project conception and design
- [x] Reference atlas construction
- [ ] Functional annotation framework
- [ ] Integration method benchmarking
- [ ] Disease dataset collection
- [ ] Function prediction model training
- [ ] Validation with independent datasets
- [ ] Documentation and tutorials
- [ ] Public release

## Citation

For now, please cite the Allen Brain Cell Atlas:
Yao, Z., van Velthoven, C.T.J., Kunst, M. et al. A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain. Nature 624, 317–332 (2023). https://doi.org/10.1038/s41586-023-06812-z

and BrainCellData:
Langlieb, J., Sachdev, N. S., Balderrama, K. S., Nadaf, N. M., Raj, M., Murray, E., Webber, J. T., Vanderburg, C., Gazestani, V., Tward, D., Mezias, C., Li, X., Flowers, K., Cable, D. M., Norton, T., Mitra, P., Chen, F., & Macosko, E. Z. (2023). The molecular cytoarchitecture of the adult mouse brain. Nature, 624(7991), 333–342. https://doi.org/10.1038/s41586-023-06818-7

## Contact

*Erik Fleischel, evz9yd@virginia.edu*
