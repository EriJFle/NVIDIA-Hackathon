def nmf_gene_programs_gpu(
    adata,
    n_programs=50,
    layer='X',
    max_iter=200,
    random_state=0
):
    """
    Version with imports inside function to avoid scope issues
    """
    # All imports inside the function to ensure they're in scope
    import cupy as cp
    import cupyx.scipy.sparse
    import scipy.sparse as sp
    import numpy as np
    
    cp.random.seed(random_state)
    
    # Get data - handle both sparse and dense
    if layer == 'X':
        X = adata.X
    else:
        X = adata.layers[layer]
    
    # Proper sparse/dense handling
    if sp.issparse(X):
        # Convert scipy sparse to cupy sparse
        X_gpu = cupyx.scipy.sparse.csr_matrix(X.astype(np.float32))
        is_sparse = True
    else:
        # Dense array
        X_gpu = cp.asarray(X.astype(np.float32))
        is_sparse = False
    
    n_cells, n_genes = X_gpu.shape
    print(f"Running NMF on {n_cells} cells × {n_genes} genes")
    
    # Initialize W and H
    W = cp.random.rand(n_cells, n_programs).astype(cp.float32) * 0.01
    H = cp.random.rand(n_programs, n_genes).astype(cp.float32) * 0.01
    
    # Multiplicative updates
    for iteration in range(max_iter):
        # Update H
        numerator = W.T @ X_gpu
        denominator = W.T @ W @ H + 1e-10
        H *= numerator / denominator
        
        # Update W
        if is_sparse:
            numerator = X_gpu @ H.T
            # Handle sparse output
            if hasattr(numerator, 'todense'):
                numerator = cp.asarray(numerator.todense())
        else:
            numerator = X_gpu @ H.T
        denominator = W @ H @ H.T + 1e-10
        W *= numerator / denominator
        
        if iteration % 50 == 0:
            print(f"Iteration {iteration}/{max_iter}")
    
    # Move back to CPU
    adata.obsm['NMF_cell_programs'] = cp.asnumpy(W)
    adata.varm['NMF_gene_loadings'] = cp.asnumpy(H.T)
    
    # Store top genes per program
    for i in range(min(10, n_programs)):
        top_gene_idx = cp.argsort(H[i])[-30:]
        adata.uns[f'program_{i}_top_genes'] = adata.var_names[cp.asnumpy(top_gene_idx)].tolist()
    
    print(f"NMF complete. Results stored in:")
    print(f"  - adata.obsm['NMF_cell_programs']: {adata.obsm['NMF_cell_programs'].shape}")
    print(f"  - adata.varm['NMF_gene_loadings']: {adata.varm['NMF_gene_loadings'].shape}")
    
    return adata