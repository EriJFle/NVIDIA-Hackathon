def nmf_gene_programs_gpu(
    adata,
    n_programs=50,
    layer='X',
    max_iter=200,
    random_state=0
):
    """
    NMF for data already on GPU
    """
    import cupy as cp
    import cupyx.scipy.sparse
    import numpy as np
    
    cp.random.seed(random_state)
    
    # Get data
    if layer == 'X':
        X_gpu = adata.X  # Already a cupyx sparse matrix!
    else:
        X_gpu = adata.layers[layer]
    
    # Check if it's already on GPU
    if isinstance(X_gpu, cupyx.scipy.sparse.csr_matrix):
        print("Data already on GPU as sparse matrix")
        is_sparse = True
    elif isinstance(X_gpu, cp.ndarray):
        print("Data already on GPU as dense array")
        is_sparse = False
    else:
        print(f"Unexpected type: {type(X_gpu)}")
        return None
    
    n_cells, n_genes = X_gpu.shape
    print(f"Running NMF on {n_cells} cells × {n_genes} genes")
    
    # Initialize W and H - smaller initialization for stability
    W = cp.random.rand(n_cells, n_programs).astype(cp.float32) * 0.01
    H = cp.random.rand(n_programs, n_genes).astype(cp.float32) * 0.01
    
    # Multiplicative updates
    for iteration in range(max_iter):
        # Update H
        numerator = W.T @ X_gpu  # This works with cupyx sparse
        denominator = W.T @ W @ H + 1e-10
        H *= numerator / denominator
        
        # Update W
        numerator = X_gpu @ H.T  # Returns dense when sparse @ dense
        denominator = W @ H @ H.T + 1e-10
        W *= numerator / denominator
        
        # Clear cache periodically to prevent memory buildup
        if iteration % 10 == 0:
            cp.get_default_memory_pool().free_all_blocks()
            
        if iteration % 50 == 0:
            print(f"Iteration {iteration}/{max_iter}")
    
    # Move results back to CPU
    W_cpu = cp.asnumpy(W)
    H_cpu = cp.asnumpy(H)
    
    # Store in adata
    adata.obsm['NMF_cell_programs'] = W_cpu
    adata.varm['NMF_gene_loadings'] = H_cpu.T
    
    # Get top genes per program
    for i in range(min(10, n_programs)):
        top_gene_idx = np.argsort(H_cpu[i])[-30:]
        adata.uns[f'program_{i}_top_genes'] = adata.var_names[top_gene_idx].tolist()
    
    print(f"NMF complete!")
    print(f"  - Cell programs shape: {W_cpu.shape}")
    print(f"  - Gene loadings shape: {H_cpu.T.shape}")
    
    return adata