import torch
def nmf_pytorch(adata, n_programs=50, device='cuda'):
    """
    NMF using PyTorch - more mature GPU implementation.
    """
    # Handle both sparse and dense matrices
    import scipy.sparse
    if scipy.sparse.issparse(adata.X):
        X_data = adata.X.toarray()
    else:
        X_data = adata.X
    X_torch = torch.tensor(X_data, dtype=torch.float32, device=device)
    # Simple PyTorch NMF implementation
    n_features, n_genes = X_torch.shape
    # Initialize W and H
    W = torch.rand(n_features, n_programs, device=device, requires_grad=True)
    H = torch.rand(n_programs, n_genes, device=device, requires_grad=True)
    # Optimization
    optimizer = torch.optim.Adam([W, H], lr=0.01)
    for epoch in range(1000):
        optimizer.zero_grad()
        # Reconstruction
        X_recon = torch.mm(W, H)
        # Loss (Frobenius norm)
        loss = torch.norm(X_torch - X_recon, 'fro') ** 2
        loss.backward()
        optimizer.step()
        # Non-negativity constraints
        with torch.no_grad():
            W.clamp_(min=0)
            H.clamp_(min=0)
        if epoch % 100 == 0:
            print(f"Epoch {epoch}, Loss: {loss.item():.4f}")
    adata.obsm['NMF_programs'] = W.detach().cpu().numpy()
    adata.varm['NMF_loadings'] = H.T.detach().cpu().numpy()
    return adata
