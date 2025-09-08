import pandas as pd

def interpret_programs(adata, n_top_genes=30):
    """
    Figure out what biological process each program captures
    """
    program_interpretations = {}
    for i in range(adata.varm['NMF_gene_loadings'].shape[1]):
        # Get top genes for this program
        loadings = adata.varm['NMF_gene_loadings'][:, i]
        top_gene_idx = np.argsort(loadings)[-n_top_genes:]
        top_genes = adata.var_names[top_gene_idx].tolist()
        print(f"\nProgram {i} top genes:")
        print(top_genes[:10])  # Print first 10
        # Check for known signatures
        if any('STMN1' in g or 'MKI67' in g for g in top_genes[:20]):
            program_interpretations[i] = "Cell Cycle"
        elif any('MT-' in g for g in top_genes[:10]):
            program_interpretations[i] = "Mitochondrial/Stress"
        elif any('RPS' in g or 'RPL' in g for g in top_genes[:10]):
            program_interpretations[i] = "Ribosomal/Translation"
        # Add more checks based on your biology
    return program_interpretations