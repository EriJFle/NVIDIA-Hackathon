# Correlate programs with cell types
if 'cell_type' in adata.obs:
    program_scores = pd.DataFrame(
        adata.obsm['NMF_cell_programs'],
        index=adata.obs_names
    )
    # Mean program usage per cell type
    program_by_celltype = program_scores.groupby(adata.obs['cell_type']).mean()
    # Heatmap
    import seaborn as sns
    plt.figure(figsize=(10, 8))
    sns.heatmap(program_by_celltype.T, cmap='RdBu_r', center=0)
    plt.xlabel('Cell Type')
    plt.ylabel('Program')
    plt.savefig("celltype_program.png")
    plt.show()