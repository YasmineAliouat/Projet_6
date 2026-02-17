import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt

def plot_gene_expression(
    adata,  # fichier qui devra être chargé en amont
    gene: str,
    *,
    layer: str | None = None,
    use_raw: bool = False,
    ensure_umap: bool = True,
    n_neighbors: int = 15,
    n_pcs: int = 50,
):
    """
    Affiche:
      - violin plot de l'expression
      - histogramme (distribution)
      - UMAP colorée par expression du gène

    """

    # Vérifier que le gène existe
    if use_raw:
        if adata.raw is None:
            raise ValueError("use_raw=True mais adata.raw est None (raw absent).")
        gene_space = adata.raw.var_names
    else:
        gene_space = adata.var_names

    if gene not in gene_space:
        # suggestions simples (contient la chaîne)
        suggestions = [g for g in gene_space if gene.upper() in str(g).upper()][:10]
        msg = f"Gène '{gene}' introuvable dans {'adata.raw.var_names' if use_raw else 'adata.var_names'}."
        if suggestions:
            msg += f" Suggestions (contient '{gene}') : {suggestions}"
        raise KeyError(msg)

    # Violin plot
    sc.pl.violin(
        adata,
        keys=[gene],
        use_raw=use_raw,
        layer=layer,
        show=True
    )

    # Histogramme simple
    expr = adata.raw[:, gene].X if use_raw else adata[:, gene].X
    if layer is not None and not use_raw:
        expr = adata[:, gene].layers[layer]

    # expr peut être sparse ou dense -> on convertit en 1D proprement
    if hasattr(expr, "toarray"):
        expr = expr.toarray()
    expr = np.asarray(expr).reshape(-1)

    plt.figure()
    plt.hist(expr, bins=50)
    plt.title(f"Distribution expression : {gene}")
    plt.xlabel("Expression")
    plt.ylabel("Nombre de cellules")
    plt.show()

    # Présence UMAP
    if ensure_umap and "X_umap" not in adata.obsm_keys():
        sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=min(n_pcs, adata.obsm["X_pca"].shape[1]))
        sc.tl.umap(adata)

    # UMAP colorée par expression
    sc.pl.umap(
        adata,  # Scanpy sait que l'UMAP est dans .obsm["X_umap"]
        color=gene,
        use_raw=use_raw,
        layer=layer,
        show=True
    )