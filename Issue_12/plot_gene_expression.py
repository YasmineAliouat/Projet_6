import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt

def plot_gene_expression(adata, gene: str, normalize: bool = False):
    """
    Affiche:
      - violin plot de l'expression
      - histogramme (distribution)
      - UMAP colorée par expression du gène

    """

    if normalize:
        # Normalisation en CPM + log1p si ce n'est pas déjà fait
        if not adata.uns.get("log1p", {}).get("base") == np.e:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
        else:
            print("Les données sont déjà normalisées (log1p détecté)")

    # Violin plot
    sc.pl.violin(
        adata,
        keys=[gene],
        show=True
    )

    # Histogramme simple
    expr = adata[:, gene].X
    if hasattr(expr, "toarray"):
        expr = expr.toarray()
    expr = np.asarray(expr).reshape(-1)

    plt.figure()
    plt.hist(expr, bins=50)
    plt.title(f"Distribution expression : {gene}")
    plt.xlabel("Expression")
    plt.ylabel("Nombre de cellules")
    plt.show()

    # UMAP colorée par expression
    sc.pl.umap(
        adata,  
        color=gene,
        color_map="viridis",
        size=50,
        title=f"Expression de {gene} (UMAP)"
    )

    plt.show()