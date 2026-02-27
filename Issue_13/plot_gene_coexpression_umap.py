import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt

def plot_gene_coexpression_umap(adata, gene_a: str, gene_b: str):

    # Paramètres
    n_neighbors = 15
    n_pcs = 50

    # Extraire l'expression des gènes A et B (log-transformée)
    expr_a = adata[:, gene_a].X.toarray().flatten()
    expr_b = adata[:, gene_b].X.toarray().flatten()

    # Seuil
    thr_a = 0.0
    thr_b = 0.0

    # Détection co-expression
    a_pos = expr_a > thr_a
    b_pos = expr_b > thr_b

    cats = np.full(adata.n_obs, "none", dtype=object)
    cats[a_pos & ~b_pos] = "A_only"
    cats[~a_pos & b_pos] = "B_only"
    cats[a_pos & b_pos] = "both"

    # Stockage dans obs
    obs_key = f"coexp_{gene_a}_{gene_b}"
    adata.obs[obs_key] = cats
    adata.obs[obs_key] = adata.obs[obs_key].astype("category")
    adata.obs[obs_key] = adata.obs[obs_key].cat.reorder_categories(
        ["none", "A_only", "B_only", "both"], ordered=True
    )

    adata.uns[f"{obs_key}_colors"] = [
        "lightgrey",
        "deepskyblue",
        "lightcoral",
        "indigo"
    ]


    adata.obs[f"{obs_key}__A"] = expr_a
    adata.obs[f"{obs_key}__B"] = expr_b
    
    # UMAP colorée par catégories
    ax = sc.pl.umap(
        adata,
        color=obs_key,
        legend_loc="right margin",
        show=False,
        return_fig=False
    )

    fig = ax.figure
    fig.set_size_inches(10,6)
    fig.subplots_adjust(right=0.78)

    plt.show()
