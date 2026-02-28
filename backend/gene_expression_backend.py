import numpy as np
import scanpy as sc
import scipy.sparse as sp
from typing import Union, List, Optional

def is_gene_normalized(adata, gene: str) -> bool:
    """
    Vérifie en interne si UN gène est déjà normalisé (log1p + CPM).
    """
    # Si les données sont déjà normalisées globalement, on suppose que le gène l'est aussi
    log1p_normalized = adata.uns.get("log1p", {}).get("base") == np.e
    norm_total_done = adata.uns.get("normalization", {}).get("target_sum") == 1e4
    return log1p_normalized and norm_total_done

def normalize_gene_subset(adata, gene_indices: List[int]) -> sc.AnnData:
    """
    Normalise uniquement un sous-ensemble de gènes (CPM + log1p).
    Retourne un nouvel AnnData avec uniquement ces gènes.
    """
    adata_sub = adata[:, gene_indices].copy()
    sc.pp.normalize_total(adata_sub, target_sum=1e4)
    sc.pp.log1p(adata_sub)
    return adata_sub

def ensure_csr_for_genes(adata, gene_indices: List[int]) -> None:
    """Convertit en CSR uniquement les gènes sélectionnés si sparse."""
    if sp.issparse(adata.X):
        # Extraire uniquement les lignes/colonnes nécessaires
        X_sub = adata.X[:, gene_indices]
        if X_sub.format != "csr":
            adata.X[:, gene_indices] = X_sub.tocsr()

def plot_gene_expression(
    adata,
    plot_types: Optional[List[str]] = None,
) -> None:
    """
    Affiche les visualisations pour un/des gène(s) avec :
    - Normalisation automatique UNIQUEMENT sur les gènes sélectionnés.
    - Optimisation sparse (CSR) UNIQUEMENT pour ces gènes.
    - Aucune sortie utilisateur (silencieux).
    """
    # Vérification des gènes
    genes = [gene] if isinstance(gene, str) else gene
    missing_genes = [g for g in genes if g not in adata.var_names]
    if missing_genes:
        raise ValueError(f"Gènes non trouvés : {', '.join(missing_genes)}")

    # Récupération des indices des gènes
    gene_indices = [adata.var_names.get_loc(g) for g in genes]

    # Vérification/normalisation pour ces gènes
    if is_gene_normalized(adata, genes[0]):  # On suppose que tous les gènes ont le même état
        adata_sub = adata[:, gene_indices].copy()
    else:
        adata_sub = normalize_gene_subset(adata, gene_indices)

    # Optimisation sparse (CSR) pour ces gènes
    ensure_csr_for_genes(adata_sub, range(len(genes)))  

    # Visualisations
    if plot_types is None or "violin" in plot_types:
        sc.pl.violin(adata_sub, keys=genes, show=False)
    if plot_types is None or "histogram" in plot_types:
        sc.pl.histogram(adata_sub, color=genes, show=False)
    if plot_types is None or "umap" in plot_types:
        sc.pl.umap(adata_sub, color=genes, show=False)

    # 6. Affichage des visualisations
    for plot_type in plot_types:
        try:
            if plot_type == "violin":
                st.subheader(f"Violin plot : {', '.join(genes)}")
                fig, ax = plt.subplots()
                sc.pl.violin(adata_sub, keys=genes, show=False, ax=ax)
                st.pyplot(fig)
                plt.close(fig)

            elif plot_type == "histogram":
                st.subheader(f"Histogramme : {', '.join(genes)}")
                for g in genes:
                    expr = adata_sub[:, g].X
                    if hasattr(expr, "toarray"):
                        expr = expr.toarray()
                    expr = np.asarray(expr).reshape(-1)

                    fig, ax = plt.subplots()
                    ax.hist(expr, bins=50)
                    ax.set_title(f"Distribution de {g}")
                    ax.set_xlabel("Expression (log1p)")
                    ax.set_ylabel("Nombre de cellules")
                    st.pyplot(fig)
                    plt.close(fig)

            elif plot_type == "umap":
                st.subheader(f"UMAP : {', '.join(genes)}")
                for g in genes:
                    fig, ax = plt.subplots()
                    sc.pl.umap(adata_sub, color=g, color_map="viridis", size=50, show=False, ax=ax)
                    ax.set_title(f"Expression de {g} (UMAP)")
                    st.pyplot(fig)
                    plt.close(fig)

        except Exception as e:
            st.error(f"Erreur lors de la génération du {plot_type} pour {gene}: {str(e)}")
