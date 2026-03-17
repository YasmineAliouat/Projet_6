import numpy as np
import scanpy as sc
import scipy.sparse as sp
from typing import Union, List, Optional

def is_gene_normalized(adata, gene: str) -> bool:
    """
    Vérifie en interne si UN gène est déjà normalisé (log1p + CPM).
    """
# Récupérer les valeurs du gène (en préservant le sparse)
    gene_data = adata[:, gene].X
    if hasattr(gene_data, "toarray"):
        gene_data = gene_data.toarray().flatten()
    else:
        gene_data = gene_data.flatten()

    # Vérifier si les valeurs sont dans une plage typique après normalisation log1p + CPM
    # - Après log1p, les valeurs sont généralement entre 0 et ~10 (rarement > 20)
    # - Après CPM, la moyenne est souvent autour de 1-5
    is_log1p = np.all(gene_data >= 0) and np.max(gene_data) < 20
    is_cpm = np.mean(gene_data) > 0.1 and np.mean(gene_data) < 10
    return is_log1p and is_cpm

def is_log1p_normalized(gene_data) -> bool:
    """Vérifie si les données sont log1pées."""
    return np.all(gene_data >= 0) and np.max(gene_data) < 20

def is_cpm_normalized(gene_data) -> bool:
    """Vérifie si les données sont CPMées."""
    return np.mean(gene_data) > 0.1 and np.mean(gene_data) < 10

def plot_gene_expression(
    adata,
    gene: Union[str, List[str]],
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
    adata_sub = adata[:, gene_indices].copy()

    # Vérifier et normaliser UNIQUEMENT les gènes non normalisés
    for g in genes:
        if not is_gene_normalized(adata_sub, g):
            # Normaliser uniquement ce gène dans adata_sub
            sc.pp.normalize_total(adata_sub[:, [g]], target_sum=1e4, inplace=True)
            sc.pp.log1p(adata_sub[:, [g]], inplace=True)

    # Visualisations
    if plot_types is None or "umap_cluster" in plot_types:
        sc.pl.umap(adata_sub, color=genes, show=False)
    if plot_types is None or "violin" in plot_types:
        sc.pl.violin(adata_sub, keys=genes, show=False)
    if plot_types is None or "histogram" in plot_types:
        sc.pl.histogram(adata_sub, keys=genes, show=False)
    if plot_types is None or "umap" in plot_types:
        sc.pl.umap(adata_sub, color=genes, show=False)

    # 6. Affichage des visualisations
    for plot_type in plot_types:
        try:
            if plot_type == "umap_cluster":
                st.subheader(f"UMAP par clusters (Louvain)")
                fig, ax = plt.subplots()
                sc.pl.umap(
                adata,
                color="louvain",  
                palette="viridis",  
                legend_loc="on data",  
                legend_fontsize=8,  
                size=50,
                )

            elif plot_type == "violin":
                st.subheader(f"Violin plot : {', '.join(genes)}")
                fig, ax = plt.subplots()
                sc.pl.violin(
                    adata_sub,
                    keys=[gene],
                    show=False,
                    ax=ax,
                    jitter=False,    
                    groupby="louvain",      
                    multi_panel=False, 
                    stripplot=False,       
                    inner=None,    
                    linewidth=0       
                    )
                ax.set_xlabel("Cluster Louvain")
                ax.set_ylabel(f"Expression normalisée de {gene}")
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
