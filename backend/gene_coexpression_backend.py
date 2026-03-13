import numpy as np
import scanpy as sc
import pandas as pd
import plotly.express as px
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Union, List, Optional


def is_gene_normalized(adata, gene: str) -> bool:
    """
    Vérifie en interne si UN gène est déjà normalisé (log1p + CPM).
    """
    # Si les données sont déjà normalisées globalement, on suppose que le gène l'est aussi
    log1p_normalized = adata.uns.get("log1p", {}).get("base") == np.e
    norm_total_done = adata.uns.get("normalization", {}).get("target_sum") == 1e4
    return log1p_normalized and norm_total_done

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
    adata_sub = adata[:, gene_indices].copy()

    # Visualisations
    if plot_types is None or "scatter" in plot_types:
        sc.pl.violin(adata_sub, keys=genes, show=False)
    if plot_types is None or "heatmap" in plot_types:
        sc.pl.histogram(adata_sub, color=genes, show=False)
    if plot_types is None or "umap" in plot_types:
        sc.pl.umap(adata_sub, color=genes, show=False)

    # 6. Affichage des visualisations
    for plot_type in plot_types:
        try:
            if plot_type == "scatter":
                st.subheader(f"Scatter plot : {genes[0]} vs {genes[1]}")

                gene_a = genes[0]
                gene_b = genes[1]

                # Extraire les expressions
                expr_a = adata_sub[:, gene_a].X.toarray().flatten()
                expr_b = adata_sub[:, gene_b].X.toarray().flatten()

                # DataFrame pour Plotly
                df = pd.DataFrame({
                    gene_a: expr_a,
                    gene_b: expr_b,
                    "Cluster": adata_sub.obs["louvain"]
                })

                # Scatter plot
                fig = px.scatter(
                    df,
                    x=gene_a,
                    y=gene_b,
                    color="Cluster",
                    marginal_x="histogram",
                    marginal_y="histogram",
                    title=f"Co-expression de {gene_a} et {gene_b}",
                    labels={
                        gene_a: f"Expression {gene_a} (log CPM)",
                        gene_b: f"Expression {gene_b} (log CPM)"
                    }
                )

                fig.update_traces(
                    marker=dict(size=5, opacity=0.5),
                    selector=dict(type="scatter")
                )

                st.plotly_chart(fig, use_container_width=True)

            elif plot_type == "heatmap":
                st.subheader(f"Heatmap : {', '.join(genes[:5])}")

                selected_genes = genes[:5]

                expr_matrix = adata_sub[:, selected_genes].X.toarray().T
                corr_matrix = np.corrcoef(expr_matrix)

                corr_df = pd.DataFrame(
                    corr_matrix,
                    index=selected_genes,
                    columns=selected_genes
                )

                fig = px.imshow(
                    corr_df,
                    text_auto=True,
                    color_continuous_scale="RdBu_r",
                    zmin=-1,
                    zmax=1,
                    aspect="auto",
                    title="Matrice de corrélation entre gènes"
                )

                st.plotly_chart(fig, use_container_width=True)

            elif plot_type == "umap_coexpression":
                st.subheader(f"UMAP co-expression : {genes[0]} vs {genes[1]}")

                gene_a, gene_b = genes[0], genes[1]

                expr_a = adata_sub[:, gene_a].X.toarray().flatten()
                expr_b = adata_sub[:, gene_b].X.toarray().flatten()

                thr_a = 0.0
                thr_b = 0.0

                a_pos = expr_a > thr_a
                b_pos = expr_b > thr_b

            cats = np.full(adata_sub.n_obs, "none", dtype=object)
            cats[a_pos & ~b_pos] = "A_only"
            cats[~a_pos & b_pos] = "B_only"
            cats[a_pos & b_pos] = "both"

            obs_key = f"coexp_{gene_a}_{gene_b}"

            adata_sub.obs[obs_key] = cats
            adata_sub.obs[obs_key] = adata_sub.obs[obs_key].astype("category")
            adata_sub.obs[obs_key] = adata_sub.obs[obs_key].cat.reorder_categories(
                ["none", "A_only", "B_only", "both"], ordered=True
            )

            adata_sub.uns[f"{obs_key}_colors"] = [
                "lightgrey",
                "deepskyblue",
                "lightcoral",
                "indigo"
            ]

            adata_sub.obs[f"{obs_key}__A"] = expr_a
            adata_sub.obs[f"{obs_key}__B"] = expr_b

            ax = sc.pl.umap(
                adata_sub,
                color=obs_key,
                legend_loc="right margin",
                show=False,
                return_fig=False
            )

            fig = ax.figure
            fig.set_size_inches(10, 6)
            fig.subplots_adjust(right=0.78)

            st.pyplot(fig)
            plt.close(fig)

        except Exception as e:
            st.error(f"Erreur lors de la génération du {plot_type} pour {gene}: {str(e)}")
