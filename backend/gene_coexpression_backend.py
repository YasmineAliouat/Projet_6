import numpy as np
import scanpy as sc
import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt
import streamlit as st
from typing import List, Optional
from scipy.stats import pearsonr

def plot_gene_coexpression(
    adata,
    genes: List[str],
    plot_types: Optional[List[str]] = None,
) -> None:

    gene_a, gene_b = genes[0], genes[1]

    for plot_type in plot_types:
        try:
            if plot_type == "scatter":

                st.subheader(f"Scatter plot : {gene_a} vs {gene_b}")

                expr_a = adata[:, gene_a].X
                expr_b = adata[:, gene_b].X

                if hasattr(expr_a, "toarray"):
                    expr_a = expr_a.toarray()
                if hasattr(expr_b, "toarray"):
                    expr_b = expr_b.toarray()

                expr_a = np.asarray(expr_a).flatten()
                expr_b = np.asarray(expr_b).flatten()

                clusters = adata.obs["louvain"].astype(str).values

                df = pd.DataFrame({
                    gene_a: expr_a,
                    gene_b: expr_b,
                    "Cluster": clusters
                })

                fig = px.scatter(
                    df,
                    x=gene_a,
                    y=gene_b,
                    color="Cluster",
                    marginal_x="histogram",
                    marginal_y="histogram",
                    opacity=0.6
                )

                # Calculer les coefficients de corrélation par cluster
                corr_text = ""
                for cluster in df["Cluster"].unique():
                    cluster_data = df[df["Cluster"] == cluster]
                    if len(cluster_data) > 2:  # Pearson nécessite au moins 3 points
                        r, _ = pearsonr(cluster_data[gene_a], cluster_data[gene_b])
                        corr_text += f"Cluster {cluster}: R = {r:.2f}<br>"
                

                fig.add_annotation(
                    xref="paper", yref="paper",
                    x=0.95, y=0.95,
                    text=corr_text,
                    showarrow=False,
                    bordercolor="black",
                    borderwidth=1,
                    borderpad=4,
                    bgcolor="white",
                    font=dict(size=12, color="black")
                )

                st.plotly_chart(fig, use_container_width=True)

            elif plot_type == "heatmap":

                st.subheader(f"Heatmap corrélation : {gene_a}, {gene_b}")

                expr_matrix = adata[:, [gene_a, gene_b]].X
                if hasattr(expr_matrix, "toarray"):
                    expr_matrix = expr_matrix.toarray()

                corr_matrix = np.corrcoef(expr_matrix.T)

                corr_df = pd.DataFrame(
                    corr_matrix,
                    index=[gene_a, gene_b],
                    columns=[gene_a, gene_b]
                )

                fig = px.imshow(
                    corr_df,
                    text_auto=True,
                    color_continuous_scale="RdBu_r",
                    zmin=-1,
                    zmax=1,
                    aspect="auto",
                )

                st.plotly_chart(fig, use_container_width=True)

            elif plot_type == "umap":

                st.subheader(f"UMAP co-expression : {gene_a} / {gene_b}")

                expr_a = adata[:, gene_a].X
                expr_b = adata[:, gene_b].X

                if hasattr(expr_a, "toarray"):
                    expr_a = expr_a.toarray()
                if hasattr(expr_b, "toarray"):
                    expr_b = expr_b.toarray()

                expr_a = np.asarray(expr_a).flatten()
                expr_b = np.asarray(expr_b).flatten()

                thr_a = 0.0
                thr_b = 0.0

                a_pos = expr_a > thr_a
                b_pos = expr_b > thr_b

                categories = np.full(adata.n_obs, "none", dtype=object)
                categories[a_pos & ~b_pos] = f"{gene_a} only"
                categories[~a_pos & b_pos] = f"{gene_b} only"
                categories[a_pos & b_pos] = "both"

                obs_key = f"coexp_{gene_a}_{gene_b}"

                adata.obs[obs_key] = pd.Categorical(
                    categories,
                    categories=["none", f"{gene_a} only", f"{gene_b} only", "both"],
                    ordered=True
                )

                adata.uns[f"{obs_key}_colors"] = [
                    "lightgrey",
                    "deepskyblue",
                    "lightcoral",
                    "indigo"
                ]

                fig, ax = plt.subplots()

                sc.pl.umap(
                    adata,
                    color=obs_key,
                    legend_loc="right margin",
                    show=False,
                    ax=ax
                )

                st.pyplot(fig)
                plt.close(fig)

        except Exception as e:
            st.error(f"Erreur lors de la génération du {plot_type} : {str(e)}")
