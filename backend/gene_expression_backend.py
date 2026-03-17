import numpy as np
import scanpy as sc
import scipy.sparse as sp
from typing import Union, List, Optional

def plot_gene_expression(
    adata,
    gene,
    plot_types: Optional[List[str]] = None,
) -> None:

    for plot_type in plot_types:
        try:
            if plot_type == "umap_clusters":
                st.subheader("UMAP par clusters (Louvain)")
                fig, ax = plt.subplots()
                sc.pl.umap(
                    adata,
                    color="louvain",
                    palette="viridis",
                    legend_loc="on data",
                    legend_fontsize=8,
                    size=50,
                    show=False,
                    ax=ax
                )
                st.pyplot(fig)
                plt.close(fig)

            elif plot_type == "violin":
                st.subheader(f"Violin plot : {gene}")
                fig, ax = plt.subplots()
                sc.pl.violin(
                    adata,
                    keys=[gene],
                    show=False,
                    ax=ax,
                    groupby="louvain",
                    stripplot=False
                )
                ax.set_xlabel("Cluster Louvain")
                ax.set_ylabel(f"Expression normalisée de {gene}")
                st.pyplot(fig)
                plt.close(fig)

            elif plot_type == "histogram":
                st.subheader(f"Histogramme : {gene}")

                expr = adata[:, gene].X
                if hasattr(expr, "toarray"):
                    expr = expr.toarray()
                expr = np.asarray(expr).reshape(-1)

                fig, ax = plt.subplots()
                ax.hist(expr, bins=50)
                ax.set_title(f"Distribution de {gene}")
                ax.set_xlabel("Expression")
                ax.set_ylabel("Nombre de cellules")
                st.pyplot(fig)
                plt.close(fig)

            elif plot_type == "umap":
                st.subheader(f"UMAP : {gene}")
                fig, ax = plt.subplots()
                sc.pl.umap(
                    adata,
                    color=gene,
                    color_map="viridis",
                    size=50,
                    show=False,
                    ax=ax
                )
                ax.set_title(f"Expression de {gene}")
                st.pyplot(fig)
                plt.close(fig)

        except Exception as e:
            st.error(f"Erreur lors de la génération du {plot_type} pour {gene}: {str(e)}")
