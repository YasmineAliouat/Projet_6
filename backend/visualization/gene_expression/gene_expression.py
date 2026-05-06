import io
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from typing import Union, List, Optional


def _display_name(adata, var_name):
    if "gene_symbol" in adata.var.columns:
        try:
            sym = adata.var.at[var_name, "gene_symbol"]
            if pd.notna(sym) and str(sym).strip():
                return str(sym)
        except KeyError:
            pass
    return var_name


def _show(fig, ratio=None):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=120)
    buf.seek(0)
    _, col, _ = st.columns(ratio or [1, 3, 1])
    with col:
        st.image(buf, use_container_width=True)


def plot_gene_expression(
    adata,
    gene,
    plot_types: Optional[List[str]] = None,
    display_name: Optional[str] = None,
):

    figs = []
    display = display_name if display_name else _display_name(adata, gene)

    for plot_type in plot_types:
        try:
            if plot_type == "umap_clusters":
                st.subheader("UMAP by cell type")
                fig, ax = plt.subplots(figsize=(4, 3))
                sc.pl.umap(
                    adata,
                    color="cell_type",
                    legend_loc="right margin",
                    legend_fontsize=7,
                    size=50,
                    show=False,
                    ax=ax
                )
                _show(fig)
                figs.append(fig)

            elif plot_type == "violin":
                st.subheader(f"Violin plot of {display}")
                fig, ax = plt.subplots(figsize=(6, 3))
                sc.pl.violin(
                    adata,
                    keys=[gene],
                    show=False,
                    ax=ax,
                    groupby="cell_type",
                    stripplot=False
                )
                ax.set_title(display)
                ax.set_xlabel("Cell type", fontweight="bold", fontsize=12)
                ax.set_ylabel(f"Normalized expression of {display} (log1p)")
                _show(fig, ratio=[0.5, 4, 0.5])
                figs.append(fig)

            elif plot_type == "histogram":
                st.subheader(f"Histogram of {display}")

                expr = adata[:, gene].X
                if hasattr(expr, "toarray"):
                    expr = expr.toarray()
                expr = np.asarray(expr).reshape(-1)

                fig, ax = plt.subplots(figsize=(4, 3))
                ax.hist(expr, bins=50)
                ax.set_title(f"Distribution of {display}")
                ax.set_xlabel(f"Normalized expression of {display} (log1p)")
                ax.set_ylabel("Number of cells")
                _show(fig)
                figs.append(fig)

            elif plot_type == "umap":
                st.subheader(f"UMAP of {display}")
                fig, ax = plt.subplots(figsize=(4, 3))
                sc.pl.umap(
                    adata,
                    color=gene,
                    color_map="inferno",
                    size=50,
                    show=False,
                    ax=ax
                )
                ax.set_title(f"Normalized expression of {display} (log1p)")
                if display != gene:
                    for ax_item in fig.get_axes():
                        if ax_item is not ax:
                            if ax_item.get_ylabel() == gene:
                                ax_item.set_ylabel(display)
                            if ax_item.get_title() == gene:
                                ax_item.set_title(display)
                _show(fig)
                figs.append(fig)

        except Exception as e:
            st.error(f"Error generating {plot_type} for {display}: {str(e)}")

    return figs