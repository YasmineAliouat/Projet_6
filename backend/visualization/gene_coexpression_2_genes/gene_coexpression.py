import io
import numpy as np
import scanpy as sc
import pandas as pd
import plotly.express as px
import streamlit as st
from typing import List, Optional
from scipy.stats import pearsonr


def _show(fig, ratio=None):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=120)
    buf.seek(0)
    _, col, _ = st.columns(ratio or [1, 3, 1])
    with col:
        st.image(buf, use_container_width=True)


def _display_name(adata, var_name):
    if "gene_symbol" in adata.var.columns:
        try:
            sym = adata.var.at[var_name, "gene_symbol"]
            if pd.notna(sym) and str(sym).strip():
                return str(sym)
        except KeyError:
            pass
    return var_name


def plot_gene_coexpression(
    adata,
    genes: List[str],
    plot_types: Optional[List[str]] = None,
    display_names: Optional[List[str]] = None,
):

    figs = []

    gene_a, gene_b = genes[0], genes[1]
    display_a = (display_names[0] if display_names and display_names[0] else None) or _display_name(adata, gene_a)
    display_b = (display_names[1] if display_names and len(display_names) > 1 and display_names[1] else None) or _display_name(adata, gene_b)

    for plot_type in plot_types:
        try:
            if plot_type == "scatter":

                st.subheader(f"Scatter plot of {display_a} vs {display_b}")

                expr_a = adata[:, gene_a].X
                expr_b = adata[:, gene_b].X

                if hasattr(expr_a, "toarray"):
                    expr_a = expr_a.toarray()
                if hasattr(expr_b, "toarray"):
                    expr_b = expr_b.toarray()

                expr_a = np.asarray(expr_a).flatten()
                expr_b = np.asarray(expr_b).flatten()

                cell_types = adata.obs["cell_type"].astype(str).values

                df = pd.DataFrame({
                    display_a: expr_a,
                    display_b: expr_b,
                    "Cell type": cell_types
                })

                fig = px.scatter(
                    df,
                    x=display_a,
                    y=display_b,
                    color="Cell type",
                    marginal_x="histogram",
                    marginal_y="histogram",
                    opacity=0.85,
                    color_discrete_sequence=px.colors.qualitative.Set1,
                    template="plotly_white"
                )
                fig.update_traces(
                    marker=dict(size=6),
                    selector=dict(mode="markers")
                )
                fig.update_layout(
                    paper_bgcolor="white",
                    plot_bgcolor="white",
                    font=dict(color="black", size=13),
                    legend=dict(
                        font=dict(color="black", size=12),
                        title=dict(font=dict(color="black", size=13)),
                    ),
                    margin=dict(r=180, t=60),
                )
                fig.update_xaxes(
                    gridcolor="#e5e5e5",
                    linecolor="#cccccc",
                    zeroline=False,
                    tickfont=dict(color="black", size=12),
                    title_font=dict(color="black", size=13),
                )
                fig.update_yaxes(
                    gridcolor="#e5e5e5",
                    linecolor="#cccccc",
                    zeroline=False,
                    tickfont=dict(color="black", size=12),
                    title_font=dict(color="black", size=13),
                )

                corr_text = ""
                for cell_type in df["Cell type"].unique():
                    cluster_data = df[df["Cell type"] == cell_type]
                    if len(cluster_data) > 2:
                        r, _ = pearsonr(cluster_data[display_a], cluster_data[display_b])
                        corr_text += f"{cell_type}: R = {r:.2f}<br>"

                fig.add_annotation(
                    xref="paper", yref="paper",
                    x=0.95, y=0.95,
                    text=corr_text,
                    showarrow=False,
                    bordercolor="gray",
                    borderwidth=1,
                    borderpad=4,
                    bgcolor="rgba(128,128,128,0.15)",
                    font=dict(size=12)
                )

                st.plotly_chart(fig, use_container_width=True)
                figs.append(fig)

            elif plot_type == "umap":

                st.subheader(f"UMAP co-expression of {display_a} vs {display_b}")

                expr_a = adata[:, gene_a].X
                expr_b = adata[:, gene_b].X

                if hasattr(expr_a, "toarray"):
                    expr_a = expr_a.toarray()
                if hasattr(expr_b, "toarray"):
                    expr_b = expr_b.toarray()

                expr_a = np.asarray(expr_a).flatten()
                expr_b = np.asarray(expr_b).flatten()

                a_pos = expr_a > 0
                b_pos = expr_b > 0

                categories = np.full(adata.n_obs, "none", dtype=object)
                categories[a_pos & ~b_pos] = f"{display_a} only"
                categories[~a_pos & b_pos] = f"{display_b} only"
                categories[a_pos & b_pos] = "both"

                obs_key = f"coexp_{gene_a}_{gene_b}"

                adata.obs[obs_key] = pd.Categorical(categories)

                custom_palette = {
                    "none": "darkgray",
                    f"{display_a} only": "mediumseagreen",
                    f"{display_b} only": "darkorange",
                    "both": "purple"
                }

                fig, ax = plt.subplots(figsize=(4, 3))

                sc.pl.umap(
                    adata,
                    color=obs_key,
                    palette=custom_palette,
                    show=False,
                    ax=ax
                )

                ax.set_title(f"coexpression of {display_a} and {display_b}")
                legend = ax.get_legend()
                if legend is not None:
                    legend.set_title(f"{display_a} / {display_b}")

                _show(fig)
                figs.append(fig)

        except Exception as e:
            st.error(f"Error during {plot_type}: {str(e)}")

    return figs