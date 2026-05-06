import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from typing import List, Optional


def _display_name(adata, var_name):
    if "gene_symbol" in adata.var.columns:
        try:
            sym = adata.var.at[var_name, "gene_symbol"]
            if pd.notna(sym) and str(sym).strip():
                return str(sym)
        except KeyError:
            pass
    return var_name


def coexp_heatmap(adata, gene_list: List[str], display_names: Optional[List[str]] = None):

    if display_names and len(display_names) == len(gene_list):
        display_list = [d if d else _display_name(adata, g) for d, g in zip(display_names, gene_list)]
    else:
        display_list = [_display_name(adata, g) for g in gene_list]

    # Utiliser une copie minimale en renomant var_name pour que la légende de le heatmap montre les noms matchés
    if display_list != gene_list and len(set(display_list)) == len(display_list):
        try:
            adata_plot = adata[:, gene_list].copy()
            adata_plot.var_names = display_list
        except Exception:
            adata_plot = adata
            display_list = gene_list
    else:
        adata_plot = adata
        display_list = gene_list

    try:
        axes_dict = sc.pl.heatmap(
            adata_plot,
            var_names=display_list,
            groupby="cell_type",
            cmap="RdBu_r",
            standard_scale="var",
            dendrogram=False,
            show_gene_labels=True,
            show=False,
            figsize=(5, 4),
        )

        if isinstance(axes_dict, dict):
            groupby_ax = axes_dict.get("groupby_ax")
            if groupby_ax is not None:
                groupby_ax.set_xticklabels([])
                groupby_ax.set_yticklabels([])

            if "cell_type_colors" in adata.uns:
                categories = adata.obs["cell_type"].cat.categories
                colors = adata.uns["cell_type_colors"]
                patches = [
                    mpatches.Patch(color=colors[i], label=cat)
                    for i, cat in enumerate(categories)
                ]
                plt.gcf().legend(
                    handles=patches,
                    loc="upper left",
                    bbox_to_anchor=(1.02, 1.0),
                    ncol=1,
                    fontsize=6,
                    title="Cell type",
                    title_fontsize=7,
                    frameon=True,
                    borderaxespad=0,
                )

    except Exception as e:
        raise RuntimeError(f"Error generating heatmap: {str(e)}")

    return axes_dict