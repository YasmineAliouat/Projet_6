import scanpy as sc
import matplotlib.pyplot as plt
from typing import List

def coexp_heatmap(adata, gene_list: List[str]):

    try:
        fig = sc.pl.heatmap(
            adata,
            var_names=gene_list,
            groupby="louvain",
            cmap="RdBu_r",
            standard_scale="var",
            dendrogram=False,
            show_gene_labels=True,
            show=False,
            figsize=(10, 8),
        )

    except Exception as e:
        raise RuntimeError(f"Erreur lors de la génération de la heatmap: {str(e)}")

    return fig