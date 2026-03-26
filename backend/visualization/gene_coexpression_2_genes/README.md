# Gene_coexpression_2_genes

## Description

Script Python permettant de visualiser la **co‑expression de deux gènes** à partir d’un objet **AnnData (Scanpy)**.

Ce module génère des visualisations interactives utilisées en analyse **single-cell RNA-seq** pour explorer les relations entre deux gènes au sein des clusters (ou pas).

## Objectifs

La fonction `plot_gene_coexpression()` permet de générer :

- Scatter plot interactif avec coloration par cluster
- Histogrammes marginaux de distribution
- Calcul du coefficient de corrélation de Pearson par cluster
- UMAP catégorisé selon la co‑expression (aucun / gène A / gène B / les deux)

## Dépendances

- `numpy`
- `scanpy`
- `pandas`
- `plotly`
- `streamlit`
- `scipy`

Installation recommandée (si nécessaire) :

```bash
pip install scanpy numpy pandas plotly streamlit scipy
```

## Fonction principale

```python
plot_gene_coexpression(adata, genes, plot_types)
```

### Paramètres

| Paramètre | Type | Description |
|-----------|------|------------|
| `adata` | AnnData | Objet AnnData contenant les données single-cell |
| `genes` | list[str] | Liste contenant exactement deux gènes `[gene_a, gene_b]` |
| `plot_types` | list[str] | Liste des visualisations à générer |

### Pré-requis sur les données

L’objet `adata` doit contenir :

- Coordonnées UMAP (`adata.obsm["X_umap"]`)
- Clusters Louvain (`adata.obs["louvain"]`)
- Données d’expression normalisées
- Les deux gènes présents dans `adata.var_names`

### Exemple d’utilisation

```python
figs = plot_gene_coexpression(
    adata=adata,
    genes=["CD3E", "CD8A"],
    plot_types=["scatter", "umap"]
)
```

## Retour de la fonction

La fonction retourne :

```python
List[Figure]
```

Liste des figures générées :

- `plotly.graph_objects.Figure` pour le scatter plot
- `matplotlib.figure.Figure` pour l’UMAP
