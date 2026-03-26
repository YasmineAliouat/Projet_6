# Gene_expression

## Description

Script Python permettant de visualiser l’expression d’un gène à partir d’un objet **AnnData (Scanpy)**.

Ce module génère différents types de visualisations utilisées en analyse **single-cell RNA-seq**.

## Objectifs

La fonction `plot_gene_expression()` permet de générer :

- UMAP coloré par clusters (Louvain)
- UMAP coloré par expression d’un gène
- Violin plot par cluster
- Histogramme de distribution d’expression

## Dépendances

- `numpy`
- `scanpy`
- `scipy`
- `matplotlib`
- `streamlit`

Installation recommandée (si nécessaire) :

```bash
pip install scanpy numpy scipy matplotlib streamlit
```

## Fonction principale

```python
plot_gene_expression(adata, gene, plot_types)
```

### Paramètres

| Paramètre | Type | Description |
|-----------|------|------------|
| `adata` | AnnData | Objet AnnData contenant les données single-cell |
| `gene` | str | Nom du gène à visualiser |
| `plot_types` | list[str] | Liste des visualisations à générer |

### Pré-requis sur les données

L’objet `adata` doit contenir :

- Coordonnées UMAP (`adata.obsm["X_umap"]`)
- Clusters Louvain (`adata.obs["louvain"]`)
- Données normalisées (log1p recommandé)

### Exemple d’utilisation

```python
figs = plot_gene_expression(
    adata=adata,
    gene="MKI67",
    plot_types=["umap_", "violin", "histogram"]
)
```

### Retour de la fonction

La fonction retourne :

```python
List[matplotlib.figure.Figure]
```

Liste des figures générées :

- `matplotlib.figure.Figure` pour l’UMAP par clusters
- `matplotlib.figure.Figure` pour l’UMAP par expression d'un gène
- `matplotlib.figure.Figure` pour le violin plot
- `matplotlib.figure.Figure` pour l'histogramme

