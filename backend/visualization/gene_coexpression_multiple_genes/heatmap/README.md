# Heatmap

## Description

Script Python permettant de visualiser la **co‑expression de plusieurs gènes** sous forme de **heatmap** à partir d’un objet **AnnData (Scanpy)**.

Ce module génère une heatmap regroupant les cellules par cluster Louvain afin d’explorer les **profils d’expression relatifs entre gènes** en analyse **single-cell RNA-seq**.

## Objectifs

La fonction `coexp_heatmap()` permet de générer :

- Heatmap d’expression de plusieurs gènes
- Regroupement des cellules par cluster Louvain
- Normalisation par gène (`standard_scale="var"`)
- Visualisation comparative des niveaux relatifs d’expression

## Dépendances

- `scanpy`
- `matplotlib`
- `typing`

Installation recommandée (si nécessaire) :

```bash
pip install scanpy matplotlib
```

## Fonction principale

```python
coexp_heatmap(adata, gene_list)
```

### Paramètres

| Paramètre | Type | Description |
|-----------|------|------------|
| `adata` | AnnData | Objet AnnData contenant les données single-cell |
| `gene_list` | List[str] | Liste des gènes à afficher dans la heatmap |

### Pré-requis sur les données

L’objet `adata` doit contenir :

- Clusters Louvain (`adata.obs["louvain"]`)
- Les gènes présents dans `adata.var_names`
- Données d’expression normalisées
- Idéalement données log-transformées (`log1p`)

### Détails de la visualisation

La heatmap est générée avec :

- `groupby="louvain"` → regroupement par cluster
- `cmap="RdBu_r"` → palette divergente
- `standard_scale="var"` → normalisation par gène
- `dendrogram=False` → pas de clustering hiérarchique
- `figsize=(10, 8)`

### Exemple d’utilisation

```python
fig = coexp_heatmap(
    adata=adata,
    gene_list=["CD3E", "CD8A", "GZMB", "MKI67"]
)
```

### Retour de la fonction

La fonction retourne :

```python
matplotlib.figure.Figure
```

Figure correspondant à la heatmap générée par `scanpy.pl.heatmap()`.