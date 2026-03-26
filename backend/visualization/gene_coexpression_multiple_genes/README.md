# Gene_coexpression_multiple_genes

## Description

Ce dossier contient des modules permettant :

- Le calcul d’un **score de signature génique**
- La visualisation multi‑gènes via **heatmap**
  
Ces outils sont conçus pour l’analyse **single-cell RNA-seq** à partir d’un objet **AnnData (Scanpy)**.

## Contenu du dossier

### 1. Signature génique

Permet de :

- Charger une liste de gènes depuis un fichier texte
- Calculer un score de signature par cellule (`scanpy.tl.score_genes`)
- Ajouter le score dans `adata.obs`
- Visualiser le score sur UMAP

Fonctions principales :

```python
load_signature_genes()
compute_signature_score()
plot_signature_score()
```

Retourne :

- `list[str]`
- `AnnData`
- `matplotlib.figure.Figure`

### 2. Heatmap de co‑expression

Permet de :

- Visualiser plusieurs gènes simultanément
- Regrouper les cellules par cluster Louvain
- Comparer les niveaux d’expression relatifs

Fonction principale :

```python
coexp_heatmap()
```

Retourne :

```python
matplotlib.figure.Figure
```

## Dépendances

- `scanpy`
- `matplotlib`
- `numpy`

Installation :

```bash
pip install scanpy matplotlib numpy
```

## Pré-requis

L’objet `adata` doit contenir :

- Coordonnées UMAP (`adata.obsm["X_umap"]`)
- Clusters Louvain (`adata.obs["louvain"]`)
- Données d’expression normalisées (idéalement log-transformées)
- Les gènes présents dans `adata.var_names`


Les détails techniques de chaque fonction sont disponibles dans les README spécifiques des modules correspondants.