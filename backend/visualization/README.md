# Visualization Modules

## Description

Ce projet regroupe plusieurs modules de visualisation pour l’analyse **single-cell RNA-seq** à partir d’un objet **AnnData (Scanpy)**.

Les outils permettent d’explorer :

- L’expression d’un gène unique  
- La co‑expression de deux gènes  
- Les signatures géniques  
- Les heatmaps multi‑gènes  

Chaque fonctionnalité est organisée dans un dossier dédié avec son propre README détaillé.

## Modules disponibles

### 1. Gene_expression

Visualisation de l’expression d’un gène :

- UMAP des clusters
- UMAP d'expression du gène
- Violin plot
- Histogramme

Fonction principale :

```python
plot_gene_expression()
```

### 2. Gene_coexpression_2_genes

Analyse conjointe de deux gènes :

- Scatter plot interactif
- Corrélation de Pearson par cluster
- Histogrammes marginaux
- UMAP catégorisé selon co‑expression

Fonction principale :

```python
plot_gene_coexpression()
```

### 3. Gene_coexpression_multiple_genes 

#### Signature_utils

Calcul et visualisation d’un score basé sur une liste de gènes :

- Calcul via `scanpy.tl.score_genes`
- Ajout du score dans `adata.obs`
- Visualisation UMAP continue

Fonctions principales :

```python
load_signature_genes()
compute_signature_score()
plot_signature_score()
```

#### Heatmap

Visualisation comparative de plusieurs gènes :

- Regroupement par cluster Louvain
- Normalisation par gène
- Heatmap d’expression relative

Fonction principale :

```python
coexp_heatmap()
```

## Dépendances principales

- `scanpy`
- `numpy`
- `pandas`
- `matplotlib`
- `plotly`
- `streamlit`
- `scipy`

Installation recommandée :

```bash
pip install scanpy numpy pandas matplotlib plotly streamlit scipy
```

## Pré-requis sur les données

L’objet `adata` doit contenir :

- Coordonnées UMAP (`adata.obsm["X_umap"]`)
- Clusters Louvain (`adata.obs["louvain"]`)
- Données d’expression normalisées (idéalement log-transformées)
- Gènes présents dans `adata.var_names`

## Structure du projet

```
visualization/
│
├── gene_expression/
├── gene_coexpression_2_genes/
├── gene_coexpression_multiple_genes/
        └── gene_signature/
        └── heatmap/
```

Chaque dossier contient :

- Le script Python
- Un README détaillé
- Les spécifications des fonctions


Pour les détails techniques (paramètres, types de retour, options avancées), se référer aux README spécifiques de chaque module.