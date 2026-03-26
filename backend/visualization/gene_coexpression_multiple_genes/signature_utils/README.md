# Signature_utils

## Description

Script Python permettant de calculer et visualiser un **score de signature génique** à partir d’un objet **AnnData (Scanpy)**.

Ce module permet d’évaluer l’enrichissement d’un ensemble de gènes (signature biologique) au niveau **cellulaire**, en analyse **single-cell RNA-seq**.

---

## Objectifs

Les fonctions permettent de :

- Charger une signature génique depuis un fichier texte
- Calculer un score de signature par cellule (`scanpy.tl.score_genes`)
- Ajouter le score dans `adata.obs`
- Visualiser le score sur un UMAP

## Dépendances

- `scanpy`
- `matplotlib`

Installation recommandée (si nécessaire) :

```bash
pip install scanpy matplotlib
```

## Fonctions principales

### 1. Charger une signature

```python
load_signature_genes(path)
```

#### Paramètres

| Paramètre | Type | Description |
|-----------|------|------------|
| `path` | str | Chemin vers un fichier texte contenant un gène par ligne |

#### Retour

```python
list[str]
```

Liste des gènes composant la signature.


### 2. Calculer un score de signature

```python
compute_signature_score(adata, gene_list, score_name)
```

#### Paramètres

| Paramètre | Type | Description |
|-----------|------|------------|
| `adata` | AnnData | Objet AnnData contenant les données single-cell |
| `gene_list` | list[str] | Liste des gènes composant la signature |
| `score_name` | str | Nom du score ajouté dans `adata.obs` |

#### Fonctionnement

- Vérifie quels gènes de la signature sont présents dans `adata.var_names`
- Ignore automatiquement les gènes absents
- Lève une erreur si aucun gène n’est présent
- Calcule le score avec `scanpy.tl.score_genes()`
- Ajoute le score dans :

```python
adata.obs[score_name]
```

#### Retour

```python
AnnData
```

Objet `adata` enrichi avec le score de signature.

### 3. Visualiser le score sur UMAP

```python
plot_signature_score(adata, gene_list, score_name="signature_score")
```

#### Paramètres

| Paramètre | Type | Description |
|-----------|------|------------|
| `adata` | AnnData | Objet AnnData |
| `gene_list` | list[str] | Liste des gènes composant la signature |
| `score_name` | str | Nom du score (par défaut `"signature_score"`) |

#### Visualisation générée

- UMAP coloré selon le score de signature
- Échelle continue représentant l’enrichissement

#### Retour

```python
matplotlib.figure.Figure
```

Figure correspondant à l’UMAP coloré par score.

## Exemple d’utilisation

```python
genes = load_signature_genes("tcell_signature.txt")

fig = plot_signature_score(
    adata=adata,
    gene_list=genes,
    score_name="Tcell_score"
)
```