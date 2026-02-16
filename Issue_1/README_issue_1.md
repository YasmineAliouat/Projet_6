# Exploration des fichiers AnnData (.h5ad)

Ce script permet d’analyser rapidement la structure d’un fichier AnnData (.h5ad) afin de comprendre l’organisation des données single-cell avant le développement de l’interface web.

Il est utilisé dans le cadre de l’Issue 1 : Exploration des données AnnData (.h5ad).

---

## Objectifs

Le script permet de :

- Charger un fichier .h5ad
- Identifier les principales structures (X, obs, var, raw, layers, obsm)
- Vérifier la présence d’une UMAP
- Examiner les métadonnées disponibles
- Estimer si les données sont log-transformées
- Fournir un aperçu des gènes et des cellules

---

## Prérequis

- Python ≥ 3.9
- Bibliothèques :
  - scanpy
  - anndata
  - numpy
  - scipy
  - pandas

Installation possible via :

```bash
pip install scanpy anndata numpy scipy pandas
