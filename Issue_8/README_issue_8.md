# Chargement et gestion des fichiers AnnData

Ce document décrit la fonctionnalité de chargement et de validation des fichiers single-cell au format AnnData (.h5ad).

Il correspond au livrable de l'Issue 8 : Chargement et gestion des fichiers AnnData.

---

## Objectif

Cette fonctionnalité permet de :
- Charger des fichiers .h5ad
- Valider leur structure et leur contenu
- Extraire des métadonnées utiles pour l'analyse et l'interface utilisateur

Elle répond aux besoins techniques suivants :
- Vérifier l'intégrité des données avant traitement
- Identifier les caractéristiques clés du jeu de données
- Préparer les données pour les visualisations et analyses ultérieures

---

## Fonctionnalités principales

### Chargement du fichier
La fonction `load_anndata()` :
- Lit un fichier .h5ad avec Scanpy
- Gère les erreurs de lecture
- Retourne un objet AnnData prêt à l'emploi

### Validation des données
La fonction `validate_anndata()` vérifie :
- La présence de la matrice d'expression (X)
- Les dimensions du jeu de données (n_obs, n_vars)
- L'unicité des noms de gènes (var_names)
- L'unicité des identifiants de cellules (obs_names)

### Résumé des métadonnées
La fonction `summarize_anndata()` extrait :
- Nombre de cellules (n_obs) et gènes (n_vars)
- Format de la matrice (sparse/dense)
- Présence de données brutes (raw)
- Liste des couches disponibles (layers)
- Aperçu des métadonnées (obs/var)
- Clés des réductions dimensionnelles (obsm)
- Présence d'une projection UMAP

---

## Structure de données

Le résumé est encapsulé dans la classe `AnnDataSummary` qui contient :
- `n_obs` : Nombre de cellules
- `n_vars` : Nombre de gènes
- `is_sparse_X` : Format de la matrice d'expression
- `has_raw` : Présence de données brutes
- `layers` : Liste des couches disponibles
- `obs_columns_preview` : Aperçu des métadonnées cellules
- `var_columns_preview` : Aperçu des métadonnées gènes
- `obsm_keys` : Clés des réductions dimensionnelles
- `has_umap` : Présence d'une projection UMAP

---

## Test

Un script Python `test_anndata.py` est fourni dans le dépôt.

Ce script permet de tester les fonctionnalités de chargement et de validation :
- Charge un fichier .h5ad de test
- Valide sa structure
- Affiche un résumé des métadonnées
- Gère les cas d'erreur (fichier corrompu, structure invalide)

Il sert à valider le bon fonctionnement de la fonctionnalité pendant le développement et avant intégration dans l'interface web.
