# Expression d'un gène - visualisations

Ce document décrit la fonctionnalité permettant de visualiser l’expression d’un gène unique dans un jeu de données single-cell au format AnnData.

Il correspond au livrable de l’Issue 12 : Visualisation de l’expression d’un gène.

---

## Objectif

Cette fonctionnalité permet de répondre aux questions biologiques suivantes :

- Est-ce que ce gène est exprimé ?
- Si oui, dans quelles cellules ?

---

## Visualisations

La fonction `plot_gene_expression()` génère trois représentations :

### Violin plot  
Distribution globale de l’expression sur l’ensemble des cellules.  
Permet d’évaluer rapidement la sparsité et la variabilité.

### Histogramme  
Distribution quantitative des niveaux d’expression.  
Met en évidence la proportion de cellules à expression nulle ou faible.

### UMAP colorée  
Projection bidimensionnelle des cellules, colorées selon leur niveau d’expression.  
Permet d’identifier les sous-populations exprimantes.

---

## Test

Un script Python `test_plot.py` est fourni dans le dépôt.

Ce script permet de tester la fonction `plot_gene_expression()` indépendamment de l’interface web.

Il :
- charge un fichier .h5ad
- demande un nom de gène
- appelle la fonction
- affiche les visualisations

Il sert à valider le bon fonctionnement de la fonctionnalité pendant le développement.
