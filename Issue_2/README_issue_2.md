# Définition des cas d’usage et fonctionnalités MVP  

Ce document définit les cas d’usage principaux de l’interface web ainsi que la priorisation des fonctionnalités pour la version minimale viable (MVP).

Il correspond au livrable de l’Issue 2 : Définition des cas d’usage et fonctionnalités MVP.

---

## Objectif  

L’objectif est de cadrer précisément les besoins fonctionnels de l’outil avant le développement, afin de garantir une interface cohérente, utile et adaptée aux utilisateurs cibles.

L’interface doit permettre l’exploration intuitive de données single-cell au format AnnData, sans nécessiter de programmation.

---

## Utilisateurs cibles  

- Biologistes expérimentaux  
- Chercheurs en biologie / bio-informatique  
- Étudiants en formation  
- Toute personne souhaitant explorer un fichier AnnData sans écrire de code  

---

## 1. Explorer l’expression d’un gène  

### Questions biologiques :
- Ce gène est-il exprimé ?
- Dans quelles cellules ?

### Résultats attendus :
- Violin plot
- Histogramme
- UMAP colorée par expression

---

## 2. Analyser la co-expression de deux gènes  

### Questions biologiques :
- Les cellules exprimant le gène A expriment-elles aussi le gène B ?
- Existe-t-il une exclusion réciproque ?

### Résultats attendus :
- Catégorisation des cellules (none / A_only / B_only / both)
- UMAP colorée par catégorie
- Scatter plot A vs B

---

## 3. Visualiser une signature génique  

### Questions biologiques :
- Quelles cellules expriment une signature génique donnée ?
- Si oui, possibilité de couleur coder chaque cellule en fonction de son intensité à exprimer cette signature ?

### Résultats attendus :
- Score de signature par cellule
- UMAP colorée selon l’intensité du score

---

## Priorisation des fonctionnalités (MVP)

### Fonctionnalités essentielles (obligatoires)

1. Chargement d’un fichier AnnData  
2. Visualisation de l’expression d’un gène  
3. UMAP colorée  
4. Analyse de co-expression de deux gènes  

### Fonctionnalités secondaires (si temps disponible)

1. Analyse de signature génique  
2. Paramétrage des seuils d’expression  
3. Personnalisation des palettes de couleurs  
4. Export des figures  

---

## Vision globale du MVP  

Le MVP doit permettre à un biologiste de :

1. Charger un fichier AnnData  
2. Rechercher un gène  
3. Visualiser son expression  
4. Tester une co-expression  
5. Identifier des sous-populations cellulaires  

Le tout via une interface simple et interactive.