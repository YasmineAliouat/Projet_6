# Indicateurs biologiques visualisables

Ce document décrit les principaux indicateurs biologiques accessibles dans l’interface et leur mode de visualisation.

Il correspond au livrable de l’Issue 3 : Définition des indicateurs biologiques.

---

## Objectif

L’objectif de l’interface est de permettre aux biologistes d’explorer facilement des données single-cell stockées sous format AnnData, sans nécessiter de compétences en programmation.

Les indicateurs biologiques sélectionnés répondent aux principales questions posées lors de l’analyse des données transcriptomiques.

---

## 1. Expression d’un gène

### Définition
L’expression d’un gène correspond au niveau de transcription mesuré dans chaque cellule.

Elle permet de déterminer :
- si un gène est exprimé,
- dans quelles cellules,
- avec quelle intensité.

### Source des données
- Matrice principale : `adata.X`
- Données brutes éventuelles : `adata.raw.X`

En général, es valeurs sont normalisées et log-transformées.

### Visualisation
- Histogramme / violin plot : distribution globale
- Projection UMAP colorée par niveau d’expression

### Question biologique associée
> Le gène est-il exprimé ?  
> Dans quelles populations cellulaires ?

---

## 2. Co-expression de deux gènes

### Définition
La co-expression correspond à l’expression simultanée de deux gènes dans une même cellule.

Elle permet d’identifier des relations fonctionnelles potentielles entre gènes.

### Source des données
- Valeurs extraites depuis `adata.X`

### Visualisation
- UMAP avec code couleur spécifique pour :
  - cellules exprimant les deux gènes ;
  - cellules exprimant un seul gène ;
  - cellules n’exprimant aucun des deux

### Question biologique associée
> Les cellules exprimant le gène A expriment-elles aussi le gène B ?

---

## 3. Signatures géniques

### Définition
Une signature génique est un ensemble de gènes dont l’expression permet d’identifier un type de cellule.

### Source des données
- Liste de gènes fournie par l’utilisateur
- Expression extraite depuis `adata.X`

### Méthode de calcul
Pour chaque cellule, un score est calculé à partir de l’expression moyenne (ou pondérée) des gènes de la signature.

### Visualisation
- UMAP colorée selon le score de signature
- Histogramme des scores

### Question biologique associée
> Quelles cellules ont une expression génique commune ?

---

## 4. Projection UMAP

### Définition
L’UMAP est une méthode de réduction de dimension permettant de représenter les cellules dans un espace bidimensionnel.

### Source des données
- `adata.obsm["X_umap"]` si disponible
- Calculée automatiquement sinon

### Visualisation
- Nuage de points 2D
- Coloration selon l’indicateur sélectionné

### Question biologique associée
> Comment s’organisent les cellules selon leur profil transcriptomique ?

---

## 5. Métadonnées expérimentales

### Définition
Les métadonnées décrivent le contexte expérimental des cellules.

### Source des données
- `adata.obs`

### Visualisation
- Coloration UMAP
- Filtres

### Question biologique associée
> Les résultats sont-ils influencés par des facteurs techniques ou expérimentaux ?

---

## 6. Indicateurs optionnels (évolutions possibles)

Selon l’avancement du projet, d’autres indicateurs pourront être intégrés.

### Pseudo-temps
- Inférence de trajectoires cellulaires
- Source : Scanpy / outils dédiés

### Vélocité ARN
- Dynamique transcriptionnelle
- Source : layers spliced / unspliced

Ces fonctionnalités ne sont pas prioritaires dans la version MVP.

---

## Synthèse

| Indicateur        | Source         | Visualisation       | Priorité |
|-------------------|----------------|---------------------|----------|
| Expression gène   | adata.X / raw  | Histogramme + UMAP  | Haute    |
| Co-expression     | adata.X        | UMAP                | Haute    |
| Signature génique | adata.X        | UMAP + histogramme  | Moyenne  |
| UMAP              | adata.obsm     | Scatter plot        | Haute    |
| Métadonnées       | adata.obs      | UMAP + filtres      | Moyenne  |

---

## Lien avec le projet

Ces indicateurs constituent la base fonctionnelle de l’interface.

Ils répondent directement aux besoins exprimés par les biologistes utilisateurs et orientent les choix de développement pour les phases suivantes.
