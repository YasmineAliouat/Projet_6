# Application Web pour l’Exploration de Données Single-Cell

## Description
Ce projet vise à développer une interface web interactive permettant l’exploration et la visualisation de données de transcriptomique single-cell au format AnnData (`.h5ad`).

Il s’inscrit dans un contexte de recherche en biologie des systèmes appliquée aux neuroblastomes, et a pour objectif de rendre ces données accessibles à des utilisateurs non informaticiens, notamment des biologistes.


## Objectifs
L'application permet de :
- Visualiser les clusters mis en évidence par les manip
- Visualiser l’expression d’un gène
- Observer la distribution d’expression
- Visualiser les cellules dans l’espace UMAP
- Analyser la co-expression de deux ou plusieurs gènes
- Explorer des signatures géniques
- Gérer intelligemment les noms de gènes (fautes, alias, suggestions)

## Fonctionalités principales

### Expression d'un gène :
- Expression normalisée dans l'espace UMAP
- Violin plot
- Histogramme de distrubtion

### Co-expression de gènes :
- Scatter plot (corrélation)
- Co-expression dans l'espace UMAP
- Heatmap 

### Signature de gènes :
- Score de signature (moyenne d'expression des gènes)
- Projection du score dans l'espace UMAP
- Visualisation globale de l'activité de la signature
- Heatmap de co-expression des gènes de la signature
- Sélection et validation des gènes composant la signature 

### Gestions des noms de gènes :
- Enrechissement du Anndata (Biomart)
- Insensible à la casse
- Tolérance aux fautes de frappe
- Suggestions automatiques


## Structure du projet
projet/   
├── frontend/ # Interface Streamlit  
│ └── app.py  
│   
├── data/  
│ └── adata_3583.h5ad  
│ └── adata_3583_new.h5ad  
│   
├── backend/  
│   ├── data_loading/     # Chargement et validation des données  
│   ├── gene_names/       # Recherche et normalisation des gènes  
│   ├── visualization/    # Fonctions de visualisation   
│  
├── Dockerfile  
├── .dockerignore   
├── requirements.txt  
└── README.md

## Installation
### 1. Cloner le projet 
```bash
git clone <url_du_projet>
cd projet-6
```
### 2. Créer un environnement virtuel
```bash
python -m venv venv
source venv/bin/activate
```
### 3. Installer des dépendances
```bash
pip install -r requirements.txt
```
### 4. Lancer l'application
```bash
streamlit run frontend/app.py
```

## Données
L’application nécessite un fichier AnnData (`.h5ad`).  
Deux options sont possibles pour l'emplacement du fichier :  


**En local**  
Entrer un chemin dans l’interface :
```
/chemin/vers/fichier.h5ad
```
**Dossier `data/`**
```
data/fichier.h5ad
```

Deux types de jeux de données peuvent être utilisés :  

**Données brutes**  
```
data/adata_3583.h5ad
```
- Données initiales
- Noms de gènes non normalisés

**Données enrichies**  
```
data/adata_3583_new.h5ad
```
- Données enrichies via Biomart
- Meilleure compatibilité avec l’interface
- Support de la recherche intelligente (correction, suggestions)

> [!TIP]
> - Il est fortement conseillé d’utiliser le fichier enrichi `data/adata_3583_new.h5ad` afin de bénéficier de toutes les fonctionnalités de l’application.

## Utilisation avec Docker
### 1. Installer Docker
Si Docker n'est pas installé sur votre machine, pour pouvoir utiliser l'application avec Docker , il faut l'installer en suivant les instrustions ici:
- Linux : https://docs.docker.com/engine/install/ubuntu/
- Windows / Mac : https://www.docker.com/products/docker-desktop/

**Vérification**
```bash
docker --version
```

### 2. Construire l'image du projet
```bash
docker build -t projet6-app .
```
### 3. Lancer le conteneur
```bash
docker run -p 8501:8501  -v /chemin/vers/donnees:/app/data  projet6-app
```
### 4. Dans l'application, utiliser :
```
data/file_name.h5ad
```

## Contexte académique :
Projet réalisé dans le cadre du Master Bioinformatique (Université Lyon 1).

### Equipe
- **Etudiants** :
  - ALIOUAT Yasmine - yasmine.aliouat@etu.univ-lyon1.fr
  - CAUCHOIS Fanny - fanny.cauchois@etu.univ-lyon1.fr
  - MERAL Duzguncan - duzguncan.meral@etu.univ-lyon1.fr
  - NAVARRO Janice - janice.navarro@etu.univ-lyon1.fr
- **Maîtres d'Ouvrage** :
  - GANDRILLON Olivier (CNRS) - olivier.gandrillon@ens-lyon.fr
- **Responsable Pédagogique** :
  - MARY Arnaud - arnaud.mary@univ-lyon1.fr
