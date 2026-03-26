# Backend

Le dossier `backend/` regroupe l’ensemble des modules responsables du traitement des données, de la gestion des noms de gènes et de la génération des visualisations.

## Structure

Le backend est organisé en trois modules principaux :

- `data_loading/`  
  Chargement, validation et résumé des fichiers AnnData (`.h5ad`)

- `gene_names/`  
  Gestion des noms de gènes (normalisation, recherche, suggestions)

- `visualization/`  
  Fonctions de visualisation (expression, co-expression, signatures)

## Rôle dans l’application

Le backend est utilisé par le frontend (Streamlit) pour :
- charger et préparer les données
- interpréter les requêtes utilisateur
- produire les graphiques
