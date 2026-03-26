# Module `gene_names`

## Rôle du module

Ce module regroupe les outils liés à la gestion des noms de gènes dans les jeux de données AnnData.

Il permet :
- De préparer les noms de gènes pour la recherche,
- D’enrichir les annotations avec BioMart,
- De rechercher un gène à partir de différents identifiants,
- De proposer des suggestions en cas de faute de frappe ou de nom incomplet,
- D’intégrer ces résultats dans l’interface Streamlit.

L’objectif est de rendre la recherche de gènes plus robuste et plus ergonomique pour l’utilisateur.



## Organisation du dossier

### `prepare_names/`
Prépare les noms de gènes dans un objet AnnData en ajoutant des colonnes utiles comme :
- `base_name`
- `has_versions`
- annotations BioMart (optionnel)

Ce sous-module est utilisé pour enrichir les données avant ou pendant leur utilisation dans l’application.

### `gene_search/`
Contient la logique de recherche de gènes :
- recherche exacte,
- recherche partielle,
- gestion des versions,
- suggestions de noms proches.

### `integration_interface/`
Fait le lien entre la logique de recherche de gènes et l’interface Streamlit.

Il transforme les résultats en formats directement exploitables dans l’application :
- résultat exact,
- choix parmi plusieurs correspondances,
- suggestions.

### `adata_exploration/`
Contient des scripts d’exploration permettant d’inspecter la structure des noms de gènes et les annotations disponibles dans un fichier AnnData.

Ce sous-module sert surtout à l’analyse et au développement, pas directement à l’utilisation finale de l’interface.



## Utilisation dans l’application

Dans l’application web :
- les données sont préparées avec `prepare_names`,
- la recherche utilisateur passe par `integration_interface`,
- qui s’appuie lui-même sur les fonctions de `gene_search`.

Schéma simplifié :
```
AnnData -> prepare_names -> gene_search -> integration_interface -> Streamlit
```

## Dépendances principales

Ce module repose principalement sur :
- `scanpy`
- `pandas`
- `pybiomart` (si enrichissement BioMart)
- les annotations présentes dans `adata.var`



> [!NOTE]
> - Les annotations BioMart sont optionnelles mais améliorent fortement la qualité de la recherche.
> - Le fichier AnnData original n’est pas modifié sauf en cas de sauvegarde explicite.
