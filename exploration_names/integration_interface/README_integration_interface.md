# Integration_Interface :
## Description :
Ce module fournit une interface entre la logique de recherche de gènes et une application utilisateur (`Streamlit`).

Il permet de :

- Préparer automatiquement les données (AnnData).
- Gérer les différents types de résultats (exact, multiple, suggestions)
- Fournir des formats directement exploitables dans une interface graphique.

## Fonctions principales :
### `prepare_adata_for_gene_search`
Prépare adata en mémoire pour la recherche en ajoutant `base_name` et `has_versions` si absents et les annotations BioMart si demandé.

### `resolve_gene_for_streamlit`
Fonction principale pour la recherche de gènes. Prend en entrées `adata`, `query`,`max_hits`, `use_biomart` et `organism`.

**Sortie :**

```python
{
    "status": "exact" ou "multiple" ou "suggestions" ou "none",
    "var_name": str ou None,
    "hits": DataFrame ou None,
    "suggestions": list,
    "match_type": str ou None
}
```

### `hits_to_options`
Transforme les résultats en options lisibles pour une interface :

```python
(label, var_name)
```
### `suggestions_to_options`
Transforme les suggestions en options sélectionnables pour l’utilisateur.


> [!NOTE]
> - Les modifications sont faites uniquement en mémoire.  
> - BioMart est optionnel

## Exemple d'utilisation Streamlit:
```python
result = resolve_gene_for_streamlit(adata,query,max_hits=20,use_biomart=False,organism="hsapiens")

if result["status"] == "exact":
    selected_gene = result["var_name"]

elif result["status"] == "multiple":
    options = hits_to_options(result["hits"])
    selected_gene = st.selectbox("Choisir un gène", options)

elif result["status"] == "suggestions":
    options = suggestions_to_options(adata, result["suggestions"])
    selected_gene = st.selectbox("Suggestions", options)
```

## Dépendances

Ce module dépend de :
- `gene_search_utils.py`
- `prepare_names.py`


