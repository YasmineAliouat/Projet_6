import sys
import os
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from gene_search.gene_search_utils import resolve_gene_to_var_name
from normalise_names.normalise_names import gene_versions, add_biomart_names

def prepare_adata_for_gene_search(adata, use_biomart=False, organism="hsapiens"):
    """
    Cette fonction permet de préparer adata pour la recherche de gènes sans obliger Streamlit
    à faire la normalisation lui-même.
    """
    if adata is None:
        return adata

    if "_gene_search_prepared" not in adata.uns:
        adata.uns["_gene_search_prepared"] = False
    if "_gene_search_biomart" not in adata.uns:
        adata.uns["_gene_search_biomart"] = False

    need_basic = (
        "base_name" not in adata.var.columns
        or "has_versions" not in adata.var.columns
    )

    need_biomart = use_biomart and not adata.uns["_gene_search_biomart"]

    if adata.uns["_gene_search_prepared"] and not need_basic and not need_biomart:
        return adata

    if need_basic:
        adata = gene_versions(adata)

    if need_biomart and "gene_symbol" not in adata.var.columns:
        try:
            adata = add_biomart_names(adata, organism=organism)
            adata.uns["_gene_search_biomart"] = True
        except Exception:
            pass

    adata.uns["_gene_search_prepared"] = True
    return adata

def resolve_gene_for_streamlit(adata,query,max_hits=20,use_biomart=False,organism="hsapiens"):    
    """
    Cette fonction permet d'adapter la recherche de gènes à l'interface Streamlit.
    Elle utilise la fonction 'resolve_gene_to_var_name' du script 'gene_search_utils'
    et retourne un dictionaire simple avec:
        - status : "exact", "partial", "multiple", "suggestions", "none"
        - var_name : str ou None
        - hits : DataFrame ou None
        - suggestions : list
    """

    if adata is None:
        return {
            "status": "none",
            "var_name": None,
            "hits": None,
            "suggestions": [],
            "match_type": None
        }

    #Netoyage de la requette
    query = str(query).strip()
    if query == "":
        return {
            "status": "none",
            "var_name": None,
            "hits": None,
            "suggestions": [],
            "match_type": None
        }

    #préparation de adata en mémoire pour la recherche
    adata = prepare_adata_for_gene_search(adata,use_biomart=use_biomart,organism=organism)

    var_name, hits, suggestions, match_type = resolve_gene_to_var_name(
        adata,
        query,
        choice=None,
        max_hits=max_hits
    )

    if var_name is not None:
        #Si recherche exacte:
        if match_type == "exact":
            return {
                "status": "exact",
                "var_name": var_name,
                "hits": hits,
                "suggestions": [],
                "match_type": match_type
            }

        #Si recharche partielle, on demande la confirmation utilisateur
        if match_type == "partial":
            return {
                "status": "multiple",
                "var_name": None,
                "hits": hits,
                "suggestions": [],
                "match_type": match_type
            }

    #Si plusieurs résulatst possible, l'utilisateur doit choisir
    if hits is not None and not hits.empty:
        return {
            "status": "multiple",
            "var_name": None,
            "hits": hits,
            "suggestions": [],
            "match_type": match_type
        }

    # Si on a pas de match direct, mais des suggestions proches
    if suggestions:
        return {
            "status": "suggestions",
            "var_name": None,
            "hits": None,
            "suggestions": suggestions,
            "match_type": match_type
        }

    #Si rien trouvé
    return {
        "status": "none",
        "var_name": None,
        "hits": None,
        "suggestions": [],
        "match_type": match_type
    }

def hits_to_options(hits):
    """
    Cette fonction permet de transformer le DataFrame de résultats (hits)
    en options lisibles pour un selectbox Streamlit.

    Chaque option est un tuple :
        (label_affiché, var_name)
    """

    if hits is None or hits.empty:
        return []
    options = []
    for _, row in hits.iterrows():
        #récupération des informations utiles
        var_name = row.get("__var_name__", "")
        gene_symbol = row.get("gene_symbol", "")
        gene_id = row.get("gene_ids", "")
        match_value = row.get("__match_value__", "")

        #on nettoie les valeurs manquantes
        if pd.isna(var_name):
            var_name = ""
        if pd.isna(gene_symbol):
            gene_symbol = ""
        if pd.isna(gene_id):
            gene_id = ""
        if pd.isna(match_value):
            match_value = ""

        var_name = str(var_name)
        gene_symbol = str(gene_symbol)
        gene_id = str(gene_id)
        match_value = str(match_value)

        #si pas de var_name on ignore
        if var_name == "":
            continue

        #construction du label affiché dans Streamlit
        label = var_name

        #ajout du symbole du gène
        if gene_symbol and gene_symbol != var_name:
            label += f" | {gene_symbol}"

        #ajout de l'identifiant Ensembl
        if gene_id and gene_id not in [var_name, gene_symbol]:
            label += f" | {gene_id}"

        #indique la valeur exacte qui a matché la requête
        if match_value and match_value not in [var_name, gene_symbol, gene_id]:
            label += f" (matched: {match_value})"
        options.append((label, var_name))

    return options
def suggestions_to_options(adata, suggestions, use_biomart=False, organism="hsapiens"):
    """
    Transforme une liste de suggestions textuelles en options lisibles pour Streamlit.
    On ajoute aussi le terme qui a servi à proposer ce gène.
    """
    options = []
    for sugg in suggestions:
        result = resolve_gene_for_streamlit(adata, sugg, max_hits=10, use_biomart=use_biomart, organism=organism)

        #Si la recherche est exacte
        if result["status"] == "exact" and result["var_name"] is not None:
            var_name = result["var_name"]
            hits = result["hits"]

            gene_symbol = ""
            gene_id = ""

            if hits is not None and not hits.empty:
                row = hits.iloc[0]
                if "gene_symbol" in hits.columns and row["gene_symbol"] == row["gene_symbol"]:
                    gene_symbol = row["gene_symbol"]
                if "gene_ids" in hits.columns and row["gene_ids"] == row["gene_ids"]:
                    gene_id = row["gene_ids"]

            #Construction du label si recherche exacte
            label = f"{var_name}"
            if gene_symbol and gene_symbol != var_name:
                label += f" | {gene_symbol}"
            if gene_id:
                label += f" | {gene_id}"

            label += f"  (matched: {sugg})"
            options.append((label, var_name))

        #Si recherche partiel et plusieurs résultats sont trouvé
        elif result["status"] == "multiple" and result["hits"] is not None:
            for _, row in result["hits"].iterrows():
                var_name = row["__var_name__"]
                gene_symbol = row["gene_symbol"] if "gene_symbol" in row and row["gene_symbol"] == row["gene_symbol"] else ""
                gene_id = row["gene_ids"] if "gene_ids" in row and row["gene_ids"] == row["gene_ids"] else ""

                #Construction du label si résutats multiples
                label = f"{var_name}"
                if gene_symbol and gene_symbol != var_name:
                    label += f" | {gene_symbol}"
                if gene_id:
                    label += f" | {gene_id}"

                label += f"  (matched: {sugg})"
                options.append((label, var_name))

    # enlever les doublons sur var_name
    unique_options = []
    seen = set()
    for label, var_name in options:
        if var_name not in seen:
            unique_options.append((label, var_name))
            seen.add(var_name)

    return unique_options