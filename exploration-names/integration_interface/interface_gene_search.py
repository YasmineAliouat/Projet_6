import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from gene_search.gene_search_utils import resolve_gene_to_var_name


def resolve_gene_for_streamlit(adata, query, max_hits=20):
    """
    Cette fonction permet d'adapter la recherche de gènes à l'interface Streamlit.
    Elle utilise la fonction 'resolve_gene_to_var_name' du script 'gene_search_utils'
    et retourne un dictionaire simple avec:
        - status : "exact", "partial", "multiple", "suggestions", "none"
        - var_name : str ou None
        - hits : DataFrame ou None
        - suggestions : list
    """
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
    Cette fonction transforme les résultats en options lisibles pour un selectbox Streamlit.
    liste de tuples du style :
            ("AL390719.1 | ENSG00000217801", "AL390719.1")
    """
    
    if hits is None or hits.empty:
        return []

    options = []
    # On parcourt les résultats pour fabriquer un label lisible qui contient:
    for _, row in hits.iterrows():
        # 1-le nom exact dans adata.var_names
        var_name = row["__var_name__"]
        # 2-l'identifiant Ensembl
        gene_id = row["gene_ids"] if "gene_ids" in row and row["gene_ids"] == row["gene_ids"] else ""
        # 3-le symbole du gène
        gene_symbol = row["gene_symbol"] if "gene_symbol" in row and row["gene_symbol"] == row["gene_symbol"] else ""
        match_value = row["__match_value__"] if "__match_value__" in row and row["__match_value__"] == row["__match_value__"] else ""

        label = f"{var_name}"

        if gene_symbol and gene_symbol != var_name:
            label += f" | {gene_symbol}"

        if gene_id:
            label += f" | {gene_id}"

        if match_value and match_value != var_name and match_value != gene_symbol and match_value != gene_id:
            label += f" (matched: {match_value})"

        options.append((label, var_name))
    return options

def suggestions_to_options(adata, suggestions):
    """
    Transforme une liste de suggestions textuelles en options lisibles pour Streamlit.
    On ajoute aussi le terme qui a servi à proposer ce gène.
    """
    options = []
    for sugg in suggestions:
        result = resolve_gene_for_streamlit(adata, sugg, max_hits=10)

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