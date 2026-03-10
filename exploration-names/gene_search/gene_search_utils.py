import pandas as pd
import difflib
import re

#Colonnes dans lesquelles on peut chercher un gène.
DEFAULT_SEARCH_COLS = [
    "gene_ids",
    "gene_symbol",
    "hgnc_symbol",
    "alias_symbol",
    "entrez_id",
    "refseq_mrna",
    "uniprot_swissprot",
    "base_name",
]


def ensure_versions_cols(adata):
    """
    Cette fonction ajoute les colonnes base_name et has_versions
    elles sont absentes.
    - base_name : var_name sans suffixe de version (.X)
    - has_versions : True si var_name finit par (.X)
    """
    if "base_name" not in adata.var.columns or "has_versions" not in adata.var.columns:
        names = pd.Index(adata.var_names.astype(str))
        adata.var["base_name"] = names.str.replace(r"\.\d+$", "", regex=True)
        adata.var["has_versions"] = names.str.contains(r"\.\d+$", regex=True)
    return adata


def _as_str(s: pd.Series) -> pd.Series:
    """
    Cette fonction convertit les valeurs en str et remplace les valeurs 
    manquantes par une chaine vide pour faciliter la recherche.
    """
    return s.fillna("").astype(str)


def _match_alias_token(alias_series: pd.Series, query: str) -> pd.Series:
    """
    Cette fonction permet de rechercher un alias exact dans alias_symbol qui peut contenir "A|B|C"
    On cherche query comme token entier (entre |), insensible à la casse.
    """
    a = _as_str(alias_series)
    q = str(query).strip()
    if q == "":
        return pd.Series([False] * len(a), index=a.index)
    #pattern token exact : (^|\|)QUERY(\||$)
    pattern = rf"(?:^|\|){re.escape(q)}(?:\||$)"
    return a.str.contains(pattern, case=False, regex=True)

def get_name_columns(adata, keywords=None):
    """
    Cette fonction récupère uniquement les colonnes de adata.var qui ressemblent
     à des colonnes de noms ou d’identifiants.
    """
    if keywords is None:
        keywords = ("gene", "hgnc", "alias", "synonym", "entrez", "refseq", "uniprot", "ensembl", "symbol", "name", "id", "base")
    cols = [c for c in adata.var.columns if any(k in c.lower() for k in keywords)]
    return cols


def search_gene_hits(adata, query, search_cols=None, max_hits=50):
    """
    Cette fonction permet de Recherche un gène dans différentes colonnes d’annotation 
    (Ensembl, symbole, alias, UniProt, etc.) et retourne les résultats correspondants
    sous forme d'un DataFrame de hits (index = var_names).
    Cherche par :
    - ENSG (gene_ids)
    - var_names exact
    - base_name exact (avec versions)
    - gene_symbol / hgnc_symbol / entrez_id / refseq_mrna / uniprot_swissprot (exact)
    - alias_symbol (token exact)
    - sinon contains et suggestions proches
    """
    adata = ensure_versions_cols(adata)
    var = adata.var.copy()
    var["__var_name__"] = var.index.astype(str)

    q = str(query).strip()
    if q == "":
        return var.iloc[0:0]

    if search_cols is None:
        search_cols = DEFAULT_SEARCH_COLS

    cols_present = [c for c in search_cols if c in var.columns]
    view_cols = ["__var_name__"] + cols_present

    #ENSG exact
    if "gene_ids" in var.columns and q.startswith("ENSG"):
        hits = var[_as_str(var["gene_ids"]) == q]
        if not hits.empty:
            return hits[view_cols].head(max_hits)

    #var_names exact
    if q in adata.var_names:
        return var.loc[[q], view_cols].head(max_hits)

    #base_name exact
    if "base_name" in var.columns:
        hits = var[_as_str(var["base_name"]) == q]
        if not hits.empty:
            return hits[view_cols].sort_index().head(max_hits)

    #match sur colonnes simples
    exact_cols = [c for c in ["gene_symbol", "hgnc_symbol", "entrez_id", "refseq_mrna", "uniprot_swissprot"] if c in var.columns]
    for c in exact_cols:
        hits = var[_as_str(var[c]).str.upper() == q.upper()]
        if not hits.empty:
            return hits[view_cols].sort_index().head(max_hits)

    #alias token exact
    if "alias_symbol" in var.columns:
        mask = _match_alias_token(var["alias_symbol"], q)
        hits = var[mask]
        if not hits.empty:
            return hits[view_cols].sort_index().head(max_hits)

    #fallback contains
    masks = []
    for c in ["gene_symbol", "hgnc_symbol", "base_name", "uniprot_swissprot", "refseq_mrna", "entrez_id"]:
        if c in var.columns:
            masks.append(_as_str(var[c]).str.contains(q, case=False, regex=False))
    if "alias_symbol" in var.columns:
        masks.append(_as_str(var["alias_symbol"]).str.contains(q, case=False, regex=False))

    if masks:
        mask_any = masks[0]
        for m in masks[1:]:
            mask_any = mask_any | m
        hits = var[mask_any]
        if not hits.empty:
            return hits[view_cols].sort_index().head(max_hits)

    #suggestions proches
    candidates = []
    for c in ["gene_symbol", "hgnc_symbol", "base_name"]:
        if c in var.columns:
            candidates += [x for x in _as_str(var[c]).tolist() if x != ""]
    candidates = sorted(set(candidates))

    close = difflib.get_close_matches(q, candidates, n=10, cutoff=0.75)
    if close:
        best = close[0]
        for c in ["gene_symbol", "hgnc_symbol", "base_name"]:
            if c in var.columns:
                hits = var[_as_str(var[c]).str.upper() == best.upper()]
                if not hits.empty:
                    return hits[view_cols].sort_index().head(max_hits)

    return var.iloc[0:0]


def format_hits(hits_df, max_lines=30):
    """
    Cette fonction permet de transformer le DataFrame de résultats 
    en texte lisible pour le terminal.
    """
    if hits_df.empty:
        return "Aucun résultat."

    #on garde des colonnes utiles si présentes
    useful = [c for c in [
        "gene_ids", "gene_symbol", "hgnc_symbol", "alias_symbol",
        "entrez_id", "refseq_mrna", "uniprot_swissprot", "base_name"
    ] if c in hits_df.columns]

    lines = []
    for i, (_, row) in enumerate(hits_df.iterrows(), start=1):
        if i > max_lines:
            lines.append(f"... ({len(hits_df) - max_lines} autres résultats)")
            break

        var_name = row.get("__var_name__", "")
        parts = [f"{i}. {var_name}"]
        for c in useful:
            v = row.get(c, "")
            if pd.isna(v) or str(v).strip() == "":
                continue
            parts.append(f"{c}={v}")
        lines.append(" | ".join(parts))

    return "\n".join(lines)


def choose_var_name(hits_df, choice):
    """
    Quand plusieurs résultats existent, cette fonction permet de récupérer le var_name 
    correspondant au numéro choisi par l’utilisateur   
    """
    if hits_df.empty:
        raise ValueError("Aucun hit.")
    k = int(choice)
    if k < 1 or k > len(hits_df):
        raise ValueError(f"Choix invalide: {k} (attendu 1..{len(hits_df)})")
    return hits_df.iloc[k - 1]["__var_name__"]


def resolve_gene_to_var_name(adata, query, choice=None, max_hits=50):
    """
    Cette fonction permet de résoudre une requête utilisateur en un var_name unique si possible, 
    sinon retourne les résultats possibles.
    """
    hits = search_gene_hits(adata, query, max_hits=max_hits)
    if hits.empty:
        return None, hits
    if len(hits) == 1:
        return hits.iloc[0]["__var_name__"], hits
    if choice is not None:
        return choose_var_name(hits, choice), hits
    return None, hits