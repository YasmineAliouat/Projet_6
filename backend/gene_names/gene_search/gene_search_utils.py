import pandas as pd
import difflib
import re

#Colonnes dans lesquelles on peut chercher un gène.
DEFAULT_SEARCH_COLS = [
    "gene_ids",
    "gene_symbol",
    "hgnc_symbol",
    "hgnc_id",
    "alias_symbol",
    "NCBI_symbol",
    "refseq_mrna",
    "uniprot_swissprot",
    "uniprot_sptrembl",
    "base_name",
]


def ensure_versions_cols(adata):
    """
    Cette fonction ajoute les colonnes base_name et has_versions
    si elles sont absentes.
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

########Fonctionalitées en plus

def _first_matching_token(value: str, query: str):
    """
    Cette fonction permet , si une cellule contient plusieurs valeurs séparées par '|',
    de retourne le premier token qui correspond à la requête
    (exactement ou en contains, insensible à la casse).
    """
    if value is None:
        return None

    text = str(value)
    q = str(query).strip().upper()

    if text == "" or q == "":
        return None

    tokens = text.split("|")

    # priorité au match exact
    for tok in tokens:
        if tok.strip().upper() == q:
            return tok.strip()

    # sinon match partiel
    for tok in tokens:
        if q in tok.strip().upper():
            return tok.strip()

    return None


def _annotate_hits(hits, field_name, query, tokenized=False):
    """
    Cette fonction ajoute aux hits la colonne qui a matché et la valeur exacte qui a matché.
    """
    hits = hits.copy()
    hits["__match_field__"] = field_name

    if tokenized:
        hits["__match_value__"] = hits[field_name].apply(lambda x: _first_matching_token(x, query))
    else:
        hits["__match_value__"] = hits[field_name]

    return hits

def _best_partial_value(row, query):
    """
    Cette fonction retourne la première vraie valeur du row qui contient la requête.
    """
    q = str(query).strip().upper()

    for c in [
        "gene_symbol",
        "hgnc_symbol",
        "hgnc_id",
        "base_name",
        "gene_ids",
        "uniprot_swissprot",
        "uniprot_sptrembl",
        "refseq_mrna",
        "NCBI_symbol",
        "alias_symbol",
    ]:
        if c in row.index:
            value = str(row[c]) if row[c] == row[c] else ""
            if value:
                token = _first_matching_token(value, query)
                if token is not None:
                    return c, token
                if q in value.upper():
                    return c, value
    return None, None

def _is_subsequence(query: str, target: str) -> bool:
    """
    Cette fonction permet une meilleur recherche non exacte des noms de gènes.
    """
    i = 0
    for ch in target:
        if i < len(query) and ch == query[i]:
            i += 1
    return i == len(query)

#########

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
    """
    adata = ensure_versions_cols(adata)

    q = str(query).strip()
    if q == "":
        return pd.DataFrame()

    if search_cols is None:
        search_cols = DEFAULT_SEARCH_COLS

    # Copie minimale , seulement les colonnes utiles, pas tout adata.var
    cols_present = [c for c in search_cols if c in adata.var.columns]
    var = adata.var[cols_present].copy()
    var.insert(0, "__var_name__", adata.var.index.astype(str))
    view_cols = ["__var_name__"] + cols_present

    # ENSG exact
    if "gene_ids" in var.columns and q.startswith("ENSG"):
        hits = var[_as_str(var["gene_ids"]) == q]
        if not hits.empty:
            hits = _annotate_hits(hits, "gene_ids", q, tokenized=False)
            return hits[view_cols + ["__match_field__", "__match_value__"]].head(max_hits)

    # var_names exact
    if q in adata.var_names:
        hits = var.loc[[q]]
        hits = _annotate_hits(hits, "__var_name__", q, tokenized=False)
        return hits[view_cols + ["__match_field__", "__match_value__"]].sort_index().head(max_hits)

    # base_name exact
    if "base_name" in var.columns:
        hits = var[_as_str(var["base_name"]) == q]
        if not hits.empty:
            hits = _annotate_hits(hits, "base_name", q, tokenized=False)
            return hits[view_cols + ["__match_field__", "__match_value__"]].sort_index().head(max_hits)

    # match exact sur colonnes simples
    exact_cols = [c for c in [
        "gene_symbol", "hgnc_symbol", "hgnc_id",
        "NCBI_symbol", "refseq_mrna",
        "uniprot_swissprot", "uniprot_sptrembl"
    ] if c in var.columns]

    for c in exact_cols:
        values = _as_str(var[c]).str.upper()
        q_upper = q.upper()

        hits = var[values == q_upper]

        # cas particulier pour les IDs NCBI stockés comme 54973.0
        if hits.empty and c == "NCBI_symbol":
            hits = var[values.str.replace(".0", "", regex=False) == q_upper]

        if not hits.empty:
            tokenized = c == "uniprot_sptrembl"
            hits = _annotate_hits(hits, c, q, tokenized=tokenized)
            return hits[view_cols + ["__match_field__", "__match_value__"]].sort_index().head(max_hits)

    # alias exact
    if "alias_symbol" in var.columns and len(q) >= 3:
        mask = _match_alias_token(var["alias_symbol"], q)
        hits = var[mask]
        if not hits.empty:
            hits = _annotate_hits(hits, "alias_symbol", q, tokenized=True)
            return hits[view_cols + ["__match_field__", "__match_value__"]].sort_index().head(max_hits)

    # uniprot trembl exact (token)
    if "uniprot_sptrembl" in var.columns and len(q) >= 3:
        mask = _match_alias_token(var["uniprot_sptrembl"], q)
        hits = var[mask]
        if not hits.empty:
            hits = _annotate_hits(hits, "uniprot_sptrembl", q, tokenized=True)
            return hits[view_cols + ["__match_field__", "__match_value__"]].sort_index().head(max_hits)

    return var.iloc[0:0]


def format_hits(hits_df, max_lines=30):
    """
    Cette fonction permet de transformer le DataFrame de résultats 
    en texte lisible pour le terminal.
    """
    if hits_df.empty:
        return "Aucun résultat."

    #on garde des colonnes utiles si présentes
    useful = [c for c in ["gene_ids","gene_symbol","hgnc_symbol","hgnc_id","alias_symbol"
    ,"NCBI_symbol", "refseq_mrna","uniprot_swissprot","uniprot_sptrembl","base_name"
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


def suggest_gene_names(adata, query, n=5):
    """
    Cette fonction permet de proposer des noms de gènes proches si aucun résultat exact n'est trouvé.
    La recherche de suggestions se fait sur les noms principaux les plus utiles.
    """
    adata = ensure_versions_cols(adata)
    var = adata.var  # lecture seule, pas de copie

    candidates = []
    for c in ["gene_symbol", "hgnc_symbol", "base_name", "NCBI_symbol"]:
        if c in var.columns:
            candidates += [x for x in _as_str(var[c]).tolist() if x != ""]

    if "alias_symbol" in var.columns:
        alias_values = _as_str(var["alias_symbol"]).tolist()
        for val in alias_values:
            if val != "":
                candidates += val.split("|")

    # doublons
    candidates = sorted(set(candidates))

    q = str(query).strip().upper()
    if len(q) < 3:
        return []

    # Pré-filtrage : on collecte les candidats par trois critères distincts.
    # le filtre 2 ne tourne que si le 1 est vide.
    # Le filtre 3 (fautes de frappe / substitution) tourne toujours et s'ajoute aux résultats
    pre_set = set()
    substring_hits = [c for c in candidates if q in str(c).upper()]
    if substring_hits:
        pre_set.update(substring_hits)
    else:
        pre_set.update(c for c in candidates if _is_subsequence(q, str(c).upper()))

    # Toujours ajouter les candidats de même longueur ±1 avec ratio suffisant
    pre_set.update(
        c for c in candidates
        if abs(len(str(c)) - len(q)) <= 1
        and difflib.SequenceMatcher(None, q, str(c).upper()).ratio() >= 0.6
    )
    pre_filtered = list(pre_set)

    scored = []
    for cand in pre_filtered:
        cand_up = str(cand).strip().upper()
        if cand_up == "":
            continue

        score = None
        cand_len = len(cand_up)
        q_len = len(q)

        # 1. similarité globale
        ratio = difflib.SequenceMatcher(None, q, cand_up).ratio()
        if ratio >= 0.5:
            score = ratio * 100

        # 2. bonus préfixe
        if cand_up.startswith(q):
            score = (score or 0) + 50

        # 3. bonus si contient
        elif q in cand_up:
            score = (score or 0) + 20

        # 4. bonus si sous-séquence
        elif _is_subsequence(q, cand_up):
            score = (score or 0) + 10

        # 5. proximité de taille
        if score is not None:
            score += max(0, 50 - abs(cand_len - q_len) * 10)
            if cand_up.isalpha() and cand_len <= 6:
                score += 20

        if score is not None:
            scored.append((score, cand))

    scored.sort(key=lambda x: (-x[0], len(x[1])))
    return [cand for _, cand in scored[:n]]

def search_gene_partial(adata, query, max_hits=50):
    """
    Cette fonction permet de rechercher un gène avec un morceau de nom.
    Par exemple : Taper 'TTL' peut retrouver 'TTLL10'.
    """
    adata = ensure_versions_cols(adata)

    q = str(query).strip()
    if len(q) < 2:
        return pd.DataFrame()

    # Copie minimale : seulement les colonnes utiles
    cols_present = [c for c in DEFAULT_SEARCH_COLS if c in adata.var.columns]
    var = adata.var[cols_present].copy()
    var.insert(0, "__var_name__", adata.var.index.astype(str))
    view_cols = ["__var_name__"] + cols_present

    masks = []
    for c in ["gene_symbol", "hgnc_symbol", "hgnc_id", "base_name", "gene_ids",
              "uniprot_swissprot", "uniprot_sptrembl", "refseq_mrna", "NCBI_symbol"]:
        if c in var.columns:
            masks.append(_as_str(var[c]).str.contains(q, case=False, regex=False))

    if "alias_symbol" in var.columns:
        masks.append(_as_str(var["alias_symbol"]).str.contains(q, case=False, regex=False))

    if not masks:
        return pd.DataFrame()

    mask_any = masks[0]
    for m in masks[1:]:
        mask_any = mask_any | m

    hits = var[mask_any].copy()

    # correspondances floues pour les typos.
    # On limite aux noms de longueur len(q)+1 ou +2 pour rester précis et rapide.
    if len(q) >= 3:
        q_up_check = q.upper()
        for c in ["gene_symbol", "hgnc_symbol", "base_name"]:
            if c not in var.columns:
                continue
            col_str = _as_str(var[c])
            len_mask = col_str.str.len().between(len(q_up_check) + 1, len(q_up_check) + 2)
            if not len_mask.any():
                continue
            subset_vals = col_str[len_mask]
            subseq_mask = subset_vals.apply(lambda x: _is_subsequence(q_up_check, x.upper()))
            new_idx = subseq_mask[subseq_mask].index.difference(hits.index)
            if len(new_idx) > 0:
                hits = pd.concat([hits, var.loc[new_idx]])

    # Assignation vectorisée du champ et de la valeur qui ont matché (par priorité de colonne).
    priority_cols = [
        "gene_symbol", "hgnc_symbol", "hgnc_id", "base_name", "gene_ids",
        "uniprot_swissprot", "uniprot_sptrembl", "refseq_mrna", "NCBI_symbol", "alias_symbol",
    ]
    match_field_s = pd.Series([None] * len(hits), index=hits.index, dtype=object)
    match_value_s = pd.Series([""]  * len(hits), index=hits.index, dtype=object)

    for c in priority_cols:
        if c not in hits.columns:
            continue
        still_unset = match_field_s.isna()
        if not still_unset.any():
            break
        col_str = _as_str(hits.loc[still_unset, c])
        matched = col_str.str.contains(q, case=False, regex=False)
        idx = matched[matched].index
        match_field_s.loc[idx] = c
        # extraire le token qui matche
        if c in ("alias_symbol", "uniprot_sptrembl", "refseq_mrna"):
            match_value_s.loc[idx] = col_str.loc[idx].apply(
                lambda x: _first_matching_token(x, q) or str(x)
            )
        else:
            match_value_s.loc[idx] = col_str.loc[idx]

    # Pour les lignes floues (pas de sous-chaîne exacte trouvée), utiliser le symbole principal
    for c in ["gene_symbol", "hgnc_symbol", "base_name"]:
        if c not in hits.columns:
            continue
        remaining = match_field_s.isna()
        if not remaining.any():
            break
        col_vals = _as_str(hits.loc[remaining, c])
        has_value = col_vals != ""
        idx = has_value[has_value].index
        match_field_s.loc[idx] = c
        match_value_s.loc[idx] = col_vals.loc[idx]

    hits["__match_field__"] = match_field_s
    hits["__match_value__"] = match_value_s

    # Score vectorisé : préfixe > contient > longueur proche
    mv_up = hits["__match_value__"].str.upper().fillna("")
    q_up = q.upper()

    hits["__score__"] = (
        mv_up.str.startswith(q_up).astype(float) * 5.0
        + (mv_up == q_up).astype(float) * 3.0
        + (2 - mv_up.str.len().sub(len(q_up)).abs()).clip(lower=0)
    )

    # Bonus de similarité pour les typos
    fuzzy_mask = ~mv_up.str.contains(q_up, regex=False) & (mv_up != "")
    if fuzzy_mask.any():
        hits.loc[fuzzy_mask, "__score__"] += hits.loc[fuzzy_mask, "__match_value__"].apply(
            lambda x: difflib.SequenceMatcher(None, q_up, str(x).upper()).ratio() * 4.0
        )

    hits = hits.sort_values("__score__", ascending=False)
    return hits[view_cols + ["__match_field__", "__match_value__"]].head(max_hits)


def resolve_gene_to_var_name(adata, query, choice=None, max_hits=50):
    """
    Cette fonction permet de résoudre une requête utilisateur en un var_name unique si possible, 
    sinon fait une recherche patielle, puis donne des suggestions proches si rien est trouvé.
    """
    #recherche exacte
    hits = search_gene_hits(adata, query, max_hits=max_hits)

    if not hits.empty:
        if len(hits) == 1:
            return hits.iloc[0]["__var_name__"], hits, [], "exact"

        if choice is not None:
            return choose_var_name(hits, choice), hits, [], "exact"

        return None, hits, [], "exact"

    #recherche partielle
    partial_hits = search_gene_partial(adata, query, max_hits=max_hits)

    if not partial_hits.empty:
        if len(partial_hits) == 1:
            return partial_hits.iloc[0]["__var_name__"], partial_hits, [], "partial"

        return None, partial_hits, [], "partial"

    #suggestions proches
    suggestions = suggest_gene_names(adata, query)

    if suggestions:
        return None, None, suggestions, "suggestion"

    return None, None, [], "none"


