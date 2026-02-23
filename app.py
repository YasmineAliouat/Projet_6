import streamlit as st
import scanpy as sc #librairie principale pour lire et manipuler des objets AnnData
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt

#MISE EN PAGE
st.set_page_config(page_title="scRNA-seq Explorer", layout="wide")
st.title("Single-cell RNA-seq data explorer")

st.markdown(
    "Web interface to explore single-cell RNA-seq data stored in **AnnData (.h5ad)**."
)

# SIDEBAR POUR CHOISR LE FICHIER DE DONNÉES (on peut chosir entre local en donnant le chemin ou uploader un fichier)
st.sidebar.header("Dataset")
input_mode = st.sidebar.radio("Dataset input", ["Local path", "Upload file"], index=0)

#Vide au début, à nous de le remplir
data_path = None
uploaded_file = None

# Si local : on affiche un champ texte pour le chemin
# Si upload : on affiche un file_uploader
if input_mode == "Local path":
    data_path = st.sidebar.text_input("Path to .h5ad file", value="data/adata_3583.h5ad")
else:
    uploaded_file = st.sidebar.file_uploader("Upload a .h5ad file", type=["h5ad"])

st.sidebar.divider()
load_clicked = st.sidebar.button("Load dataset", type="primary") #charger au clic

#MÉMORISER LE DATASET CHARGÉ DANS LA SESSION (pour éviter de le recharger à chaque interaction)
if "dataset" not in st.session_state:
    st.session_state.dataset = None
if "report" not in st.session_state:
    st.session_state.report = None
if "dataset_source" not in st.session_state:
    st.session_state.dataset_source = None

# FONCTIONS DE CHARGEMENT (cela permet de cacher le résultat pour éviter de recharger le même fichier plusieurs fois)
@st.cache_resource
def cached_load_dataset_from_path(path: str):
    adata = sc.read_h5ad(path)
    report = {
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "has_umap": ("X_umap" in adata.obsm.keys()),
        "layers": list(adata.layers.keys()),
        "note": None,
    }
    return adata, report


@st.cache_resource
def cached_load_dataset_from_bytes(file_bytes: bytes):
    import tempfile
    with tempfile.NamedTemporaryFile(suffix=".h5ad") as tmp:
        tmp.write(file_bytes)
        tmp.flush()
        adata = sc.read_h5ad(tmp.name)
    report = {
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "has_umap": ("X_umap" in adata.obsm.keys()),
        "layers": list(adata.layers.keys()),
        "note": None,
    }
    return adata, report

# LOGIQUE DU BOUTON DE CHARGEMENT : on essaie de charger le dataset selon le mode choisi, et on stocke le résultat dans la session
if load_clicked:
    try:
        if input_mode == "Upload file":
            if uploaded_file is None:
                st.sidebar.error("Please upload a .h5ad file.")
            else:
                ds, rep = cached_load_dataset_from_bytes(uploaded_file.getvalue())
                st.session_state.dataset = ds
                st.session_state.report = rep
                st.session_state.dataset_source = f"upload:{uploaded_file.name}"
                st.sidebar.success("Dataset successfully loaded.")
        else:
            if not data_path:
                st.sidebar.error("Please provide a dataset path.")
            else:
                ds, rep = cached_load_dataset_from_path(data_path)
                st.session_state.dataset = ds
                st.session_state.report = rep
                st.session_state.dataset_source = data_path
                st.sidebar.success("Dataset successfully loaded.")
    except Exception as e:
        st.sidebar.error(f"Error loading dataset: {e}")

#Affiche le dataset actif
if st.session_state.dataset_source:
    st.sidebar.info(f"Active dataset: {st.session_state.dataset_source}")

# DATASET REPORT (nombre de cellules, de gènes, présence d'UMAP, etc.) dans un expander pour ne pas encombrer l'interface
with st.expander("Dataset report", expanded=False):
    rep = st.session_state.report or {}
    cols = st.columns(4)
    cols[0].metric("Cells", rep.get("n_cells", "—"))
    cols[1].metric("Genes", rep.get("n_genes", "—"))
    cols[2].metric("UMAP available", rep.get("has_umap", "—"))
    cols[3].metric("Layers", len(rep.get("layers", []) or []))
    if rep.get("note"):
        st.caption(rep["note"])

# FONCTIONS UTILES POUR L'EXPLORATION DES GÈNES (retourne la liste des gènes, suggère des gènes similaires, vérifie l'existence d'un gène, récupère le vecteur d'expression d'un gène)
def _var_names(adata):
    return list(adata.var_names.astype(str))


def suggest_genes(adata, query: str, k: int = 10) -> list[str]:
    if adata is None:
        return []
    q = (query or "").strip()
    if not q:
        return []
    q_up = q.upper()
    names = _var_names(adata)

    exact = [g for g in names if g.upper() == q_up]
    if exact:
        return exact[:k]

    starts = [g for g in names if g.upper().startswith(q_up)]
    contains = [g for g in names if q_up in g.upper()]

    if "." in q:
        q_base = q.split(".")[0].upper()
    else:
        q_base = q_up
    base_contains = [g for g in names if g.upper().split(".")[0] == q_base]

    out = []
    for group in (starts, base_contains, contains):
        for g in group:
            if g not in out:
                out.append(g)
            if len(out) >= k:
                return out
    return out[:k]

#Vérifie si le gène existe dans le dataset (en comparant avec les var_names de l'adata)
def gene_exists(adata, gene: str) -> bool:
    if adata is None:
        return False
    return str(gene) in set(adata.var_names.astype(str))

# Récupère le vecteur d'expression d'un gène donné (en gérant le cas sparse et en s'assurant que c'est un vecteur 1D)
def get_gene_vector(adata, gene: str) -> np.ndarray:
    x = adata[:, gene].X
    if sp.issparse(x):
        x = x.toarray()
    x = np.asarray(x).reshape(-1)
    return x

# transformations possibles pour les valeurs d'expression (logarithmique ou linéaire) avant de les passer aux fonctions de tracé
def maybe_transform(values: np.ndarray, scale_mode: str) -> np.ndarray:
    if scale_mode == "log":
        values = np.log1p(np.clip(values, a_min=0, a_max=None))
    return values

# FONCTION PLOT (histogramme, violon, UMAP) : prend les valeurs d'expression et le titre, et retourne une figure matplotlib)
def plot_hist(values: np.ndarray, title: str):
    fig, ax = plt.subplots()
    ax.hist(values, bins=50)
    ax.set_title(title)
    ax.set_xlabel("Expression")
    ax.set_ylabel("Cell count")
    return fig


def plot_violin(values: np.ndarray, title: str):
    fig, ax = plt.subplots()
    ax.violinplot(values, showmeans=True)
    ax.set_title(title)
    ax.set_ylabel("Expression")
    ax.set_xticks([])
    return fig


def plot_umap(adata, color_values: np.ndarray, title: str):
    xy = adata.obsm["X_umap"]
    fig, ax = plt.subplots()
    sca = ax.scatter(xy[:, 0], xy[:, 1], c=color_values, s=6)
    ax.set_title(title)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    fig.colorbar(sca, ax=ax, label="Expression")
    return fig

# CRÉATION DES TABS POUR LES DIFFÉRENTES ANALYSES (expression d'un gène, co-expression de plusieurs gènes, score de signature)
tab1, tab2, tab3 = st.tabs(["Single gene", "Co-expression (multi genes)", "Signature score"])

# TAB 1 : SINGLE GENE EXPRESSION : on affiche un champ pour entrer le nom du gène, des options pour choisir les types de plots à afficher, et on affiche les plots correspondants (histogramme, violon, UMAP) ainsi que des statistiques de base sur l'expression du gène. On gère aussi les cas d'erreur (gène non trouvé, dataset non chargé) et on suggère des gènes similaires si le gène entré n'est pas trouvé. 
with tab1:
    st.subheader("Single gene expression")

    adata = st.session_state.dataset

    top = st.columns([2, 1, 1, 1])
    gene = top[0].text_input("Gene name", placeholder="e.g. MYCN", key="gene_single")
    scale_mode = top[1].selectbox("Color scale", ["auto", "linear", "log"], index=0)
    run_single = top[2].button("Run", key="run_single")
    show_stats = top[3].checkbox("Show stats", value=True)

    st.markdown("**Plots to display:**")
    pcols = st.columns(3)
    show_hist = pcols[0].checkbox("Histogram", value=True)
    show_violin = pcols[1].checkbox("Violin plot", value=False)
    show_umap_opt = pcols[2].checkbox("UMAP", value=True)

    left, right = st.columns(2)
    with left:
        hist_placeholder = st.empty() if show_hist else None
        violin_placeholder = st.empty() if show_violin else None
    with right:
        umap_placeholder = st.empty() if show_umap_opt else None

    if gene and adata is not None and not run_single:
        sugg = suggest_genes(adata, gene)
        if sugg:
            st.caption("Suggestions: " + ", ".join(sugg))

    if run_single:
        if adata is None:
            st.error("Load a dataset first.")
        elif not gene:
            st.error("Please enter a gene name.")
        elif not gene_exists(adata, gene):
            st.error("Gene not found in dataset.")
            sugg = suggest_genes(adata, gene)
            if sugg:
                st.write("Did you mean:", ", ".join(sugg))
        else:
            v = get_gene_vector(adata, gene)
            v_plot = maybe_transform(v, scale_mode)

            if show_hist and hist_placeholder is not None:
                fig = plot_hist(v_plot, f"{gene} — Histogram ({scale_mode})")
                hist_placeholder.pyplot(fig, clear_figure=True)

            if show_violin and violin_placeholder is not None:
                fig = plot_violin(v_plot, f"{gene} — Violin ({scale_mode})")
                violin_placeholder.pyplot(fig, clear_figure=True)

            if show_umap_opt and umap_placeholder is not None:
                if "X_umap" not in adata.obsm.keys():
                    umap_placeholder.warning("UMAP not found in this dataset (obsm['X_umap'] missing).")
                else:
                    fig = plot_umap(adata, v_plot, f"{gene} — UMAP colored by expression")
                    umap_placeholder.pyplot(fig, clear_figure=True)

            if show_stats:
                st.markdown("### Summary statistics")
                stats = {
                    "Mean": float(np.mean(v)),
                    "Min": float(np.min(v)),
                    "Max": float(np.max(v)),
                    "% cells > 0": float(100.0 * np.mean(v > 0)),
                }
                st.table({"Metric": list(stats.keys()), "Value": list(stats.values())})


# TAB 2 : CO-EXPRESSION DE PLUSIEURS GÈNES : on affiche un champ pour entrer une liste de gènes, une option pour choisir la logique (AND/OR) pour définir les cellules co-exprimant ces gènes, des options pour choisir les types de plots à afficher (UMAP avec les cellules co-exprimant les gènes mises en évidence, scatter plot des expressions des 2 gènes si exactement 2 gènes sont fournis), et on affiche les plots correspondants ainsi qu'un résumé du nombre et du pourcentage de cellules co-exprimant les gènes selon la logique choisie. On gère aussi les cas d'erreur (pas de dataset chargé, moins de 2 gènes fournis, certains gènes non trouvés) et on suggère des gènes similaires pour ceux qui ne sont pas trouvés.
with tab2:
    st.subheader("Co-expression (multiple genes)")
    adata = st.session_state.dataset

    st.markdown("Enter a list of genes (one per line), then choose a logical rule to highlight cells.")

    ctop = st.columns([2, 1, 1])
    genes_text = ctop[0].text_area("Genes", placeholder="PHOX2B\nMYCN\nTH", height=130, key="genes_multi")
    logic = ctop[1].selectbox("Logic", ["AND (all genes)", "OR (any gene)"], index=0)
    run_multi = ctop[2].button("Run", key="run_multi")

    genes = [g.strip() for g in (genes_text or "").splitlines() if g.strip()]

    st.markdown("**Plots to display:**")
    pc = st.columns(2)
    show_umap_multi = pc[0].checkbox("UMAP co-expression", value=True)
    show_scatter = pc[1].checkbox("Scatter (2 genes only)", value=True)

    left, right = st.columns(2)
    umap_multi_placeholder = left.empty() if show_umap_multi else None
    scatter_placeholder = right.empty() if show_scatter else None

    if run_multi:
        if adata is None:
            st.error("Load a dataset first.")
        elif len(genes) < 2:
            st.error("Please enter at least 2 genes.")
        else:
            not_found = [g for g in genes if not gene_exists(adata, g)]
            if not_found:
                st.error("Some genes were not found: " + ", ".join(not_found))
                for g in not_found[:5]:
                    sugg = suggest_genes(adata, g)
                    if sugg:
                        st.write(f"Suggestions for {g}:", ", ".join(sugg))
            else:
                exprs = [get_gene_vector(adata, g) for g in genes]
                on = [(e > 0) for e in exprs]
                if logic.startswith("AND"):
                    mask = np.logical_and.reduce(on)
                else:
                    mask = np.logical_or.reduce(on)

                if show_umap_multi and umap_multi_placeholder is not None:
                    if "X_umap" not in adata.obsm.keys():
                        umap_multi_placeholder.warning("UMAP not found in this dataset (obsm['X_umap'] missing).")
                    else:
                        xy = adata.obsm["X_umap"]
                        fig, ax = plt.subplots()
                        ax.scatter(xy[:, 0], xy[:, 1], s=5, alpha=0.4)
                        ax.scatter(xy[mask, 0], xy[mask, 1], s=8)
                        ax.set_title(f"Co-expression ({logic}) — highlighted cells: {int(mask.sum())}")
                        ax.set_xlabel("UMAP1")
                        ax.set_ylabel("UMAP2")
                        umap_multi_placeholder.pyplot(fig, clear_figure=True)

                if show_scatter and scatter_placeholder is not None:
                    if len(genes) != 2:
                        scatter_placeholder.info("Scatter is available only when exactly 2 genes are provided.")
                    else:
                        x = exprs[0]
                        y = exprs[1]
                        fig, ax = plt.subplots()
                        ax.scatter(x, y, s=6, alpha=0.5)
                        ax.set_title(f"{genes[0]} vs {genes[1]}")
                        ax.set_xlabel(genes[0])
                        ax.set_ylabel(genes[1])
                        scatter_placeholder.pyplot(fig, clear_figure=True)

                st.markdown("### Summary")
                st.table(
                    {
                        "Metric": ["Genes", "Logic", "Cells highlighted", "% highlighted"],
                        "Value": [
                            ", ".join(genes),
                            logic,
                            int(mask.sum()),
                            float(100.0 * mask.mean()),
                        ],
                    }
                )

# TAB 3 : GENE SIGNATURE SCORE : on affiche un champ pour entrer une liste de gènes, une option pour choisir la méthode de calcul du score (moyenne, somme, etc.), des options pour choisir les types de plots à afficher (UMAP avec le score de signature comme couleur), et on affiche les plots correspondants ainsi qu'un résumé du nombre de gènes trouvés, du nombre de gènes manquants, et une description de la méthode de calcul du score. On gère aussi les cas d'erreur (pas de dataset chargé, aucun gène fourni, certains gènes non trouvés) et on suggère des gènes similaires pour ceux qui ne sont pas trouvés.
with tab3:
    st.subheader("Gene signature score")
    adata = st.session_state.dataset

    st.markdown("Paste a list of genes (one per line) to compute a per-cell signature score.")

    ctop = st.columns([2, 1, 1])
    sig_text = ctop[0].text_area("Signature genes", placeholder="MYCN\nPHOX2B\nTH", height=130, key="sig_list")
    score_method = ctop[1].selectbox("Scoring method", ["mean"], index=0)
    run_sig = ctop[2].button("Run", key="run_sig")

    sig_genes = [g.strip() for g in (sig_text or "").splitlines() if g.strip()]

    st.markdown("**Plots to display:**")
    show_sig_umap = st.checkbox("UMAP signature", value=True)
    sig_placeholder = st.empty() if show_sig_umap else None

    if run_sig:
        if adata is None:
            st.error("Load a dataset first.")
        elif len(sig_genes) < 1:
            st.error("Please provide at least one gene.")
        else:
            not_found = [g for g in sig_genes if not gene_exists(adata, g)]
            found = [g for g in sig_genes if gene_exists(adata, g)]

            if not found:
                st.error("None of the provided genes were found in the dataset.")
                for g in not_found[:5]:
                    sugg = suggest_genes(adata, g)
                    if sugg:
                        st.write(f"Suggestions for {g}:", ", ".join(sugg))
            else:
                mat = np.vstack([get_gene_vector(adata, g) for g in found])
                sig = mat.mean(axis=0)

                if show_sig_umap and sig_placeholder is not None:
                    if "X_umap" not in adata.obsm.keys():
                        sig_placeholder.warning("UMAP not found in this dataset (obsm['X_umap'] missing).")
                    else:
                        fig = plot_umap(adata, sig, f"Signature score — mean({len(found)} genes)")
                        sig_placeholder.pyplot(fig, clear_figure=True)

                st.markdown("### Summary")
                st.table(
                    {
                        "Metric": ["Genes provided", "Genes found", "Genes missing", "Score (mean)"],
                        "Value": [
                            len(sig_genes),
                            len(found),
                            len(not_found),
                            "mean expression",
                        ],
                    }
                )

                if not_found:
                    st.warning("Missing genes: " + ", ".join(not_found[:20]))
                    if len(not_found) > 20:
                        st.caption("… (truncated)")

# 

st.caption("Prototype – scRNA-seq data exploration tool")