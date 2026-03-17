import importlib
import importlib.util
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import scanpy as sc
import streamlit as st

import backend.gene_coexpression_backend as gcb
import backend.signature_utils as su

from Issue_8.load_anndata import load_anndata, summarize_anndata, validate_anndata

# Classe pour capturer les figures matplotlib, en remplaçant temporairement plt.show pour éviter l'affichage immédiat des figures, et en collectant les nouvelles figures créées pendant la période de capture. Cela permet de contrôler quand les figures sont affichées dans Streamlit, en les rendant disponibles pour un affichage ultérieur via st.pyplot. 
class _PlotCapture:
    def __enter__(self):
        self._old_show = plt.show
        plt.show = lambda *args, **kwargs: None
        self._before = set(plt.get_fignums())
        return self

    def __exit__(self, exc_type, exc, tb):
        plt.show = self._old_show
        self._after = set(plt.get_fignums())

    def figures(self):
        new_nums = sorted(list(self._after - self._before))
        return [plt.figure(n) for n in new_nums]

# Cette fonction charge un fichier AnnData à partir du chemin spécifié, valide sa structure, et génère un résumé de ses caractéristiques. En utilisant le décorateur @st.cache_resource, les résultats de cette fonction sont mis en cache par Streamlit, ce qui permet d'éviter de recharger et de traiter le même fichier plusieurs fois, améliorant ainsi les performances lors de l'exploration interactive des données.
@st.cache_resource
def cached_load_from_path(path: str):
    adata = load_anndata(path)
    validate_anndata(adata)
    rep = summarize_anndata(adata)
    return adata, rep

# Cette fonction charge dynamiquement un module Python à partir d'un chemin de fichier spécifique, en utilisant les fonctionnalités d'importation de Python. Le module chargé contient des utilitaires pour la recherche de gènes, et en le mettant en cache avec @st.cache_resource, on s'assure que le module est chargé une seule fois, même si la fonction est appelée plusieurs fois, ce qui améliore les performances lors de l'exploration interactive des données.
@st.cache_resource
def load_gene_search_utils():
    base_dir = Path(__file__).resolve().parent
    module_path = base_dir / "exploration names" / "gene_search" / "gene_search_utils.py"
    spec = importlib.util.spec_from_file_location("gene_search_utils_dynamic", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def gene_resolution(adata, query: str):
    if not query.strip():
        return None

    gsu = load_gene_search_utils()
    var_name, hits, suggestions, match_type = gsu.resolve_gene_to_var_name(
        adata, query, choice=None, max_hits=30
    )

    st.markdown("**Gene matches:**")

    if hits is not None and len(hits) > 0:
        lines = [f"{idx+1}. {row['__var_name__']}" for idx, (_, row) in enumerate(hits.iterrows())]
        st.code("\n".join(lines), language="text")

    elif suggestions:
        st.warning("No exact match. Suggestions:")
        for s in suggestions:
            st.write(f"- {s}")

    else:
        st.error("No result found.")

    return var_name

# Configuration de la page Streamlit et titre principal de l'application, avec une description de son objectif. La configuration définit le titre de la page et le layout, tandis que le titre et la description expliquent que l'application est une interface web pour explorer des données de scRNA-seq stockées au format AnnData (.h5ad).
st.set_page_config(page_title="scRNA-seq Explorer", layout="wide")
st.title("Single-cell RNA-seq data explorer")
st.markdown("Web interface to explore single-cell RNA-seq data stored in **AnnData (.h5ad)**.")

# Sidebar pour le chargement du dataset, avec des champs pour spécifier le chemin du fichier .h5ad, la source d'expression à utiliser (adata.X ou adata.raw), et un bouton pour lancer le chargement. La sidebar affiche également des informations sur le dataset actuellement chargé, si disponible.
st.sidebar.header("Dataset")
data_path = st.sidebar.text_input("Path to .h5ad file", value="data/adata_3583.h5ad")
expr_source = st.sidebar.selectbox("Expression source", ["adata.X", "adata.raw"], index=0)
load_clicked = st.sidebar.button("Load dataset", type="primary")
st.sidebar.divider()

# Session state initialization pour stocker le dataset chargé, le rapport résumé, et la source du dataset. Cela permet de conserver ces informations entre les interactions de l'utilisateur avec l'application, évitant ainsi de devoir recharger ou recalculer les données à chaque action.
if "dataset" not in st.session_state:
    st.session_state.dataset = None
if "report" not in st.session_state:
    st.session_state.report = None
if "dataset_source" not in st.session_state:
    st.session_state.dataset_source = None

# Load logic : lorsque l'utilisateur clique sur le bouton de chargement, le code vérifie si un chemin de dataset a été fourni. Si c'est le cas, il tente de charger le fichier .h5ad en utilisant la fonction cached_load_from_path, qui est mise en cache pour améliorer les performances. Si l'utilisateur a choisi d'utiliser adata.raw comme source d'expression, le code vérifie si adata.raw est présent et le convertit en AnnData si nécessaire. Ensuite, le dataset chargé, le rapport résumé, et la source du dataset sont stockés dans st.session_state pour une utilisation ultérieure dans l'application. Si une erreur survient lors du chargement, un message d'erreur est affiché dans la sidebar.
if load_clicked:
    try:
        if not data_path:
            st.sidebar.error("Please provide a dataset path.")
        else:
            ds, rep = cached_load_from_path(data_path)

            if expr_source == "adata.raw" and ds.raw is not None:
                ds = ds.raw.to_adata()
            elif expr_source == "adata.raw":
                st.sidebar.error("adata.raw is absent.")

            st.session_state.dataset = ds
            st.session_state.report = rep
            st.session_state.dataset_source = data_path

            st.sidebar.success("Dataset loaded.")

    except Exception as e:
        st.sidebar.error(f"Error: {e}")

if st.session_state.dataset_source:
    st.sidebar.info(f"Active dataset: {st.session_state.dataset_source}")

# Dataset report : un bloc extensible qui affiche un résumé du dataset actuellement chargé, y compris le nombre de cellules, le nombre de gènes, la disponibilité d'une projection UMAP, et le nombre de couches d'expression. Si aucun dataset n'est chargé, des métriques vides sont affichées. Si un dataset est chargé, les métriques sont remplies avec les informations extraites du rapport résumé généré lors du chargement du dataset.
with st.expander("Dataset report", expanded=False):
    rep = st.session_state.report
    adata = st.session_state.dataset
    cols = st.columns(4)

    if rep is None or adata is None:
        cols[0].metric("Cells", "—")
        cols[1].metric("Genes", "—")
        cols[2].metric("UMAP available", "—")
        cols[3].metric("Layers", "—")
    else:
        cols[0].metric("Cells", rep.n_obs)
        cols[1].metric("Genes", rep.n_vars)
        cols[2].metric("UMAP available", "Yes" if "X_umap" in adata.obsm else "No")
        cols[3].metric("Layers", len(rep.layers))

# Tabs for different analyses : trois onglets sont créés pour différentes analyses de l'expression génique. Le premier onglet est dédié à l'expression d'un seul gène, le deuxième à la co-expression de deux gènes, et le troisième au calcul d'un score de signature basé sur une liste de gènes. Chaque onglet contient des champs d'entrée spécifiques pour les paramètres de l'analyse, des options pour sélectionner les types de graphiques à afficher, et un bouton pour lancer l'analyse correspondante.
tab1, tab2, tab3 = st.tabs(["Single gene", "Co-expression (2 genes)", "Signature score"])

# ---------------- TAB 1 ----------------
with tab1:
    st.subheader("Single gene expression")
    adata = st.session_state.dataset

    gene_query = st.text_input("Gene name", placeholder="e.g. MYCN", key="gene_single")

    resolved_gene = gene_resolution(adata, gene_query) if adata is not None and gene_query.strip() else None

    st.markdown("**Plots to display:**")
    c1, c2, c3, c4 = st.columns(4)

    show_umap_clusters = c1.checkbox("UMAP_clusters", value=True, key="sg_show_umap_clusters")
    show_hist = c2.checkbox("Histogram", value=True, key="sg_show_hist")
    show_violin = c3.checkbox("Violin", value=True, key="sg_show_violin")
    show_umap = c4.checkbox("UMAP", value=True, key="sg_show_umap")

    run_single = st.button("Run single gene", key="btn_run_single")

    plot_types = (
        ["violin"] * show_violin +
        ["histogram"] * show_hist +
        ["umap"] * show_umap +
        ["umap_clusters"] * show_umap_clusters
    )

    if run_single:
        if adata is None:
            st.error("Load a dataset first.")
        elif not gene_query.strip():
            st.error("Please enter a gene.")
        elif resolved_gene is None:
            st.error("Invalid gene.")
        elif not plot_types:
            st.warning("Select at least one plot.")
        else:
            try:
                geb = importlib.import_module("backend.gene_expression_backend")
                geb.st = st
                geb.plt = plt
                geb.gene = resolved_gene

                geb.plot_gene_expression(adata, resolved_gene, plot_types=plot_types)

            except Exception as e:
                st.error(f"Error: {e}")

# ---------------- TAB 2 ----------------
with tab2:
    st.subheader("Co-expression (2 genes)")
    adata = st.session_state.dataset

    c1, c2 = st.columns(2)
    with c1:
        gene_a_query = st.text_input("Gene A", key="gene_a")
    with c2:
        gene_b_query = st.text_input("Gene B", key="gene_b")

    resolved_gene_a = gene_resolution(adata, gene_a_query) if adata is not None and gene_a_query.strip() else None
    resolved_gene_b = gene_resolution(adata, gene_b_query) if adata is not None and gene_b_query.strip() else None

    st.markdown("**Plots to display:**")
    p1, p2, p3 = st.columns(3)

    show_umap_coexp = p1.checkbox("UMAP", value=True, key="coexp_show_umap")
    show_scatter = p2.checkbox("Scatter", value=True, key="coexp_show_scatter")
    show_heatmap = p3.checkbox("Heatmap", value=False, key="coexp_show_heatmap")

    run_coexp = st.button("Run co-expression", key="btn_run_coexp")

    plot_types = (
        ["scatter"] * show_scatter +
        ["heatmap"] * show_heatmap +
        ["umap"] * show_umap_coexp
    )

    if run_coexp:
        if adata is None:
            st.error("Load a dataset first.")
        elif not (gene_a_query.strip() and gene_b_query.strip()):
            st.error("Please enter two genes.")
        elif resolved_gene_a is None or resolved_gene_b is None:
            st.error("Please resolve both genes first.")
        elif not plot_types:
            st.warning("Select at least one plot type.")
        else:
            try:
                gcb.st = st
                gcb.plt = plt
                gcb.gene = [resolved_gene_a, resolved_gene_b]

                with _PlotCapture():
                    gcb.plot_gene_expression(adata, plot_types=plot_types)

            except Exception as e:
                st.error(f"Error running co-expression backend: {e}")

# ---------------- TAB 3 ----------------
with tab3:
    st.subheader("Signature score")
    adata = st.session_state.dataset
    st.caption("Compute a gene signature score and visualize it on UMAP.")

    sig_text = st.text_area(
        "Signature genes (one per line)",
        placeholder="MYCN\nPHOX2B\nTH",
        height=140,
        key="sig_text"
    )

    run_sig = st.button("Run signature", key="run_sig")
    # Lorsque l'utilisateur clique sur le bouton pour lancer le calcul du score de signature, le code vérifie que les conditions nécessaires sont remplies (dataset chargé). Si tout est en ordre, il compile la liste des gènes de la signature à partir du champ de texte, vérifie que la liste n'est pas vide, et appelle la fonction plot_signature_score du module signature_utils pour calculer le score de signature et générer une visualisation UMAP correspondante. Si une erreur survient lors de l'exécution du backend, un message d'erreur est affiché.
    if run_sig:
        if adata is None:
            st.error("Load a dataset first.")
        else:
            gene_list = [g.strip() for g in sig_text.split("\n") if g.strip()]

            if len(gene_list) == 0:
                st.error("Please enter at least one gene.")
            else:
                try:
                    su.plot_signature_score(adata, gene_list)
                    st.pyplot(plt.gcf())
                except Exception as e:
                    st.error(f"Error running signature backend: {e}")

st.caption("Prototype – scRNA-seq data exploration tool")