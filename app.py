import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt

from Issue_8.load_anndata import load_anndata, validate_anndata, summarize_anndata
from Issue_12.plot_gene_expression import plot_gene_expression
from Issue_13.plot_gene_coexpression_umap import plot_gene_coexpression_umap


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


def ensure_umap(adata, *, auto_compute: bool):
    if "X_umap" in adata.obsm:
        return
    if not auto_compute:
        return
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)


st.set_page_config(page_title="scRNA-seq Explorer", layout="wide")
st.title("Single-cell RNA-seq data explorer")
st.markdown("Web interface to explore single-cell RNA-seq data stored in **AnnData (.h5ad)**.")


# Sidebar: dataset
st.sidebar.header("Dataset")
data_path = st.sidebar.text_input("Path to .h5ad file", value="data/adata_3583.h5ad")
expr_source = st.sidebar.selectbox("Expression source", ["adata.X", "adata.raw"], index=0)
auto_umap = st.sidebar.checkbox("Compute UMAP if missing", value=True)
load_clicked = st.sidebar.button("Load dataset", type="primary")
st.sidebar.divider()


# Session state
if "dataset" not in st.session_state:
    st.session_state.dataset = None
if "report" not in st.session_state:
    st.session_state.report = None
if "dataset_source" not in st.session_state:
    st.session_state.dataset_source = None

#Cache pour éviter de recharger à chaque fois
@st.cache_resource
def cached_load_from_path(path: str):
    adata = load_anndata(path)
    validate_anndata(adata)
    rep = summarize_anndata(adata)
    return adata, rep


# Load logic
if load_clicked:
    if not data_path:
        st.sidebar.error("Please provide a dataset path.")
    else:
        try:
            ds, rep = cached_load_from_path(data_path)

            # Choix X vs raw
            if expr_source == "adata.raw":
                if ds.raw is None:
                    st.sidebar.error("adata.raw is absent. Choose adata.X.")
                    st.session_state.dataset = None
                    st.session_state.report = None
                    st.session_state.dataset_source = None
                else:
                    ds = ds.raw.to_adata()

            # Vérif UMAP + compute si besoin
            ensure_umap(ds, auto_compute=auto_umap)

            st.session_state.dataset = ds
            st.session_state.report = rep
            st.session_state.dataset_source = data_path
            st.sidebar.success("Dataset successfully loaded.")
        except Exception as e:
            st.sidebar.error(f"Error loading dataset: {e}")


if st.session_state.dataset_source:
    st.sidebar.info(f"Active dataset: {st.session_state.dataset_source}")


# Dataset report
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
        cols[2].metric("UMAP available", ("X_umap" in adata.obsm))
        cols[3].metric("Layers", len(rep.layers))


# Tabs (mise en page en 3 tabs)
tab1, tab2, tab3 = st.tabs(["Single gene", "Co-expression (2 genes)", "Signature score"])

# ---- TAB 1
with tab1:
    st.subheader("Single gene expression")
    adata = st.session_state.dataset

    gene = st.text_input("Gene name", placeholder="e.g. MYCN", key="gene_single")

    st.markdown("**Plots to display:**")
    c1, c2, c3 = st.columns(3)
    show_hist = c1.checkbox("Histogram", value=True, key="sg_show_hist")
    show_violin = c2.checkbox("Violin plot", value=True, key="sg_show_violin")
    show_umap_plot = c3.checkbox("UMAP", value=True, key="sg_show_umap")

    run = st.button("Run single gene", key="btn_run_single")

    if run:
        if adata is None:
            st.error("Load a dataset first.")
        elif not gene:
            st.error("Please enter a gene name.")
        elif gene not in adata.var_names:
            st.error("Gene not found in dataset.")
        else:
            plot_types = []
            if show_violin:
                plot_types.append("violin")
            if show_hist:
                plot_types.append("histogram")
            if show_umap_plot:
                plot_types.append("umap")

            if not plot_types:
                st.warning("Select at least one plot type.")
            else:
                import backend.gene_expression_backend as geb

                geb.st = st
                geb.plt = plt
                geb.gene = gene  

                import scanpy as sc
                import numpy as np
                import scipy.sparse as sp
                import matplotlib.pyplot as plt

                def _sc_pl_histogram(adata, color=None, show=False, ax=None, **kwargs):
                    genes = color if isinstance(color, (list, tuple)) else [color]
                    if ax is None:
                        _, ax = plt.subplots()

                    for g in genes:
                        x = adata[:, g].X
                        if sp.issparse(x):
                            x = x.toarray()
                        x = np.asarray(x).reshape(-1)
                        ax.hist(x, bins=50, alpha=0.6, label=str(g))

                    ax.set_title("Histogram")
                    ax.set_xlabel("Expression")
                    ax.set_ylabel("Cell count")
                    if len(genes) > 1:
                        ax.legend()
                    return ax

                if not hasattr(sc.pl, "histogram"):
                    sc.pl.histogram = _sc_pl_histogram

                try:
                    geb.plot_gene_expression(adata, plot_types=plot_types)
                except Exception as e:
                    st.error(f"Error running gene expression backend: {e}")

# TAB 2
with tab2:
    st.subheader("Co-expression (2 genes)")
    adata = st.session_state.dataset

    c1, c2 = st.columns(2)
    with c1:
        gene_a = st.text_input("Gene A", key="gene_a")
    with c2:
        gene_b = st.text_input("Gene B", key="gene_b")

    st.markdown("**Plots to display:**")
    p1, p2, p3 = st.columns(3)
    show_umap_coexp = p1.checkbox("UMAP", value=True, key="coexp_show_umap")
    show_scatter = p2.checkbox("Scatter", value=True, key="coexp_show_scatter")
    show_heatmap = p3.checkbox("Heatmap", value=False, key="coexp_show_heatmap")

    run = st.button("Run co-expression")

    if run:
        if adata is None:
            st.error("Load a dataset first.")
        elif not gene_a or not gene_b:
            st.error("Please enter two genes.")
        elif gene_a not in adata.var_names or gene_b not in adata.var_names:
            st.error("One gene (or both) not found.")
        else:
            with _PlotCapture() as cap:
                plot_gene_coexpression_umap(adata, gene_a, gene_b)

            figs = cap.figures()

            idx = 0
            if show_umap_coexp and idx < len(figs):
                st.pyplot(figs[idx])
            idx += 1

            if show_scatter and idx < len(figs):
                st.pyplot(figs[idx])
            idx += 1

            if show_heatmap and idx < len(figs):
                st.pyplot(figs[idx])


# TAB 3
with tab3:
    st.subheader("Signature score")
    adata = st.session_state.dataset
    st.caption("Cet onglet doit appeler le module de scoring de signature de tes camarades.")

    sig_text = st.text_area(
        "Signature genes (one per line)",
        placeholder="MYCN\nPHOX2B\nTH",
        height=140,
        key="sig_text"
    )
    run = st.button("Run signature", key="run_sig")

    if run:
        if adata is None:
            st.error("Load a dataset first.")
        else:
            st.warning(
                "Missing"
            )

st.caption("Prototype – scRNA-seq data exploration tool")