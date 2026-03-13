import importlib
import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc
import streamlit as st

from Issue_8.load_anndata import load_anndata, summarize_anndata, validate_anndata


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


@st.cache_resource
def cached_load_from_path(path: str):
    adata = load_anndata(path)
    validate_anndata(adata)
    rep = summarize_anndata(adata)
    return adata, rep


@st.cache_resource
def load_gene_search_utils():
    base_dir = Path(__file__).resolve().parent
    module_path = base_dir / "exploration names" / "gene_search" / "gene_search_utils.py"

    spec = importlib.util.spec_from_file_location("gene_search_utils_dynamic", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def patch_histogram_if_missing():
    import numpy as np
    import scipy.sparse as sp

    if hasattr(sc.pl, "histogram"):
        return

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

    sc.pl.histogram = _sc_pl_histogram


def resolve_gene_ui(adata, query: str, choice: int | None = None, max_hits: int = 30):
    gsu = load_gene_search_utils()
    return gsu.resolve_gene_to_var_name(adata, query, choice=choice, max_hits=max_hits)


def render_gene_resolution_block(title: str, adata, query: str, state_prefix: str):
    gsu = load_gene_search_utils()

    if not query.strip():
        return None

    choice_key = f"{state_prefix}_choice"

    if choice_key not in st.session_state:
        st.session_state[choice_key] = None

    var_name, hits, suggestions, match_type = resolve_gene_ui(
        adata,
        query,
        choice=st.session_state[choice_key],
        max_hits=30,
    )

    with st.container():
        st.markdown(f"**{title}**")

        if var_name is not None:
            if match_type == "exact":
                st.success(f"Resolved gene: `{var_name}` (exact match)")
            elif match_type == "partial":
                st.success(f"Resolved gene: `{var_name}` (partial match)")
            return var_name

        if hits is not None and len(hits) > 0:
            if match_type == "exact":
                st.info("Several exact matches were found.")
            elif match_type == "partial":
                st.info("Several partial matches were found.")

            display_lines = []
            for i, (_, row) in enumerate(hits.iterrows(), start=1):
                var_label = row.get("__var_name__", "")
                symbol = row.get("gene_symbol", "")
                hgnc = row.get("hgnc_symbol", "")
                base_name = row.get("base_name", "")

                parts = [f"{i}. {var_label}"]
                if symbol:
                    parts.append(f"gene_symbol={symbol}")
                if hgnc:
                    parts.append(f"hgnc_symbol={hgnc}")
                if base_name:
                    parts.append(f"base_name={base_name}")

                display_lines.append(" | ".join(parts))

            st.code("\n".join(display_lines), language="text")

            option_labels = [f"{i}. {hits.iloc[i-1]['__var_name__']}" for i in range(1, len(hits) + 1)]
            selected_label = st.selectbox(
                "Choose one result",
                options=option_labels,
                key=f"{state_prefix}_selectbox",
            )

            selected_idx = int(selected_label.split(".", 1)[0])
            if st.button("Confirm gene choice", key=f"{state_prefix}_confirm"):
                st.session_state[choice_key] = selected_idx
                st.rerun()

            return None

        if suggestions:
            st.warning("No exact result found.")
            st.write("Suggestions:")
            for s in suggestions:
                st.write(f"- `{s}`")
            return None

        st.error("No result found.")
        return None


st.set_page_config(page_title="scRNA-seq Explorer", layout="wide")
st.title("Single-cell RNA-seq data explorer")
st.markdown("Web interface to explore single-cell RNA-seq data stored in **AnnData (.h5ad)**.")

# Sidebar
st.sidebar.header("Dataset")
data_path = st.sidebar.text_input("Path to .h5ad file", value="data/adata_3583.h5ad")
expr_source = st.sidebar.selectbox("Expression source", ["adata.X", "adata.raw"], index=0)
load_clicked = st.sidebar.button("Load dataset", type="primary")
st.sidebar.divider()

# Session state
if "dataset" not in st.session_state:
    st.session_state.dataset = None
if "report" not in st.session_state:
    st.session_state.report = None
if "dataset_source" not in st.session_state:
    st.session_state.dataset_source = None

# Load logic
if load_clicked:
    if not data_path:
        st.sidebar.error("Please provide a dataset path.")
    else:
        try:
            ds, rep = cached_load_from_path(data_path)

            if expr_source == "adata.raw":
                if ds.raw is None:
                    st.sidebar.error("adata.raw is absent. Choose adata.X.")
                    st.session_state.dataset = None
                    st.session_state.report = None
                    st.session_state.dataset_source = None
                else:
                    ds = ds.raw.to_adata()

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
        cols[2].metric("UMAP available", "Yes" if "X_umap" in adata.obsm else "No")
        cols[3].metric("Layers", len(rep.layers))

# Tabs
tab1, tab2, tab3 = st.tabs(["Single gene", "Co-expression (2 genes)", "Signature score"])

# ---------------- TAB 1 ----------------
with tab1:
    st.subheader("Single gene expression")
    adata = st.session_state.dataset

    gene_query = st.text_input("Gene name", placeholder="e.g. MYCN", key="gene_single")

    resolved_gene = None
    if adata is not None and gene_query.strip():
        resolved_gene = render_gene_resolution_block(
            "Gene resolution",
            adata,
            gene_query,
            "single_gene",
        )

    st.markdown("**Plots to display:**")
    c1, c2, c3 = st.columns(3)
    show_hist = c1.checkbox("Histogram", value=True, key="sg_show_hist")
    show_violin = c2.checkbox("Violin plot", value=True, key="sg_show_violin")
    show_umap_plot = c3.checkbox("UMAP", value=True, key="sg_show_umap")

    run_single = st.button("Run single gene", key="btn_run_single")

    if run_single:
        if adata is None:
            st.error("Load a dataset first.")
        elif not gene_query.strip():
            st.error("Please enter a gene name.")
        elif resolved_gene is None:
            st.error("Please resolve the gene first.")
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
                try:
                    patch_histogram_if_missing()
                    geb = importlib.import_module("backend.gene_expression_backend")
                    geb.st = st
                    geb.plt = plt
                    geb.gene = resolved_gene
                    geb.plot_gene_expression(adata, plot_types=plot_types)
                except Exception as e:
                    st.error(f"Error running gene expression backend: {e}")

# ---------------- TAB 2 ----------------
with tab2:
    st.subheader("Co-expression (2 genes)")
    adata = st.session_state.dataset

    c1, c2 = st.columns(2)
    with c1:
        gene_a_query = st.text_input("Gene A", key="gene_a")
    with c2:
        gene_b_query = st.text_input("Gene B", key="gene_b")

    resolved_gene_a = None
    resolved_gene_b = None

    if adata is not None and gene_a_query.strip():
        resolved_gene_a = render_gene_resolution_block(
            "Gene A resolution",
            adata,
            gene_a_query,
            "coexp_gene_a",
        )

    if adata is not None and gene_b_query.strip():
        resolved_gene_b = render_gene_resolution_block(
            "Gene B resolution",
            adata,
            gene_b_query,
            "coexp_gene_b",
        )

    st.markdown("**Plots to display:**")
    p1, p2, p3 = st.columns(3)
    show_umap_coexp = p1.checkbox("UMAP", value=True, key="coexp_show_umap")
    show_scatter = p2.checkbox("Scatter", value=True, key="coexp_show_scatter")
    show_heatmap = p3.checkbox("Heatmap", value=False, key="coexp_show_heatmap")

    run_coexp = st.button("Run co-expression", key="btn_run_coexp")

    if run_coexp:
        if adata is None:
            st.error("Load a dataset first.")
        elif not gene_a_query.strip() or not gene_b_query.strip():
            st.error("Please enter two genes.")
        elif resolved_gene_a is None or resolved_gene_b is None:
            st.error("Please resolve both genes first.")
        else:
            plot_types = []
            if show_scatter:
                plot_types.append("scatter")
            if show_heatmap:
                plot_types.append("heatmap")
            if show_umap_coexp:
                plot_types.append("umap_coexpression")

            if not plot_types:
                st.warning("Select at least one plot type.")
            else:
                try:
                    gcb = importlib.import_module("backend.gene_coexpression_backend")
                    gcb.st = st
                    gcb.plt = plt
                    gcb.gene = [resolved_gene_a, resolved_gene_b]

                    with _PlotCapture():
                        gcb.plot_gene_expression(adata, plot_types=plot_types)

                except Exception as e:
                    st.error(f"Error running co-expression backend: {e}")
                    st.info(
                        "The frontend calls the backend only. "
                        "If this fails, the issue is currently in backend/gene_coexpression_backend.py."
                    )

# ---------------- TAB 3 ----------------
with tab3:
    st.subheader("Signature score")
    adata = st.session_state.dataset
    st.caption("Compute a gene signature score and visualize it on UMAP.")

    import backend.signature_utils as su

    sig_text = st.text_area(
        "Signature genes (one per line)",
        placeholder="MYCN\nPHOX2B\nTH",
        height=140,
        key="sig_text"
    )

    run_sig = st.button("Run signature", key="run_sig")

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