import importlib
from exploration_names.integration_interface.integration_interface import (
    resolve_gene_for_streamlit,
    hits_to_options,
    suggestions_to_options
)
from exploration_names.prepare_names.prepare_names import gene_versions
from typing import Optional
import matplotlib.pyplot as plt
import scanpy as sc
import streamlit as st
import backend.gene_coexpression_backend as gcb
import backend.signature_utils as su
import backend.heatmap_backend as hb
from Issue_8.load_anndata import load_anndata, summarize_anndata, validate_anndata
import io
import zipfile

# Class to capture matplotlib figures, temporarily overriding plt.show to prevent immediate display of figures, and collecting new figures created during the capture period. This allows control over when figures are displayed in Streamlit, making them available for later display via st.pyplot.
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

#This function loads an AnnData file from the specified path, validates its structure, and generates a summary of its characteristics. By using the @st.cache_resource decorator, the results of this function are cached by Streamlit, which avoids reloading and processing the same file multiple times, thus improving performance during interactive data exploration.
@st.cache_resource
def cached_load_from_path(path: str):
    adata = load_anndata(path)
    validate_anndata(adata)
    adata = gene_versions(adata)
    adata.var["base_name"] = adata.var["base_name"].astype(str)
    rep = summarize_anndata(adata)
    return adata, rep


# Streamlit page configuration and main application title, with a description of its purpose. The configuration defines the page title and layout, while the title and description explain that the application is a web interface for exploring scRNA-seq data stored in AnnData (.h5ad) format.
st.set_page_config(page_title="scRNA-seq Explorer", layout="wide")
st.title("Single-cell RNA-seq data explorer")
st.markdown("Web interface to explore single-cell RNA-seq data stored in **AnnData (.h5ad)**.")

# Sidebar for loading the dataset, with fields to specify the path to the .h5ad file, the expression source to use (adata.X or adata.raw), and a button to start loading. The sidebar also displays information about the currently loaded dataset, if available.
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

# Load logic: When the user clicks the load button, the code checks if a dataset path has been provided. If so, it attempts to load the .h5ad file using the cached_load_from_path function, which is cached to improve performance. If the user has chosen to use adata.raw as the expression source, the code checks if adata.raw is present and converts it to AnnData if necessary. The loaded dataset, summary report, and dataset source are then stored in st.session_state for later use in the application. If an error occurs during loading, an error message is displayed in the sidebar.
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

# Dataset report: an expandable block that displays a summary of the currently loaded dataset, including the number of cells, the number of genes, the availability of a UMAP projection, and the number of expression layers. If no dataset is loaded, empty metrics are displayed. If a dataset is loaded, the metrics are populated with information extracted from the summary report generated when the dataset was loaded.
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

# Tabs for different analyses: Three tabs are created for different gene expression analyses. The first tab is dedicated to the expression of a single gene, the second to the co-expression of two genes, and the third to the calculation of a signature score based on a list of genes. Each tab contains specific input fields for the analysis parameters, options to select the types of graphs to display, and a button to launch the corresponding analysis.
tab1, tab2, tab3 = st.tabs(["Single gene", "Co-expression (2 genes)", "Co-expression (multiple genes)"])

# ---------------- TAB 1 ----------------
with tab1:
    st.subheader("Single gene expression")
    adata = st.session_state.dataset

    gene_query = st.text_input("Gene name", placeholder="e.g. MYCN", key="gene_single")

    resolved_gene = None

    if adata is not None and gene_query.strip():
        result = resolve_gene_for_streamlit(adata, gene_query)

        if result["status"] == "exact":
            resolved_gene = result["var_name"]

        elif result["status"] == "multiple":
            options = hits_to_options(result["hits"])
            labels = [x[0] for x in options]
            label_to_var = {x[0]: x[1] for x in options}

            selected = st.selectbox("Select a gene", labels, key="sg_select")
            resolved_gene = label_to_var[selected]

        elif result["status"] == "suggestions":
            options = suggestions_to_options(adata, result["suggestions"])
            labels = [x[0] for x in options]
            label_to_var = {x[0]: x[1] for x in options}

            selected = st.selectbox("Suggestions", labels, key="sg_suggest")
            resolved_gene = label_to_var[selected]

        else:
            st.error("No result found.")
    
    st.markdown("**Plots to display:**")
    c1, c2, c3, c4 = st.columns(4)

    show_umap_clusters = c1.checkbox("UMAP_clusters", value=True, key="sg_show_umap_clusters")
    show_umap = c2.checkbox("UMAP_gene", value=True, key="sg_show_umap")
    show_violin = c3.checkbox("Violin_plot", value=True, key="sg_show_violin")    
    show_hist = c4.checkbox("Histogram", value=True, key="sg_show_hist")

    run_single = st.button("Run single gene", key="btn_run_single")

    plot_types = (
        ["umap_clusters"] * show_umap_clusters +
        ["umap"] * show_umap +
        ["violin"] * show_violin +
        ["histogram"] * show_hist
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

                figs = geb.plot_gene_expression(adata, resolved_gene, plot_types=plot_types)

                if figs:
                    zip_buffer = io.BytesIO()

                    with zipfile.ZipFile(zip_buffer, "w") as zf:
                        for i, fig in enumerate(figs):
                            img = io.BytesIO()
                            fig.savefig(img, format="png", bbox_inches="tight")
                            img.seek(0)
                            zf.writestr(f"plot_{i+1}.png", img.read())

                    zip_buffer.seek(0)

                    st.download_button(
                        "Download all plots",
                        zip_buffer,
                        file_name="plots.zip",
                        mime="application/zip"
                    )

            except Exception as e:
                st.error(f"Error: {e}")

# ---------------- TAB 2 ----------------
with tab2:
    st.subheader("Co-expression (2 genes)")
    adata = st.session_state.dataset

    c1, c2 = st.columns(2)
    with c1:
        gene_a = st.text_input("Gene A")
    with c2:
        gene_b = st.text_input("Gene B")

    genes_list = [gene_a, gene_b]

    resolved_genes = []

    if adata:
        for gene, key in [(gene_a, "A"), (gene_b, "B")]:
            if not gene.strip():
                resolved_genes.append(None)
                continue

            result = resolve_gene_for_streamlit(adata, gene)

            if result["status"] == "exact":
                resolved_genes.append(result["var_name"])

            elif result["status"] == "multiple":
                options = hits_to_options(result["hits"])
                labels = [x[0] for x in options]
                label_to_var = {x[0]: x[1] for x in options}

                selected = st.selectbox(f"Select gene {key}", labels, key=f"tab2_select_{key}")
                resolved_genes.append(label_to_var[selected])

            elif result["status"] == "suggestions":
                options = suggestions_to_options(adata, result["suggestions"])
                labels = [x[0] for x in options]
                label_to_var = {x[0]: x[1] for x in options}

                selected = st.selectbox(f"Suggestions gene {key}", labels, key=f"tab2_suggest_{key}")
                resolved_genes.append(label_to_var[selected])

            else:
                st.error(f"No result for gene {key}")
                resolved_genes.append(None)

    st.markdown("**Plots to display:**")
    p1, p2 = st.columns(2)

    show_umap = p1.checkbox("UMAP", True)
    show_scatter = p2.checkbox("Scatter_plot", True)

    run_coexp = st.button("Run co-expression")

    plot_types = (
        ["scatter"] * show_scatter +
        ["umap"] * show_umap
    )

    if run_coexp:
        if adata is None:
            st.error("Load a dataset first.")
        elif not (gene_a.strip() and gene_b.strip()):
            st.error("Enter 2 genes.")
        elif any(g is None for g in resolved_genes):
            st.error("Invalid gene(s).")
        elif not plot_types:
            st.warning("Select at least one plot.")
        else:
            try:
                gcb.st = st
                gcb.plt = plt

                figs = gcb.plot_gene_coexpression(
                    adata,
                    genes=resolved_genes,
                    plot_types=plot_types
                )

                if figs:
                    zip_buffer = io.BytesIO()

                    with zipfile.ZipFile(zip_buffer, "w") as zf:
                        for i, fig in enumerate(figs):

                            img = io.BytesIO()

                            if hasattr(fig, "write_image"):
                                fig.write_image(img, format="png")
                            else:
                                fig.savefig(img, format="png", bbox_inches="tight")

                            img.seek(0)
                            zf.writestr(f"plot_{i+1}.png", img.read())

                    zip_buffer.seek(0)

                    st.download_button(
                        "Download all plots",
                        zip_buffer,
                        file_name="coexpression_plots.zip",
                        mime="application/zip"
                    )

            except Exception as e:
                st.error(f"Error: {e}")

# ---------------- TAB 3 ----------------
with tab3:
    st.subheader("Co-expression (multiple genes)")
    adata = st.session_state.dataset

    sig_text = st.text_area(
        "Genes (one per line)",
        placeholder="MYCN\nPHOX2B\nTH",
        key="sig_text"
    )

    st.markdown("**Plots to display:**")
    c1, c2 = st.columns(2)

    show_umap = c1.checkbox("UMAP_signature", True, key="tab3_umap")
    show_heatmap = c2.checkbox("Heatmap", True, key="tab3_heatmap")

    if "tab3_run" not in st.session_state:
        st.session_state.tab3_run = False

    if st.button("Run analysis", key="run_sig"):
        st.session_state.tab3_run = True

    plot_types = (
        ["umap"] * show_umap +
        ["heatmap"] * show_heatmap
    )

    if st.session_state.tab3_run:
        if adata is None:
            st.error("Load a dataset first.")
        else:
            raw_genes = [g.strip() for g in sig_text.split("\n") if g.strip()]

            resolved_genes = []
            invalid_genes = []
            selection_needed = False

            for i, g in enumerate(raw_genes):
                result = resolve_gene_for_streamlit(adata, g)

                if result["status"] == "exact":
                    resolved_genes.append(result["var_name"])

                elif result["status"] == "multiple":
                    options = hits_to_options(result["hits"])
                    labels = [x[0] for x in options]
                    label_to_var = {x[0]: x[1] for x in options}

                    selected = st.selectbox(
                        f"Select match for '{g}'",
                        labels,
                        key=f"tab3_select_{i}"
                    )

                    resolved_genes.append(label_to_var[selected])
                    selection_needed = True

                elif result["status"] == "suggestions":
                    options = suggestions_to_options(adata, result["suggestions"])
                    labels = [x[0] for x in options]
                    label_to_var = {x[0]: x[1] for x in options}

                    selected = st.selectbox(
                        f"Suggestions for '{g}'",
                        labels,
                        key=f"tab3_suggest_{i}"
                    )

                    resolved_genes.append(label_to_var[selected])
                    selection_needed = True

                else:
                    invalid_genes.append(g)

            if invalid_genes:
                st.warning(f"Ignored (invalid): {', '.join(invalid_genes)}")

            if selection_needed:
                st.info("Please validate all gene selections to run the analysis.")
                st.stop()

            if len(resolved_genes) == 0:
                st.error("No valid genes found.")
            elif not plot_types:
                st.warning("Select at least one plot.")
            else:
                try:
                    resolved_genes = [g for g in resolved_genes if g in adata.var_names]

                    if len(resolved_genes) == 0:
                        st.error("None of the input genes are present in the dataset.")
                    else:
                        st.write(f"{len(resolved_genes)} valid genes used")

                        all_figs = []

                        if show_umap:
                            st.subheader("UMAP signature score")

                            with _PlotCapture() as cap_sig:
                                su.plot_signature_score(adata, resolved_genes)

                            figs = cap_sig.figures()
                            for fig in figs:
                                st.pyplot(fig)
                            all_figs.extend(figs)

                        if show_heatmap:
                            st.subheader("Co-expression heatmap")

                            with _PlotCapture() as cap_heat:
                                hb.coexp_heatmap(adata, resolved_genes)

                            figs = cap_heat.figures()
                            for fig in figs:
                                st.pyplot(fig)
                            all_figs.extend(figs)

                        if all_figs:
                            zip_buffer = io.BytesIO()

                            with zipfile.ZipFile(zip_buffer, "w") as zf:
                                for i, fig in enumerate(all_figs):
                                    img = io.BytesIO()
                                    fig.savefig(img, format="png", bbox_inches="tight")
                                    img.seek(0)
                                    zf.writestr(f"tab3_plot_{i+1}.png", img.read())

                            zip_buffer.seek(0)

                            st.download_button(
                                "Download all plots",
                                zip_buffer,
                                file_name="tab3_plots.zip",
                                mime="application/zip"
                            )

                except Exception as e:
                    st.error(f"Error running analysis: {e}")

st.caption("Prototype – scRNA-seq data exploration tool")