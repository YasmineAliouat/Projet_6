import streamlit as st

st.set_page_config(page_title="scRNA-seq Explorer", layout="wide")
st.title("Single-cell RNA-seq data explorer")

st.markdown(
    """
    Web interface to explore single-cell RNA-seq data stored in **AnnData (.h5ad)**.
    """
)

st.sidebar.header("Dataset")

input_mode = st.sidebar.radio("Dataset input", ["Local path", "Upload file"], index=0)

data_path = None
uploaded_file = None

if input_mode == "Local path":
    data_path = st.sidebar.text_input("Path to .h5ad file", value="data/example.h5ad")
else:
    uploaded_file = st.sidebar.file_uploader("Upload a .h5ad file", type=["h5ad"])

st.sidebar.divider()

# Dataset status (UI only)
if uploaded_file is not None:
    st.sidebar.success("✅ File uploaded (backend not connected yet).")
elif data_path:
    st.sidebar.info("ℹ️ Path set (backend not connected yet).")
else:
    st.sidebar.warning("⚠️ No dataset selected.")

# Backend hooks 
# When we start coding, we will implement these functions in a backend module (like bbackend/data_access.py, backend/plots.py)
def backend_available() -> bool:
    return False  # will become True when backend is plugged

def get_gene_suggestions(_query: str) -> list[str]:
    return ["MYCN", "PHOX2B", "TH", "ALK", "SOX11"]

tab1, tab2, tab3 = st.tabs(["Single gene", "Co-expression", "Signature score"])

with tab1:
    st.subheader("Single gene expression")

    col1, col2, col3 = st.columns([2, 1, 1])
    with col1:
        gene = st.text_input("Gene name", placeholder="e.g. MYCN", key="gene_single")
    with col2:
        scale_mode = st.selectbox("Color scale", ["auto", "linear", "log"], index=0)
    with col3:
        show_stats = st.checkbox("Show stats", value=True)

    if not gene:
        st.info("Enter a gene name to display expression plots.")
    else:
        st.caption(f"Selected gene: **{gene}** | Scale: **{scale_mode}**")

        if not backend_available():
            st.warning("Backend not connected yet. Showing placeholders.")
            suggestions = get_gene_suggestions(gene)
            st.write("Example gene suggestions:", ", ".join(suggestions))

            left, right = st.columns(2)
            with left:
                st.markdown("### Expression distribution")
                st.caption("Histogram / violin plot will appear here.")
                plot_container = st.container()
                with plot_container:
                    st.empty()


            with right:
                st.markdown("### UMAP projection")
                st.caption("UMAP colored by gene expression will appear here.")
                plot_container = st.container()
                with plot_container:
                    st.empty()


            if show_stats:
                st.markdown("### Summary statistics")
                st.table(
                    {
                        "Metric": ["Mean", "% cells > 0", "Min", "Max"],
                        "Value": ["—", "—", "—", "—"],
                    }
                )
        else:
            st.success("Backend connected. (This block will run later.)")

# Co-expression
with tab2:
    st.subheader("Co-expression of two genes")

    c1, c2, c3 = st.columns([2, 2, 1])
    with c1:
        gene_x = st.text_input("Gene X", placeholder="e.g. PHOX2B", key="gene_x")
    with c2:
        gene_y = st.text_input("Gene Y", placeholder="e.g. MYCN", key="gene_y")
    with c3:
        rule = st.selectbox("Rule", ["both > 0", "top quantile"], index=0)

    if not gene_x or not gene_y:
        st.info("Enter two gene names to display co-expression outputs.")
    else:
        st.caption(f"Co-expression: **{gene_x}** + **{gene_y}** | Rule: **{rule}**")

        if not backend_available():
            st.warning("Backend not connected yet. Showing placeholders.")
            left, right = st.columns(2)
            with left:
                st.markdown("### UMAP co-expression map")
                st.caption("Cells co-expressing X and Y will be highlighted here.")
                plot_container = st.container()
                with plot_container:
                    st.empty()


            with right:
                st.markdown("### X vs Y scatter")
                st.caption("Scatter plot of expression values will appear here.")
                plot_container = st.container()
                with plot_container:
                    st.empty()

        else:
            st.success("Backend connected. (This block will run later.)")

# Signature score
with tab3:
    st.subheader("Gene signature score")

    st.markdown("Paste a list of genes (one per line).")
    genes_text = st.text_area(
        "Gene list",
        placeholder="MYCN\nPHOX2B\nTH",
        height=140,
        key="signature_list",
    )

    score_method = st.selectbox("Scoring method", ["mean expression", "score_genes (later)"], index=0)

    genes = [g.strip() for g in genes_text.splitlines() if g.strip()]

    if not genes:
        st.info("Provide at least one gene to compute a signature score.")
    else:
        st.caption(f"Signature genes: **{len(genes)}** | Method: **{score_method}**")

        if not backend_available():
            st.warning("Backend not connected yet. Showing placeholders.")
            st.markdown("### UMAP signature projection")
            st.caption("UMAP colored by signature score will appear here.")
            st.empty()
        else:
            st.success("Backend connected. (This block will run later.)")

with st.expander("About / Notes"):
    st.markdown(
        """
        - This is the **UI skeleton**.  
        - Data loading and plots will be connected by the backend team.  
        - The UI is designed to support: single gene, co-expression, and signatures.
        """
    )

st.markdown("---")
st.caption("Prototype – scRNA-seq data exploration tool")
