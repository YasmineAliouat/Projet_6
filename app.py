import streamlit as st

# Page configuration

st.set_page_config(
    page_title="scRNA-seq Explorer",
    layout="wide"
)

st.title("Single-cell RNA-seq data explorer")

st.markdown(
    """
    This application allows interactive exploration of single-cell RNA-seq data
    stored in an AnnData object.
    """
)

# Sidebar
st.sidebar.header("Dataset")

data_path = st.sidebar.text_input(
    "Path to AnnData file (.h5ad)",
    value="data/example.h5ad"
)

st.sidebar.markdown(
    """
    The dataset is expected to be provided in AnnData format.
    """
)

# Gene selection
st.header("Gene expression")

gene_name = st.text_input(
    "Enter a gene name",
    placeholder="e.g. MYCN"
)

if gene_name:
    st.info(
        f"Gene selected: {gene_name}\n\n"
        "Expression values and UMAP projection will be displayed here."
    )
else:
    st.warning("Please enter a gene name to visualize its expression.")

# UMAP 
st.header("UMAP projection")

st.markdown(
    """
    A 2D projection (e.g. UMAP) colored by gene expression
    will be displayed here once the data is loaded.
    """
)

st.empty()

# Advanced analysis
with st.expander("Advanced analysis (coming soon)"):
    st.markdown(
        """
        - Co-expression analysis between two genes  
        - Correlation with pseudotime  
        - Filtering by cell annotations  
        """
    )

# Footer
st.markdown("---")
st.caption("Prototype – scRNA-seq data exploration tool")
