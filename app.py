import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import os
import io
from backend.signature_utils import load_signature_genes, compute_signature_score

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
st.sidebar.header("Dataset & Signature")

data_path = st.sidebar.text_input(
    "Path to AnnData file (.h5ad)",
    value="data/adata_3583.h5ad"
)

st.sidebar.markdown(
    """
    Provide valid paths for the dataset and gene signature.
    """
)

signature_path = st.sidebar.text_input(
    "Path to signature file (.txt)",
    value="data/signature.txt"
)

#Vérification de l'existance des fichiers
if not os.path.exists(data_path):
    st.error(f"Le fichier AnnData n'existe pas: {data_path}")
if not os.path.exists(signature_path):
    st.error(f"Le fichier de signature n'existe pas: {signature_path}")

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

#Signature calculation 
if os.path.exists(data_path) and os.path.exists(signature_path):
    st.header("Signature analysis")

    #Charger le dataset et la signature
    adata = sc.read_h5ad(data_path)
    signature = load_signature_genes(signature_path)

    #Vérification des gènes présents et manquants
    genes_present = [g for g in signature if g in adata.var_names]
    genes_missing = [g for g in signature if g not in adata.var_names]

    st.write("Gènes présents:", genes_present)
    st.write("Gènes manquants:", genes_missing)

    if genes_present:
        #Calculer le score
        adata = compute_signature_score(adata, genes_present, "signature_score")
        st.success("Score de signature calculé")

        #Aperçu du score
        st.subheader("Aperçu du score par cellule")
        st.dataframe(adata.obs[["signature_score"]].head())

        #Histogramme rapide
        st.subheader("Distribution du score")
        st.bar_chart(adata.obs["signature_score"])

        #UMAP coloré par le score si UMAP existe
        if "X_umap" in adata.obsm.keys():
            st.subheader("UMAP coloré par le score")
            fig, ax = plt.subplots(figsize=(6,5))
            sc.pl.umap(adata, color="signature_score", show=False, ax=ax)
            st.pyplot(fig)

        #Violin plot par type de cellule si la colonne existe
        if "cell_type" in adata.obs.columns:
            st.subheader("Violin plot du score par type de cellule")
            fig2, ax2 = plt.subplots(figsize=(8,5))
            sc.pl.violin(adata, keys="signature_score", groupby="cell_type", rotation=90, show=False, ax=ax2)
            st.pyplot(fig2)

        #Télécharger le dataset avec le score
        buffer = io.BytesIO()
        adata.write_h5ad(buffer)
        st.dowload_button(
            label="Télécharger dataset avec le score",
            data = buffer.getvalue(),
            file_name="adata_with_signature.h5ad"
        )



# Footer
st.markdown("---")
st.caption("Prototype – scRNA-seq data exploration tool")
