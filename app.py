import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import tempfile
from backend.signature_utils import compute_signature_score

# -------------------------
# Page configuration
# -------------------------
st.set_page_config(
    page_title="scRNA-seq Explorer",
    layout="wide"
)

st.title("Single-cell RNA-seq data explorer")
st.markdown(
    """
    Explore your single-cell RNA-seq dataset interactively.
    """
)

# -------------------------
# Sidebar : Dataset
# -------------------------
st.sidebar.header("Dataset")

data_path = st.sidebar.text_input(
    "Path to AnnData file (.h5ad)",
    value="data/adata_3583.h5ad"
)

# Vérification du fichier
if not os.path.exists(data_path):
    st.error(f"Le fichier AnnData n'existe pas: {data_path}")

# Charger le dataset
adata = None
genes_list = []
if os.path.exists(data_path):
    adata = sc.read_h5ad(data_path)
    genes_list = list(adata.var_names)  # Liste de tous les gènes

# -------------------------
# Expression d'un gène (champ texte avec vérification)
# -------------------------
st.header("Gene expression")

gene_name_input = st.text_input(
    "Tapez le nom d'un gène (insensible à la casse)",
    placeholder="e.g. WASH9P"
)

if gene_name_input and adata is not None:
    # Vérifier que le gène existe (insensible à la casse)
    matches = [g for g in genes_list if g.lower() == gene_name_input.lower()]
    if matches:
        gene_name = matches[0]
        expr = adata[:, gene_name].X
        try:
            expr = expr.toarray()
        except:
            pass
        df = pd.DataFrame(expr, columns=[gene_name])

        st.subheader(f"Expression du gène {gene_name}")
        st.dataframe(df.head(10))  # Affiche les 10 premières lignes

        st.subheader("Distribution de l'expression")
        st.bar_chart(df)

        # UMAP coloré par le gène si disponible
        if "X_umap" in adata.obsm.keys():
            st.subheader("UMAP coloré par le gène")
            fig, ax = plt.subplots(figsize=(6,5))
            sc.pl.umap(adata, color=gene_name, show=False, ax=ax)
            st.pyplot(fig)
    else:
        st.warning(f"Gène '{gene_name_input}' non trouvé dans le dataset.")

# -------------------------
# Signature analysis (multi-gene)
# -------------------------
st.header("Signature analysis (multi-gene)")

if genes_list:
    selected_genes = st.multiselect(
        "Sélectionnez les gènes pour la signature",
        options=genes_list,
        default=genes_list[:5]  # 5 gènes par défaut
    )

    if selected_genes:
        adata = compute_signature_score(adata, selected_genes, "signature_score")
        st.success("Score de signature calculé")

        # Aperçu du score
        st.subheader("Aperçu du score par cellule")
        st.dataframe(adata.obs[["signature_score"]].head())

        # Histogramme
        st.subheader("Distribution du score")
        st.bar_chart(adata.obs["signature_score"])

        # UMAP coloré par le score
        if "X_umap" in adata.obsm.keys():
            st.subheader("UMAP coloré par le score")
            fig, ax = plt.subplots(figsize=(6,5))
            sc.pl.umap(adata, color="signature_score", show=False, ax=ax)
            st.pyplot(fig)

        # Violin plot par type de cellule si colonne existe
        if "cell_type" in adata.obs.columns:
            st.subheader("Violin plot du score par type de cellule")
            fig2, ax2 = plt.subplots(figsize=(8,5))
            sc.pl.violin(adata, keys="signature_score", groupby="cell_type", rotation=90, show=False, ax=ax2)
            st.pyplot(fig2)

        # Télécharger dataset avec le score
        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp_file:
            adata.write_h5ad(tmp_file.name)
            tmp_file_path = tmp_file.name

        st.download_button(
            label="Télécharger dataset avec le score",
            data=open(tmp_file_path, "rb").read(),
            file_name="adata_with_signature.h5ad"
        )

# -------------------------
# Footer
# -------------------------
st.markdown("---")
st.caption("Prototype – scRNA-seq data exploration tool")