import argparse
import scanpy as sc
import pandas as pd
import os
import tempfile
from scanpy.queries import biomart_annotations

def gene_versions(adata):
    """
    Cette fonction rajoute les colonne base_name et has_versions au adata. 
    Si un gène a plus qu'une version on aura True dans has_versions,
    et le nom sans sufix (.X) dans base_name
    """
    names= pd.Index(adata.var_names.astype(str))
    adata.var["base_name"]=names.str.replace(r"\.\d+$", "", regex=True) #garde pas le (.X)
    adata.var["has_versions"]=names.str.contains(r"\.\d+$", regex=True) # a plus qu'une version ou pas
    return adata

def make_groups(annot: pd.DataFrame)-> pd.DataFrame:
    """
    Cette fonction permet de regrouper les veleurs différente pour un 
    même gene_id.
    """
    return(annot.groupby("ensembl_gene_id", sort=False)
    .agg(lambda x: "|". join(x.dropna().astype(str).unique()))
    ) #Les noms seront séparer par un '|'

def mapping(adata, data, grouped_names, mapping_dict):
    """
    Cette fonction permet de faire la correspondance entre le nom de la 
    colonne dans les AnnData et le nom présent dans BioMart
    """
    for biomart_col, out_col in mapping_dict.items():
        if biomart_col in grouped_names: #Car BioMart ne renvoie pas toujours toutes les colonnes
            adata.var[out_col]=data.map(grouped_names[biomart_col]) 
    return adata

def biomart_query_with_fallback(organism, cols):
    """
    Cette fonction interroge BioMart avec un ou plusieurs hosts de secours.
    """
    hosts = [
        "http://www.ensembl.org",
        "http://useast.ensembl.org",
        "http://asia.ensembl.org",
    ]

    last_error = None

    for host in hosts:
        try:
            return biomart_annotations(organism, cols, host=host)
        except Exception as e:
            last_error = e

    raise last_error

def add_biomart_names(adata, organism="hsapiens"):
    """
    Cette fonction ajoute plusierurs colonnes si disponible via BioMart:
    gene_symbol; hgnc_symbol; hgnc_id, alias_symbol, entrezgene_id (NCBI)
    refseq_mrna, refseq_peptide, uniprotswissprot, uniprotsptrembl
    Séparer en plusieurs requetes car BioMart permet la rechrche de 4 attribus à la fois.
    """
    if "gene_ids" not in adata.var.columns:
        print("gene_ids manquant : utilisation de adata.var_names comme fallback")
        adata.var["gene_ids"] = adata.var_names.astype(str)  #Si on utilise un autre fichier

    ids= adata.var["gene_ids"].astype(str)

    # Requete 1:
    cols1= [
        "ensembl_gene_id", "external_gene_name", "hgnc_symbol", "hgnc_id",
    ]
    grp1 = make_groups(biomart_query_with_fallback(organism, cols1))
    adata= mapping(adata, ids, grp1, {"external_gene_name": "gene_symbol", 
    "hgnc_symbol":"hgnc_symbol" , "hgnc_id": "hgnc_id"})

    # Requete 2:
    cols2= [
        "ensembl_gene_id", "external_synonym", "entrezgene_id",
    "refseq_mrna",
    ]
    grp2 = make_groups(biomart_query_with_fallback(organism, cols2))
    adata= mapping(adata, ids, grp2, {
    "external_synonym":"alias_symbol" , "entrezgene_id": "NCBI_symbol", "refseq_mrna": "refseq_mrna"})

    # Requete 3:
    cols3= [
        "ensembl_gene_id", "uniprotswissprot", "uniprotsptrembl"
    ]
    grp3 = make_groups(biomart_query_with_fallback(organism, cols3))
    adata= mapping(adata, ids, grp3, {
    "uniprotswissprot": "uniprot_swissprot", "uniprotsptrembl": "uniprot_sptrembl"})

    return adata

def normalize_adata_in_memory(adata, use_biomart=False, organism="hsapiens", copy=True):
    """
    Cette fonction permet de préparer un AnnData pour la recherche de gènes sans modifier le fichier source.
    """
    if copy:
        adata = adata.copy()

    adata = gene_versions(adata)

    if use_biomart:
        adata = add_biomart_names(adata, organism=organism)

    return adata

def save_normalized_h5ad(adata, output_path, use_biomart=False, organism="hsapiens"):
    """
    Cette fonction permet de sauvegarder une copie normalisée de l'AnnData.
    """
    normalized = normalize_adata_in_memory(
        adata,
        use_biomart=use_biomart,
        organism=organism,
        copy=True
    )
    normalized.write_h5ad(output_path)
    return output_path

def export_normalized_h5ad_bytes(adata, use_biomart=False, organism="hsapiens"):
    """
    Cette fonction retourne le contenu binaire d'un .h5ad normalisé,
    prêt pour un download_button Streamlit.
    """

    normalized = normalize_adata_in_memory(
        adata,
        use_biomart=use_biomart,
        organism=organism,
        copy=True
    )

    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        normalized.write_h5ad(tmp_path)
        with open(tmp_path, "rb") as f:
            data = f.read()
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    return data




