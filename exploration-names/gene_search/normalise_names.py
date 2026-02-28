import argparse
import scanpy as sc
import pandas as pd
from scanpy.queries import biomart_annotations

def gene_versions(adata):
    """
    Cette fonction rajoute les colonne base_name et has_versions au adata. 
    Si un gène a plus qu'une version on aura True dans has_versions,
    et le nom sans sufix (.X) dans base_name
    """
    names= pd.Index(adata.var_names.astype(str))
    adata.var["base_name"]=names.str.replace(r"\.", "", regex=True) #garde pas le (.X)
    adata.var["has_versions"]=names.str.contains(r"\.", regex=True) # a plus qu'une version ou pas
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


def add_biomart_names(adata, organism="hsapiens"):
    """
    Cette fonction ajoute plusierurs colonnes si disponible via BioMart:
    gene_symbol; hgnc_symbol; hgnc_id, alias_symbol, entrezgene_id (NCBI)
    refseq_mrna, refseq_peptide, uniprotswissprot, uniprotsptrembl
    Séparer en plusieurs requetes car BioMart permet la rechrche de 4 attribus à la fois.
    """
    if "gene_ids" not in adata.var.columns:
        raise ValueError ("gene_ids manquant")  #Si on utilise un autre fichier

    ids= adata.var["gene_ids"].astype(str)

    # Requete 1:
    cols1= [
        "ensembl_gene_id", "external_gene_name", "hgnc_symbol", "hgnc_id",
    ]
    grp1= make_groups(biomart_annotations(organism, cols1))
    adata= mapping(adata, ids, grp1, {"external_gene_name": "gene_symbol", 
    "hgnc_symbol":"hgnc_symbol" , "hgnc_id": "hgnc_id"})

    # Requete 2:
    cols2= [
        "ensembl_gene_id", "external_synonym", "entrezgene_id",
    "refseq_mrna",
    ]
    grp2= make_groups(biomart_annotations(organism, cols2))
    adata= mapping(adata, ids, grp2, {
    "external_synonym":"alias_symbol" , "entrezgene_id": "NCBI_symbol", "refseq_mrna": "refseq_mrna"})

    # Requete 3:
    cols3= [
        "ensembl_gene_id", "uniprotswissprot", "uniprotsptrembl"
    ]
    grp3= make_groups(biomart_annotations(organism, cols3))
    adata= mapping(adata, ids, grp3, {
    "uniprotswissprot": "uniprot_swissprot", "uniprotsptrembl": "uniprot_sptrembl"})

    return adata


if __name__=="__main__":
    parser= argparse.ArgumentParser(
        description= "Préparer les noms des gènes disponibles"
    )
    parser.add_argument(
        "file_path", help= "Chemin vers le fichier AnnData"
    )
    parser.add_argument(
        "--biomart", action="store_true" , help= "Ajoute gene_symbol via biomart (besoin d'une connection internet)" 
    )
    parser.add_argument(
        "--organism" ,default= "hsapiens", help= "Organism étudier (défaut: hsapiens)" 
    )
    parser.add_argument(
        "--save" ,action= "store_true", help= "Créer un nouveau fichier avec les nouvelles colonnes (.prep-versions.h5ad)" 
    )
    parser.add_argument(
        "--output" ,type = str, help= "Nom du fichier de sortie" 
    )

    args = parser.parse_args()

    adata=sc.read_h5ad(args.file_path)
    adata=gene_versions(adata)

    if args.biomart:
        adata= add_biomart_names(adata, organism=args.organism)

    new_columns= [col for col in [
        "gene_ids", "gene_symbol", "alias-symbol", 
        "hgnc_symbol", "hgnc_id",
        "NCBI_symbol", "refseq_mrna", "uniprot_swissprot",
        "uniprot_sptrembl", "base_name", "has_version"
    ]if col in adata.var.columns]

    print(adata.var[new_columns]. head())

    #Sauvegarder un nouveau fichier modifier si demander
    if args.save:
        if args.output:
            output_file=args.output
        else:
            output_file=args.file_path.replace(".h5ad", ".prep_versions.h5ad")
 
        adata.write_h5ad(output_file)
        print(f"Nouveau fichier créé: {output_file}")
    else:
        print("Aucun fichier créé. Modifications uniquement en mémoire")

    




