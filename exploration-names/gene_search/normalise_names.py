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

def add_biomart_names(adata, organism="hsapiens"):
    """
    Cette fonction ajoute plusierurs colonnes si disponible via BioMart:
    gene_symbol; hgnc_symbol; hgnc_id, alias_symbol, entrezgene_id (NCBI)
    refseq_mrna, refseq_peptide, uniprotswissprot, uniprotsptrembl
    """
    if "gene_ids" not in adata.var.columns:
        raise ValueError ("gene_ids manquant")  #Si on utilise un autre fichier

    columns= [
        "ensembl_gene_id", "external_gene_name", "hgnc_symbol", "hgnc_id",
        "external_synonym", "entrezgene_id", "refseq_mrna", "refseq_peptide", 
        #"uniprotswissprot", "uniprotsptrembl",
    ]

    annot=biomart_annotations(organism, columns)

    existant= [col for col in columns if col in annot.columns ] #Ne garder que les colonnes existantes

    grouped_annot= (annot.groupby("ensembl_gene_id").agg(lambda x: "|".join (x.dropna().astype(str).unique()))) #Certain Ids ensemble pointe vers plusiers noms
    #adata.var["alias_symbol"]= adata.var["gene_ids"].map(grouped_annot["alias_symbol"])

    rename={
        "external_gene_name": "gene_symbol", "entrezgene_id":"NCBI_gene_id",
        "uniprotswissprot":"uniprot_swissprot", "uniprotsptrembl":"uniprot_trembl",

    }

    for column in existant:
        if column == "ensembl_gene_id":
            continue
        out_column= rename.get(column, column)
        adata.var[out_column]=adata.var["gene_ids"].map(grouped_annot[column])

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
        "gene_ids", "gene_symbol", #"alias-symbol", 
        "hgnc_symbol",
        "NCBI_gene_id", "refseq_mrna", "uniprot_swissprot", "base_name", "has_version"
    ]if col in adata.var.columns]
    print(adata.var[new_columns]. head())


    if args.save:
        if args.output:
            output_file=args.output
        else:
            output_file=args.file_path.replace(".h5ad", ".prep_versions.h5ad")
 
        adata.write_h5ad(output_file)
        print(f"Nouveau fichier créé: {output_file}")
    else:
        print("Aucun fichier créé. Modifications uniquement en mémoire")

    




