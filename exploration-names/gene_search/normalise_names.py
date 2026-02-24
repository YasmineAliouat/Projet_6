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
    adata.var["has_versions"]=names.str.contains(r"\.", "", regex=True) # a plus qu'une version ou pas
    return adata

def gene_symbol(adata, organism="hsapiens"):
    """
    Cette fonction ajoute une colonne gene_symbol qui contient
    les noms des gène dans biomart.
    """
    if "gene_ids" not in adata.var.columns:
        raise ValueError ("gene_ids manquant")
    annot=biomart_annotations(organism, ["ensembl_gene_id", "external_gene_name"])
    mapping= dict(zip(annot["ensembl_gene_id"], annot ["external_gene_name"]))
    adata.var["gene_symbol"]= adata.var["gene_ids"].map(mapping)
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
        adata= gene_symbol(adata, organism=args.organism)
    print(adata.var[[ col for col in ["gene_ids", "gene_symbol","base_name", "has_versions"]if col in adata.var.columns]]. head())


    if args.save:
        if args.output:
            output_file=args.output
        else:
            output_file=args.file_path.replace(".h5ad", ".prep_versions.h5ad")
 
        adata.write_h5ad(output_file)
        print(f"Nouveau fichier créé: {output_file}")
    else:
        print("Aucun fichier créé. Modifications uniquement en mémoire")

    




