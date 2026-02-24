import argparse
import scanpy as sc
import pandas as pd

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

if __name__=="__main__":
    parser= argparse.ArgumentParser(
        description= "Préparer base_name pour proposer les versions des gènes"
    )
    parser.add_argument(
        "file_path", help= "Chemin vers le fichier AnnData"
    )

    parser.add_argument(
        "--output" ,type = str, help= "Nom du fichier de sortie" 
    )
    parser.add_argument(
        "--save" ,action= "store_true", help= "Créer un nouveau fichier avec les nouvelles colonnes (.prep-versions.h5ad)" 
    )

    args = parser.parse_args()

    adata=sc.read_h5ad(args.file_path)
    adata=gene_versions(adata)
    print(adata.var[["base_name", "has_versions"]]. head())


    if args.save:
        if args.output:
            output_file=args.output
        else:
            output_file=args.file_path.replace(".h5ad", ".prep_versions.h5ad")
 
        adata.write_h5ad(output_file)
        print(f"Nouveau fichier créé: {output_file}")
    else:
        print("Aucun fichier créé. Modifications uniquement en mémoire")

    




