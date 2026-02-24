import argparse
import scanpy as sc
import pandas as pd

def gene_versions(adata):
    """
    Cette fonction rajoute les colonne base_name et has_versions au adata. 
    Si un gène a plus qu'une version on aura True dans has_versions,
    le nom complet avec le sufix (.X) dans var_name, et le nom sans sufix dans base_name
    """
    names= pd.Index(adata.var_names.astype(str))
    adata.var["var_name"]=name #garde le nom originale
    adata.var["base_name"]=name.str.replace(r"\.", "", regex=True) #garde pas le (.X)
    adata.var["has_versions"]=name.str.contains(r"\.", "", regex=True) # a plus qu'une version ou pas
    return adata

if __name__=="__main__":
    parser= argparse.ArgumentParser(
        description= "Préparer base_name pour proposer les versions des gènes"
    )
    parser.add_argument(
        "file_path", help= "Chemin vers le fichier AnnData"
    )
    parser.add_argument(
        "-o", "--out" ,type=int, default=None, help= "Nom du fichier de sortie" 
    )
    args = parser.parse_args()

    adata=sc.read_h5ad(args.file_path)
    adata=gene_versions(adata)

    output_file=args.out or args.file_path.replace(".h5ad", ".prep_versions.h5ad")
    adata.write_h5ad(output_file)
    print(output_file)




