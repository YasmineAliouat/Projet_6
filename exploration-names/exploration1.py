import scanpy as sc
import pandas as pd

def explore_gene_names(path):
    """
    Fonction pour une premiére visualisation des noms des génes et de leurs Ids
    """
    adata = sc.read_h5ad(path)

    print("Les premiers var_names")
    print(adata.var_names[:5])

    if 'gene_ids' in adata.var.columns:
        print("\nLes premiers gene_ids")
        print(adata.var['gene_ids'].head(5))

    print("\nNombre total de gènes", adata.n_vars)

    print("\nvar_names uniques?", adata.var_names.is_unique)

    has_versions = adata.var_names.str.contains("\.").sum()
    print("\nNombre de noms avec plusieurs versions", has_versions)

if __name__ == "__main__":
    explore_gene_names("/home/lili/Documents/M1/S2/Projet6/adata_3583.h5ad")

