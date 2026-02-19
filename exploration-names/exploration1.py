import scanpy as sc
import argparse

def explore_gene_names(path: str, n: int = 10) -> None:
    """
    Fonction pour afficher les n premiers var_names et gene_ids.
    """
    adata=sc.read_h5ad(path)

    print(f" Premies {n} var_names:")
    print(adata.var_names[:n])

    print(f" Premies {n} gene_ids:")
    print(adata.var['gene_ids'].head(n))

    print("var_names uniques?:", adata.var_names.is_unique)

    x_versions= adata.var_names.str.contains(r"\.").sum()
    print("Nombre de noms avec plusieurs versions (.x)", x_versions)

if __name__=="__main__":
    parser= argparse.ArgumentParser(
        description= "Explorer les noms des gènes dans un fichier AnnData"
    )
    parser.add_argument(
        "file_path", help= "Chemin vers le fichier AnnData"
    )
    parser.add_argument(
        "-n", type=int, default=10, help= "Nombre de gènes à afficher (10 par défaut)" 
    )

    args = parser.parse_args()
    explore_gene_names(args.file_path, args.n)