import scanpy as sc
import argparse
import pandas as pd

pd.set_option("display.max_columns", None)
pd.set_option("display.width", None)
pd.set_option("display.max_colwidth", None)



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

    colonnes=[ "gene_ids", "gene_symbol", "alias_symbol", "hgnc_symbol",
     "hgnc_id", "NCBI_symbol", "refseq_mrna", "uniprot_swissprot", "uniprot_sptrembl", "base_name", "has_versions" ]

    cols= [col for col in colonnes if col in adata.var.columns]
    print(adata.var[cols].head(n))

    #Pour savoir les noms disponibles pour chaque attributs
    print ( "Nombre de gènes avec gene_symbol:", int(adata.var["gene_symbol"].notna().sum()), "/", adata.n_vars)
    print ( "Nombre de gènes avec alias_symbol:", int(adata.var["alias_symbol"].notna().sum()), "/", adata.n_vars)
    print ( "Nombre de gènes avec hgnc_symbol:", int(adata.var["hgnc_symbol"].notna().sum()), "/", adata.n_vars)
    print ( "Nombre de gènes avec hgnc_id:", int(adata.var["hgnc_id"].notna().sum()), "/", adata.n_vars)
    print ( "Nombre de gènes avec NCBI_symbol:", int(adata.var["NCBI_symbol"].notna().sum()), "/", adata.n_vars)
    print ( "Nombre de gènes avec refseq_mrna:", int(adata.var["refseq_mrna"].notna().sum()), "/", adata.n_vars)
    print ( "Nombre de gènes avec uniprot_swissprot:", int(adata.var["uniprot_swissprot"].notna().sum()), "/", adata.n_vars)
    print ( "Nombre de gènes avec uniprot_sptrembl:", int(adata.var["uniprot_sptrembl"].notna().sum()), "/", adata.n_vars)


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