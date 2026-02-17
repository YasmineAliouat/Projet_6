import scanpy as sc
from plot_gene_expression import plot_gene_expression

adata = sc.read_h5ad("/Users/fannycauchois/Projet/projet-6/data/adata_3583.h5ad")
gene_test = input("Nom du gène : ")
print("Test gène:", gene_test)

if gene_test not in adata.var_names:
    print(f"Gène introuvable : {gene_test}")
    print("Exemples de gènes disponibles :", list(adata.var_names[:10]))
else:
    plot_gene_expression(adata, gene_test)
