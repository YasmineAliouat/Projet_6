import scanpy as sc 
from backend.signature_utils import load_signature_genes

#Chargement du dataset
adata =sc.read_h5ad("data/adata_3583.h5ad")

#Chargement de la signature génique 
siganture = load_signature_genes("data/signature.txt")

#Vérification des gènes présents et manquants
genes_present = [g for g in signature if g in adata.var_names]
genes_missing = [g for g in signature if g not in adata.var_names]

print("Gènes présents:", genes_present)
print("Gènes manquants:", genes_missing)



from backend.signature_utils import compute_signature_score
score_name = "signature_score"
adata = compute_signature_score(adata, genes_present, score_name)

print(adata.obs[[score_name]].head())

adata.write("data/adata_3583_with_signature.h5ad")
print("Dataset sauvegardé avec le score de signature")