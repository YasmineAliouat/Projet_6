import argparse
import scanpy as sc
from normalise_names import normalize_adata_in_memory

parser = argparse.ArgumentParser(
    description="Préparer les noms des gènes disponibles"
)
parser.add_argument("file_path", help="Chemin vers le fichier AnnData")
parser.add_argument("--biomart", action="store_true")
parser.add_argument("--organism", default="hsapiens")
parser.add_argument("--save", action="store_true")
parser.add_argument("--output", type=str)

args = parser.parse_args()

adata = sc.read_h5ad(args.file_path)
adata = normalize_adata_in_memory(
    adata,
    use_biomart=args.biomart,
    organism=args.organism,
    copy=False
)

new_columns = [col for col in [
    "gene_ids",
    "gene_symbol",
    "alias_symbol",
    "hgnc_symbol",
    "hgnc_id",
    "NCBI_symbol",
    "refseq_mrna",
    "uniprot_swissprot",
    "uniprot_sptrembl",
    "base_name",
    "has_versions",
] if col in adata.var.columns]

print(adata.var[new_columns].head())

if args.save:
    output_file = args.output or args.file_path.replace(".h5ad", ".prep_versions.h5ad")
    adata.write_h5ad(output_file)
    print(f"Nouveau fichier créé : {output_file}")
else:
    print("Aucun fichier créé. Modifications uniquement en mémoire.")