import argparse
import scanpy as sc

from gene_search_utils import (
    search_gene_hits,
    format_hits,
    choose_var_name,
)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Recherche de gènes dans un .h5ad (tous IDs + gestion versions)")
    p.add_argument("file_path", help="Chemin vers le fichier .h5ad")
    p.add_argument("query", help="Le mot à rechercher")
    p.add_argument("--max", type=int, default=30, help="Nombre max de résultats affichés (30 pas défaut)")
    p.add_argument("--choose", type=int, default=None, help="Choisir directement un résultat")
    p.add_argument("--interactive", action="store_true", help="Demander le choix si plusieurs résultats")

    args = p.parse_args()

    adata = sc.read_h5ad(args.file_path)
    hits = search_gene_hits(adata, args.query, max_hits=args.max)

    print(format_hits(hits, max_lines=args.max))

    if hits.empty:
        raise SystemExit(0)

    if len(hits) == 1:
        print("\nChoisissez:", hits.iloc[0]["__var_name__"])
        raise SystemExit(0)

    if args.choose is not None:
        print("\nChoisissez:", choose_var_name(hits, args.choose))
        raise SystemExit(0)

    if args.interactive:
        k = int(input(f"\nChoisissez une ligne (1..{len(hits)}): ").strip())
        print("\nChoisissez:", choose_var_name(hits, k))