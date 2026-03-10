import argparse
import scanpy as sc

from gene_search_utils import (
    resolve_gene_to_var_name,
    search_gene_partial,
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
    p.add_argument("--show-partial",action="store_true",help="Afficher les résultats partiels même si un match exact est trouvé",)

    args = p.parse_args()

    adata = sc.read_h5ad(args.file_path)

    var_name, hits, suggestions, match_type = resolve_gene_to_var_name(
    adata,
    args.query,
    choice=args.choose,
    max_hits=args.max
)

# Si un résultat existe
if var_name is not None:

    if match_type == "exact":
        print("Résultat exact :\n")
    elif match_type == "partial":
        print("Résultat partiel :\n")

    print(format_hits(hits, max_lines=args.max))
    print("\nGène sélectionné :", var_name)

    if args.show_partial and match_type == "exact":
        partial_hits = search_gene_partial(adata, args.query, max_hits=args.max)

        if not partial_hits.empty:
            partial_hits = partial_hits[partial_hits["__var_name__"] != var_name]

            if not partial_hits.empty:
                print("\nAutres résultats partiels possibles :\n")
                print(format_hits(partial_hits, max_lines=args.max))

    raise SystemExit(0)

# Si plusieurs résultats trouvés
if hits is not None:
    if match_type == "exact":
        print("Plusieurs résultats exacts trouvés :\n")
    elif match_type == "partial":
        print("Plusieurs résultats partiels trouvés :\n")

    print(format_hits(hits, max_lines=args.max))

    if args.interactive:
        k = int(input(f"\nChoisir une ligne (1..{len(hits)}): ").strip())
        print("\nNom sélectionné :", choose_var_name(hits, k))

    raise SystemExit(0)

# Si aucune correspondance mais suggestions proches
if suggestions:
    print("Aucun résultat exact.")
    print("\nSuggestions possibles :")
    for s in suggestions:
        print("-", s)
    raise SystemExit(0)

print("Aucun résultat.")