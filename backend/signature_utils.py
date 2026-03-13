import scanpy as sc
import matplotlib.pyplot as plt

def load_signature_genes(path):
    """
    Charge une signature génique depuis un fichier texte.
    Le fichier doit contenir un nom de gène par ligne.
    """
    with open(path, "r") as f:
        genes = [line.strip() for line in f if line.strip()]
    return genes


def compute_signature_score(adata, gene_list, score_name):
    """
    Calcule un score de signature génique par cellule.

    Paramètres
    ----------
    adata : AnnData
        Objet AnnData contenant les données scRNA-seq
    gene_list : list
        Liste des gènes composant la signature
    score_name : str
        Nom du score à ajouter dans adata.obs

    Retour
    ------
    adata : AnnData
        Objet AnnData avec le score ajouté dans adata.obs
    """

    # Vérifier quels gènes de la signature sont présents
    genes_present = [g for g in gene_list if g in adata.var_names]

    if len(genes_present) == 0:
        raise ValueError(
            f"Aucun gène de la signature '{score_name}' n'est présent dans les données"
        )

    # Calcul du score de signature avec Scanpy
    sc.tl.score_genes(
        adata,
        gene_list=genes_present,
        score_name=score_name
    )

    return adata

def plot_signature_score(adata, gene_list, score_name="signature_score"):
    """
    Calcule et affiche le score de signature sur UMAP.
    """

    # calcul du score
    compute_signature_score(adata, gene_list, score_name)

    # plot UMAP
    sc.pl.umap(
        adata,
        color=score_name,
        show=False
    )

    plt.show()