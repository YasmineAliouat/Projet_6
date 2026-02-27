import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_gene_coexpression_heatmap(adata, gene_a: str, gene_b: str, gene_c: str, gene_d: str, gene_e: str):

    ### Heatmap de corrélation

    # Sélectionner quelques gènes d'intérêt
    genes = gene_a, gene_b, gene_c, gene_d, gene_e
    expr_matrix = adata[:, genes].X.toarray().T  # Transposer pour avoir les gènes en lignes

    # Calculer la matrice de corrélation
    corr_matrix = np.corrcoef(expr_matrix)
    corr_df = pd.DataFrame(corr_matrix, index=genes, columns=genes)

    # Heatmap
    plt.figure(figsize=(8, 6))
    sns.heatmap(
        corr_df,
        annot=True,
        cmap="coolwarm",
        vmin=-1, vmax=1,
        xticklabels=genes,
        yticklabels=genes
    )
    plt.title("Matrice de corrélation entre gènes")
    plt.show()