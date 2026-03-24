import pandas as pd
import plotly.express as px

def plot_gene_coexpression_scatter_plot(adata, gene_a: str, gene_b: str):

    # Extraire l'expression des gènes A et B (log-transformée)
    expr_a = adata[:, gene_a].X.toarray().flatten()
    expr_b = adata[:, gene_b].X.toarray().flatten()

    # Créer un DataFrame pour Plotly
    df = pd.DataFrame({
        gene_a: expr_a,
        gene_b: expr_b,
        "Cluster": adata.obs["louvain"] 
    })

    ### Scatter plot avec densité (hexbin)
    fig = px.scatter(
        df,
        x=gene_a,
        y=gene_b,
        color="Cluster",  # Optionnel : colorer par cluster
        marginal_x="histogram",
        marginal_y="histogram",
        title=f"Co-expression de {gene_a} et {gene_b}",
        labels={gene_a: f"Expression {gene_a} (log CPM)", gene_b: f"Expression {gene_b} (log CPM)"}
    )
    fig.update_traces(
        marker=dict(size=5, opacity=0.5),
        selector=dict(type='scatter')
    )
    fig.show()