# Chargement des librairies
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse as sp

# Fonction pour lire un fichier .h5ad,
    # analyse sa structure,
    # afficher un rapport lisible

def explore_h5ad(path: str, n_preview: int = 10) -> None:
    adata = sc.read_h5ad(path)

    print("=== AnnData report ===")
    print(f"Fichier: {path}")
    print(f"Dimensions: {adata.n_obs} cellules x {adata.n_vars} gènes")
    print(f"X type: {type(adata.X)}")
    print(f"X sparse: {sp.issparse(adata.X)}")

    # --- var / gene names ---
    print("\n--- Genes (.var) ---")
    print(f"var columns (preview): {list(adata.var.columns)[:min(n_preview, len(adata.var.columns))]}")
    print(f"var_names preview: {list(adata.var_names[:min(n_preview, adata.n_vars)])}")

   # --- obs / metadata ---
    print("\n--- Cell metadata (.obs) ---")
    print(f"obs columns: {len(adata.obs.columns)}")
    print(f"obs columns (preview): {list(adata.obs.columns)[:min(25, len(adata.obs.columns))]}")
   
   # --- embeddings / UMAP ---
    print("\n--- Embeddings (.obsm) ---")
    obsm_keys = list(adata.obsm.keys())
    print(f"obsm keys: {obsm_keys}")
    has_umap = "X_umap" in obsm_keys
    print(f"UMAP présente (obsm['X_umap']): {has_umap}")
    if has_umap:
        print(f"UMAP shape: {adata.obsm['X_umap'].shape}")

    # --- raw / layers ---
        print("\n--- raw & layers ---")
        has_raw = adata.raw is not None # vérifie si raw existe
        print(f"raw présent: {has_raw}")
        # Si raw existe, afficher infos
        if has_raw:
            print(f"raw dimensions: {adata.raw.n_obs} x {adata.raw.n_vars}")
            preview = min(n_preview, adata.raw.n_vars)
            print(f"raw var_names preview: {list(adata.raw.var_names[:preview])}")

        layer_keys = list(adata.layers.keys())
        print(f"layers: {layer_keys}")

        for lk in layer_keys[:min(n_preview, len(layer_keys))]:
            print(f"layer '{lk}' type: {type(adata.layers[lk])}")
        # affiche le type de chaque layer
    
    # --- quick check: data looks log? ---
    print("\n--- Quick sanity checks ---")
    try:
        if hasattr(adata.X, "data"):  # sparse
            sample = adata.X.data[:1000] # sélection d'un échantillon
        else:
            sample = np.ravel(adata.X)[:1000]
        sample = np.asarray(sample, dtype=float)

        max_val = np.max(sample)
        median_val = np.median(sample)
        print(f"X sample: median={median_val:.2f}, max={max_val:.2f}")

        looks_logged = max_val < 30

        print(f"X semble log-transformée: {looks_logged}")

    except Exception as e:
        print(f"Vérification impossible: {e}")

    print("\n=== End report ===")

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("h5ad_path", help="Chemin vers le fichier .h5ad")
    args = p.parse_args()
    explore_h5ad(args.h5ad_path)