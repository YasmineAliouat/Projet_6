from anndata import AnnData
from scipy import sparse

import scanpy as sc


class AnnDataSummary:
    def __init__(
        self,
        n_obs,
        n_vars,
        is_sparse_X,
        has_raw,
        layers,
        obs_columns_preview,
        var_columns_preview,
        obsm_keys,
        has_umap,
    ):
        self.n_obs = n_obs
        self.n_vars = n_vars
        self.is_sparse_X = is_sparse_X
        self.has_raw = has_raw
        self.layers = layers
        self.obs_columns_preview = obs_columns_preview
        self.var_columns_preview = var_columns_preview
        self.obsm_keys = obsm_keys
        self.has_umap = has_umap

    def __repr__(self):
        return (
            f"AnnDataSummary("
            f"n_obs={self.n_obs}, "
            f"n_vars={self.n_vars}, "
            f"is_sparse_X={self.is_sparse_X}, "
            f"has_raw={self.has_raw}, "
            f"layers={self.layers}, "
            f"has_umap={self.has_umap}"
            f")"
        )


def load_anndata(h5ad_path: str) -> AnnData:
    """
    Charge un fichier .h5ad et renvoie un objet AnnData.
    """
    try:
        adata = sc.read_h5ad(h5ad_path)
    except Exception as e:
        raise RuntimeError(
            f"Cannot read the .h5ad file: {h5ad_path}\nDetail: {e}"
        ) from e

    return adata

def validate_anndata(adata: AnnData, *, require_var_names: bool = True) -> None:
    """
    Vérifications.
    """
    if adata is None:
        raise ValueError("AnnData is None.")

    if adata.X is None:
        raise ValueError("AnnData.X is missing (no expression matrix available).")

    if adata.n_obs <= 0 or adata.n_vars <= 0:
        raise ValueError(f"Invalid dimensions: n_obs={adata.n_obs}, n_vars={adata.n_vars}")

    if require_var_names:
        if adata.var_names is None or len(adata.var_names) == 0:
            raise ValueError("adata.var_names is empty.")
        if not adata.var_names.is_unique:
            raise ValueError("adata.var_names contains duplicates (must be unique).")

    if adata.obs_names is not None and not adata.obs_names.is_unique:
        raise ValueError("adata.obs_names contains duplicates (must be unique).")


def summarize_anndata(adata: AnnData, *, preview_n: int = 15,) -> AnnDataSummary:
    """
    Renvoie un résumé structuré (utile pour debug, logs, ou interface).
    """
    is_sparse_X = sparse.issparse(adata.X)  # sparse ou dense ?
    has_raw = adata.raw is not None # données brutes ?
    layers = list(adata.layers.keys())
    obs_cols = list(adata.obs.columns)[: min(preview_n, adata.obs.shape[1])]    # aperçu des colonnes
    var_cols = list(adata.var.columns)[: min(preview_n, adata.var.shape[1])]
    obsm_keys = list(adata.obsm.keys())
    has_umap = "X_umap" in obsm_keys

    return AnnDataSummary(
        n_obs=adata.n_obs,
        n_vars=adata.n_vars,
        is_sparse_X=is_sparse_X,
        has_raw=has_raw,
        layers=layers,
        obs_columns_preview=obs_cols,
        var_columns_preview=var_cols,
        obsm_keys=obsm_keys,
        has_umap=has_umap,
    )