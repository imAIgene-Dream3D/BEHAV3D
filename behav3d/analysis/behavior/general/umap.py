import numpy as np
import pandas as pd
import umap

from behav3d.analysis.behavior.utils import _handle_nan_in_distance_matrix


def fit_umap(
    dtw_distance_matrix,
    config=None,
    umap_n_neighbors=None,
    umap_minimal_distance=None,
    random_state=None,
):
    print("- Fitting the dynamic time warping to a UMAP")
    assert config is not None or all(
        [umap_n_neighbors, umap_minimal_distance]
    ), "Either 'config' or 'umap_n_neighbors and umap_minimal_distance' must be supplied"

    if umap_n_neighbors is None or umap_minimal_distance is None:
        assert config is not None, "Provide config or both umap_n_neighbors and umap_minimal_distance"
        umap_n_neighbors = config["umap_n_neighbors"]
        umap_minimal_distance = config["umap_minimal_distance"]

    dtw_matrix_values = dtw_distance_matrix.values.copy()
    dtw_matrix_values = _handle_nan_in_distance_matrix(
        dtw_matrix_values, context="distance matrix (pre-UMAP)"
    )

    umap_model = umap.UMAP(
        n_components=2,
        n_neighbors=umap_n_neighbors,
        min_dist=umap_minimal_distance,
        init="random",
        random_state=random_state,
        metric="precomputed",
    )
    umap_embedding = umap_model.fit_transform(dtw_matrix_values)
    umap_embedding = pd.DataFrame(umap_embedding, columns=["UMAP1", "UMAP2"])
    umap_embedding[["TrackID", "sample_name"]] = pd.DataFrame(
        [string.split("--") for string in dtw_distance_matrix.index]
    )
    umap_embedding["TrackID"] = umap_embedding["TrackID"].astype(np.int64)
    return umap_embedding
