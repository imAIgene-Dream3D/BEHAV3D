"""Plot helpers for general behaviour visualisation.

Kept deliberately import-free: the heavy cluster plotting functions live in
:mod:`.cluster_plots` (scanpy, seaborn, sklearn) and are re-exported lazily via
PEP 562 ``__getattr__``, so importing a light sibling submodule such as
:mod:`.proportion_bars` no longer drags scanpy into the process. Accessing any
name below still works exactly as before -- it just pays the import at that
point rather than at package-import time.
"""

__all__ = [
    "plot_per_cluster_proportions",
    "plot_number_per_clusters",
    "plot_top_ranking_features",
    "plot_umap_feature_grid",
    "plot_clustering_feature_heatmap",
]


def __getattr__(name):
    if name in __all__:
        from . import cluster_plots

        return getattr(cluster_plots, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(list(globals()) + __all__)
