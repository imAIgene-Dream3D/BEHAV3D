# 📊 Behavioural Analysis

```{admonition} GUI coming soon
:class: warning

The napari **Analysis tab** is currently a stub ([behav3d/napari/_stubs.py](../api/behav3d.napari)). This page documents the underlying Python entry point you can call from a script today; it will be expanded with full GUI documentation once the sub-tab lands.
```

Behavioural Analysis clusters tracks by **motion + intensity + contact pattern** to discover recurring behavioural phenotypes — e.g. "static, alive, in contact with organoid" vs "fast-moving, away from organoid". The pipeline is **cell-type agnostic** (despite the historical "T-cell" name in the source).

## What it does

1. Reads the **filtered + summarised** track features written by the [Filtering](filtering) tab.
2. Selects a configurable subset of feature columns.
3. Normalises across samples / conditions.
4. Computes pairwise **DTW (Dynamic Time Warping)** distances between tracks.
5. Embeds the distance matrix with **UMAP**.
6. Clusters the UMAP embedding with **Leiden** community detection.
7. Writes per-track cluster IDs + diagnostic plots.

```{mermaid}
flowchart LR
    A["Filtered + summarised<br/>track features CSV"] --> B["Pick feature columns"]
    B --> C["Per-condition normalisation"]
    C --> D["DTW pairwise distances"]
    D --> E["UMAP embedding"]
    E --> F["Leiden clustering"]
    F --> G["Per-track cluster IDs<br/>+ heatmap + UMAP plot"]
```

## Python entry point

The function is `run_tcell_analysis` ([behav3d/analysis/tcell_analysis.py:107](../api/behav3d.analysis)) — despite the legacy name it accepts a `cell_type` parameter and works for any cell type with a filtered features CSV.

```python
from behav3d.analysis import run_tcell_analysis

run_tcell_analysis(
    output_dir="D:/data/exp1/behav3d_out",
    cell_type="tcell",
    columns_to_use=[
        "mean_square_displacement",
        "speed",
        "mean_dead_dye",
        "tcell_contact",
        "organoid_contact",
    ],
    columns_to_normalize=[
        "mean_square_displacement",
        "speed",
    ],
)
```

Common parameters (full signature in the source):

| Parameter | Description |
|---|---|
| `output_dir` | The BEHAV3D output directory (the same one you set in Data Preparation). |
| `cell_type` | The cell type to analyse (e.g. `tcell`, `macrophage`). |
| `df_tracks_path` | Optional explicit override for the filtered tracks CSV; defaults to `<output_dir>/analysis/<cell_type>/track_features/...filtered.csv`. |
| `df_tracks_summarized_path` | Optional explicit override for the summarised CSV. |
| `columns_to_use` | Feature columns fed into the DTW distance. |
| `columns_to_normalize` | Subset of `columns_to_use` normalised per condition before DTW. |

## Outputs

```
<output_dir>/analysis/<cell_type>/results/
    BEHAV3D_<cell_type>_UMAP_clusters.csv
    BEHAV3D_<cell_type>_UMAP_clusters.pdf
    BEHAV3D_<cell_type>_UMAP_cluster_feature_heatmap.pdf
```

Plus a `per_sample/` subdirectory with per-sample breakdowns.

## What the future GUI will look like

When the **Analysis → Behavioural** sub-tab is implemented, it will expose:

- The `columns_to_use` / `columns_to_normalize` checklist via the GUI.
- UMAP / DTW hyperparameters (n_neighbors, min_distance, n_clusters).
- A live preview of the UMAP after clustering.
- A queue button to add a `BEHAVIOURAL_ANALYSIS` step to the [Processing Queue](../plugin_essentials/processing_queue).

This page will be replaced with the full GUI walkthrough once that lands.

## See also

- [Filtering](filtering) — produces the input CSVs.
- [Backprojection](backprojection) — visualise clusters back onto the raw image in napari.
- Source: [behav3d/analysis/tcell_analysis.py:107](../api/behav3d.analysis).
