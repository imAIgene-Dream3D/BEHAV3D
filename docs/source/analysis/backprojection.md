# 🎨 Backprojection

```{admonition} GUI coming soon
:class: warning

The napari **Analysis tab** is currently a stub ([behav3d/napari/_stubs.py](../api/behav3d.napari)). This page documents the underlying Python entry points you can call from a script today; it will be expanded with full GUI documentation once the sub-tab lands.
```

Backprojection takes the **results of an analysis step** (behavioural cluster IDs, death classifications, interaction event flags) and paints them back onto the original tracked-segments zarr — colouring each cell by its cluster / class / event status — so you can visualise the analysis output directly in napari, frame by frame.

## What it does

For every (sample × cell type × timepoint × track_id) it looks up the analysis result (e.g. cluster ID) and writes a new label image where the pixel value is the analysis result instead of the original track ID. The result is a 4-D `(T, Z, Y, X)` integer zarr that opens in napari as a Labels layer.

```{mermaid}
flowchart LR
    A["Analysis CSV<br/>(cluster IDs, classes, events)"] --> B["Backprojection"]
    T["Tracked-segments zarr"] --> B
    B --> Z["Backprojected zarr<br/>(pixel value = cluster ID)"]
    Z --> N["napari Labels layer<br/>+ colour legend"]
```

## Python entry points

The main functions are in [behav3d/analysis/backprojection.py](../api/behav3d.analysis):

| Function | Purpose |
|---|---|
| `backproject_mean_features_behav3d` | Backproject the per-track *mean* of a feature column (e.g. mean speed). |
| `backproject_time_features_behav3d` | Backproject a per-timepoint feature column (e.g. instantaneous speed). |
| `backproject_columns` | Generic backprojection of any columns from a clustered CSV. |
| `view_napari` | Convenience opener — launches napari with the relevant raw image + backprojected layer pre-loaded. |

Example:

```python
from behav3d.analysis.backprojection import backproject_columns, view_napari

backproject_columns(
    output_dir="D:/data/exp1/behav3d_out",
    cell_type="tcell",
    df_clustered_path="…/BEHAV3D_tcell_UMAP_clusters.csv",
    columns=["umap_cluster"],
    sample_name="Sample_A1",
)

view_napari(
    output_dir="D:/data/exp1/behav3d_out",
    cell_type="tcell",
    sample_name="Sample_A1",
    column="umap_cluster",
)
```

## Outputs

```
<output_dir>/analysis/<cell_type>/backprojection/
    <sample_name>_<cell_type>_<column>_backprojected.zarr
```

One zarr per sample × column. Same `(T, Z, Y, X)` layout as the other BEHAV3D zarrs.

## What the future GUI will look like

When the **Analysis → Backprojection** sub-tab is implemented, it will expose:

- A column picker — choose which analysis output to backproject (UMAP cluster, death class, etc.).
- A sample selector.
- One-click "Backproject & Open in napari" that runs the projection and opens the result alongside the raw image with a colour legend.
- A queue button to add a `BACKPROJECT` step.

This page will be replaced with the full GUI walkthrough once that lands.

## See also

- [Behavioural Analysis](behavioural), [Death & Morphology Analysis](death_morphology), [Interaction Analysis](interaction) — produce the CSVs that backprojection consumes.
- [Visualization tab](../plugin_essentials/visualization) — used to open the backprojected zarr alongside raw + segments.
- Source: [behav3d/analysis/backprojection.py](../api/behav3d.analysis).
