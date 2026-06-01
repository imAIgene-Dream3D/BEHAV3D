# ☠ Death & Morphology Analysis

```{admonition} GUI coming soon
:class: warning

The napari **Analysis tab** is currently a stub ([behav3d/napari/_stubs.py](../api/behav3d.napari)). This page documents the underlying Python entry point you can call from a script today; it will be expanded with full GUI documentation once the sub-tab lands.
```

Death & Morphology Analysis quantifies, per condition / line, how each organoid (or other cell type) population evolves over time in terms of **death dynamics** and **morphology trajectories**. The pipeline is **cell-type agnostic**: despite the historical "organoid" name in the source, it works on any cell type with the relevant features extracted.

## What it does

1. Reads the **filtered + summarised** features for the chosen cell type.
2. Aggregates the per-track death-fraction trajectory per condition / line.
3. Fits per-track death dynamics (sigmoid / step / monotonic models).
4. Clusters morphology trajectories (volume, sphericity, elongation over time).
5. Writes per-track classifications + condition-level plots.

```{mermaid}
flowchart LR
    A["Filtered + summarised<br/>features CSV"] --> B["Death dynamics<br/>per track"]
    A --> C["Morphology trajectory<br/>clustering"]
    B --> D["Per-condition<br/>survival / death curves"]
    C --> D
    D --> E["Quality plots +<br/>per-track classifications"]
```

## Python entry point

The function is `run_organoid_analysis` ([behav3d/analysis/organoid_analysis.py:34](../api/behav3d.analysis)) — despite the name it takes an `org_type` parameter and works for any cell type with morphology + death features:

```python
from behav3d.analysis import run_organoid_analysis

run_organoid_analysis(
    output_dir="D:/data/exp1/behav3d_out",
    org_type="organoid",      # or any other cell type with features extracted
    dead_perc_threshold=0.1,  # same threshold as the Feature Extraction tab
)
```

Common parameters (full signature in the source):

| Parameter | Description |
|---|---|
| `output_dir` | BEHAV3D EXPLORER output directory. |
| `org_type` | Cell type to analyse (default `organoid`). |
| `dead_perc_threshold` | Dead-mask % threshold for the death classification (typically matched to the one used during Feature Extraction). |
| `df_tracks_path` | Optional explicit override for the input CSV. |

## Outputs

```
<output_dir>/analysis/<cell_type>/results/
    BEHAV3D_<cell_type>_death_dynamics.csv
    BEHAV3D_<cell_type>_death_curves.pdf
    BEHAV3D_<cell_type>_morphology_clusters.csv
    BEHAV3D_<cell_type>_morphology_heatmap.pdf
```

Plus per-sample plots in `per_sample/` and quality plots in `quality_control/`.

## What the future GUI will look like

When the **Analysis → Death & Morphology** sub-tab is implemented, it will expose:

- A cell-type sub-tab structure (mirroring Feature Extraction and Filtering).
- Death threshold spinner — pre-filled from the Feature Extraction setting.
- Morphology feature checklist.
- Live preview of death curves per condition.
- A queue button to add a `DEATH_MORPHOLOGY_ANALYSIS` step.

This page will be replaced with the full GUI walkthrough once that lands.

## See also

- [Feature Extraction](feature_extraction) — must include `death` and `morphology` features for this analysis.
- [Filtering](filtering) — produces the input CSV.
- [Backprojection](backprojection) — visualise death classifications back onto the raw image.
- Source: [behav3d/analysis/organoid_analysis.py:34](../api/behav3d.analysis).
