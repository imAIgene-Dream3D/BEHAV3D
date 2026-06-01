# 🤝 Interaction Analysis

```{admonition} GUI coming soon
:class: warning

The napari **Analysis tab** is currently a stub ([behav3d/napari/_stubs.py](../api/behav3d.napari)). This page documents the underlying Python entry point you can call from a script today; it will be expanded with full GUI documentation once the sub-tab lands.
```

Interaction Analysis quantifies how a focal cell type (typically an organoid) responds to **physical interactions** with one or more partner cell types (typically immune cells). It produces per-interaction event tables and time-resolved interaction statistics.

## What it does

1. Reads the filtered features for the focal cell type and each interacting cell type.
2. Identifies **interaction events** — periods during which a focal cell is in contact with one or more partner cells, using the contact threshold from Feature Extraction.
3. Aligns each interaction event to a common time origin (e.g. first contact).
4. Computes the focal cell's response — death progression, morphology changes, intensity transfer — across the interaction window.
5. Writes per-event tables and condition-level summary plots.

## Python entry point

The function is `run_interaction_analysis` ([behav3d/analysis/interaction_analysis.py:26](../api/behav3d.analysis)):

```python
from behav3d.analysis import run_interaction_analysis

run_interaction_analysis(
    output_dir="D:/data/exp1/behav3d_out",
    cell_type="organoid",                         # focal cell type
    interacting_cell_types=["tcell", "macrophage"],
    dead_threshold=0.02,
)
```

Parameters:

| Parameter | Description |
|---|---|
| `output_dir` | BEHAV3D output directory. |
| `cell_type` | The **focal** cell type whose response is being measured (typically the target — e.g. `organoid`). |
| `interacting_cell_types` | List of partner cell types whose contact with the focal cells defines an interaction event. |
| `dead_threshold` | Dead-fraction threshold used when scoring focal-cell death within the interaction window. |
| `df_tracks_path` | Optional explicit input CSV override. |
| `show_plots` | Whether to display matplotlib plots interactively (default `True`). |

## Outputs

```
<output_dir>/analysis/<cell_type>/interaction/
    BEHAV3D_<cell_type>_interaction_events.csv
    BEHAV3D_<cell_type>_interaction_response.pdf
```

Plus per-condition breakdowns under `per_sample/`.

## What the future GUI will look like

When the **Analysis → Interaction** sub-tab is implemented, it will expose:

- Focal cell type selector.
- Multi-select for interacting cell types.
- Contact threshold (linked to the Feature Extraction tab's contact threshold).
- Live preview of interaction events on one sample.
- A queue button to add an `INTERACTION_ANALYSIS` step.

This page will be replaced with the full GUI walkthrough once that lands.

## See also

- [Feature Extraction](feature_extraction) — must compute the contact + death features used here.
- [Behavioural Analysis](behavioural) — clusters the partner cell type's behaviour around interactions.
- [Backprojection](backprojection) — visualise interaction events back onto the raw image.
- Source: [behav3d/analysis/interaction_analysis.py:26](../api/behav3d.analysis).
