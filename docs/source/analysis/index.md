# Analysis

The Analysis section covers everything that happens **after tracking**: turning the tracked-segments + raw image into per-track features, filtering for quality, and running the higher-level analyses (behavioural clustering, death dynamics, interactions, backprojection).

## Pipeline stage in the dock widget

```{mermaid}
flowchart LR
    T["📍 Tracked segments + tracks.csv"] --> FE["🧪 Feature Extraction"]
    FE --> F["🧹 Filtering"]
    F --> A["📊 Analysis<br/>(behavioural, death &amp; morphology, interaction, backprojection)"]
```

| Page | Tab in napari | Implementation status |
|---|---|---|
| [Feature Extraction](feature_extraction) | 🧪 | ✅ Implemented in the GUI |
| [Filtering](filtering) | 🧹 | ✅ Implemented in the GUI |
| [Behavioural Analysis](behavioural) | 📊 (sub) | 🚧 GUI coming soon — function callable from Python |
| [Death & Morphology Analysis](death_morphology) | 📊 (sub) | 🚧 GUI coming soon — function callable from Python |
| [Interaction Analysis](interaction) | 📊 (sub) | 🚧 GUI coming soon |
| [Backprojection](backprojection) | 📊 (sub) | 🚧 GUI coming soon |

```{important}
The **napari Analysis tab is currently a stub** ([behav3d/napari/_stubs.py](../api/behav3d.napari)). The four sub-pages above describe what the future GUI will look like and document the underlying Python functions you can call today from a script / notebook. Each page will get its full GUI walkthrough as the corresponding sub-tab lands.
```

```{toctree}
:hidden:
:maxdepth: 1

feature_extraction
filtering
behavioural
death_morphology
interaction
backprojection
```
