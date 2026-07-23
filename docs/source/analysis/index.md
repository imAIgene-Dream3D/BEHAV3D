# Analysis

The Analysis section covers everything that happens **after tracking**: turning the tracked segments + raw image into per-cell features, filtering for quality, and running the higher-level analyses (death dynamics, cell–cell interactions, behavioural-state classification) and finally visualising those results back on the images.

In the BEHAV3D EXPLORER dock widget these steps are spread across four tabs:

```{mermaid}
flowchart LR
    T["📍 Tracked segments"] --> FE["🧪 Feature Extraction"]
    FE --> F["🧹 Filtering"]
    F --> A["📊 Analysis"]
```

| Page | Tab in napari | Status |
|---|---|---|
| [Feature Extraction](feature_extraction) | 🧪 Feature Extraction | ✅ Documented |
| [Filtering](filtering.md) | 🧹 Filtering | ✅ Documented |
| [Death Dynamics, Interaction & Invasiveness](death_dynamics) | 📊 Analysis → 💀 Death Dynamics | ✅ Documented |
| [Single Cell](single_cell/index) | 📊 Analysis → 🧬 Single Cell | ✅ State Classification · ✅ Track Classification |

```{note}
The **📊 Analysis** tab has two sub-tabs:

- **💀 Death Dynamics** — population death dynamics, target–effector interaction analysis, and immune-cell invasiveness analysis (organoid/other vs immune/other).
- **🧬 Single Cell** — per-cell behavioural classification: **🔬 State Classification** (per-timepoint HMM states) and **🛤️ Track Classification** (whole-trajectory DTW clustering with a trainable classifier).

**Backprojection** — painting state / track-cluster labels back onto the raw images — is  the final step inside each Single Cell workflow (**State Classification → Step 4** and **Track Classification → Step 5**).
```

```{toctree}
:hidden:
:maxdepth: 2

feature_extraction
filtering
death_dynamics
single_cell/index
```
