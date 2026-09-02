# Analysis

The Analysis section covers everything that happens **after tracking**: turning the tracked segments + raw image into per-cell features, filtering for quality, and running the higher-level analyses (population signal dynamics, cell–cell interactions, behavioural-state classification) and finally visualising those results back on the images.

```{note}
**Which route fits depends on your readout, not on your cell types.**

- A reporter that **switches on and stays on** — a death dye, a differentiation or activation marker — belongs in [Population Analysis](population_analysis/index), which measures how that signal rises across a population and whether the rise follows contact.
- A **fluctuating** reporter that goes on and off, such as calcium, belongs in [Single Cell](single_cell/index), where its intensity becomes a behavioural feature.
- **No reporter at all** is fine: [Single Cell](single_cell/index) profiles cells on motility, morphology, or contact alone.

Nothing here requires an immune-versus-target design. A single population with one readout is a valid experiment.
```

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
| [Population Analysis](population_analysis/index) — Death Dynamics, Interaction, Invasiveness, Active Killing | 📊 Analysis → 💀 Death Dynamics (Active Killing runs in 🧪 Feature Extraction) | ✅ Documented |
| [Single Cell](single_cell/index) | 📊 Analysis → 🧬 Single Cell | ✅ State Classification · ✅ Track Classification |

```{note}
The **📊 Analysis** tab has two sub-tabs:

- **💀 Death Dynamics** — the population analyses: signal dynamics across target populations, target–effector interaction, and effector invasiveness (organoid/other vs immune/other). **Active Killing** completes this group but is configured in the 🧪 Feature Extraction tab.
- **🧬 Single Cell** — per-cell behavioural classification: **🔬 State Classification** (per-timepoint HMM states) and **🛤️ Track Classification** (whole-trajectory DTW clustering with a trainable classifier).

**Backprojection** — painting state / track-cluster labels back onto the raw images — is  the final step inside each Single Cell workflow (**State Classification → Step 4** and **Track Classification → Step 5**).
```

```{toctree}
:hidden:
:maxdepth: 2

feature_extraction
filtering
population_analysis/index
single_cell/index
```
