# 🧬 Single Cell

The second sub-tab of the **📊 Analysis** tab. Where [Death Dynamics & Interaction](../death_dynamics) work at the **population** level, Single Cell classifies the **behaviour of individual cells** over time. It has two inner sub-tabs:

| Inner sub-tab | What it does | Status |
|---|---|---|
| [🔬 State Classification](state_classification.md) | Assigns each cell, at each timepoint, to a recurring **behavioural state** (a motion / intensity / contact "mode") using a Hidden Markov Model. | ✅ Implemented |
| [🛤️ Track Classification](track_classification.md) | Groups whole **trajectories** into clusters by their behaviour over time (DTW), with cluster renaming, a trainable classifier, and exemplar / diagnostic plots. | ✅ Implemented |

## Cell-type scope and selector

At the top of the Single Cell sub-tab there is a single **Cell type** dropdown that both inner sub-tabs share. It lists **immune** and **other** cell types only (multicolor channel splits are excluded) — these are the motile cells whose behaviour the analysis is designed for. Organoid types are not offered here; they are analysed at the population level in [Death Dynamics](../death_dynamics).

Whatever cell type you pick applies to whichever inner sub-tab you are on. A status line next to the dropdown tells you how many cell types are available, and the inner sub-tabs re-check the output directory whenever you change the selection so their buttons reflect what already exists on disk.

```{important}
Single Cell reads the **per-track features CSV** produced by [Feature Extraction](../feature_extraction) for the chosen cell type. If you haven't run Feature Extraction (and usually [Filtering](../filtering)) for that cell type yet, the feature checklists will be empty and the run buttons will stay disabled.
```

```{toctree}
:hidden:
:maxdepth: 1

state_classification
track_classification
```
