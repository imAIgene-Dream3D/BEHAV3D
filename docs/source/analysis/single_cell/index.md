# 🧬 Single Cell

The second sub-tab of the **📊 Analysis** tab. Where the [Population Analyses](../population_analysis/index) work at the **population** level, Single Cell classifies the **behaviour of individual cells** over time. It has two inner sub-tabs:

| Inner sub-tab | What it does | Status |
|---|---|---|
| [🔬 State Classification](state_classification.md) | Assigns each cell, at each timepoint, to a recurring **behavioural state** (a motion / intensity / contact "mode") using a Hidden Markov Model. | ✅ Implemented |
| [🛤️ Track Classification](track_classification.md) | Groups whole **trajectories** into clusters by their behaviour over time (DTW), with cluster renaming, a trainable classifier, and exemplar / diagnostic plots. | ✅ Implemented |

## First decide what you are profiling on

Before configuring either sub-tab, settle the question the feature selection answers: **what makes two cells behave differently in your experiment?** The choice determines which feature families you extract and which you feed to the model.

| Profile on | Typical features | Ask this when |
|---|---|---|
| **Motility** | speed, displacement, cumulative displacement, directional persistence, straightness | Cells move, and how they move is the readout |
| **Interaction** | contact with another population, or with their own type; contact duration | Engagement is the readout, whether or not anything follows from it |
| **Morphology** | volume, sphericity, elongation and related shape descriptors | Shape change carries the signal — spreading, polarising, rounding |
| **Channel intensity** | mean or quantile intensity of a given channel, and its fold change | A reporter's level is the readout, including fluctuating reporters such as calcium |

These combine, but starting with a small, biologically interpretable set is better than selecting everything: states are easier to name and the model is easier to trust. Add families once you can explain the states you already have.

## Cell-type scope and selector

At the top of the Single Cell sub-tab there is a single **Cell type** dropdown that both inner sub-tabs share. It lists **immune** and **other** cell types only (multicolor channel splits are excluded) — these are the individually-tracked populations the analysis is designed for. Organoid types are not offered here; they are analysed at the population level in [Population Analysis](../population_analysis/index).

Nothing about these analyses is specific to immune cells. Any population tracked as individual objects can be profiled here — declare it as an *other* cell type if it is not immune, and it becomes available.

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
