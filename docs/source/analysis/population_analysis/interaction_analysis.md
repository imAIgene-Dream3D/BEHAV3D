# 🤝 Interaction Analysis

Relates target–effector contacts to target death and effector engagement. Unlike Death Dynamics, it **can run without a dead channel** — it simply skips the death-aware panels.

| Button | What it does | Enabled when |
|---|---|---|
| **▶ Run Interaction Analysis (per target)** | For each selected target, analyses its interactions with each selected effector. | ≥ 1 target **and** ≥ 1 effector selected. |
| **▶ Run Interaction Overview** | Produces the cumulative-contacts violin, the cumulative-to-death curves, and the active-killing dashboard across the selected target(s). Works with one or more targets. | ≥ 1 target **and** ≥ 1 effector selected. |

```{note}
When **no dead channel** is configured, Interaction Analysis still runs, but the following are skipped automatically:

- Alive-vs-dead comparison panels (overall and per-sample)
- Fate-based statistics (how many targets survive vs die, and their contact counts)
- Cumulative-to-death curves (Interaction Overview)
- The active-killing dashboard (Interaction Overview)
```

## Interaction settings (two collapsible panels)

Interaction Analysis has **two** collapsible settings panels — one per run — and their controls affect **different plots**. Both are collapsed by default and are saved to `behav3d_parameters.yml` when you run the corresponding step.

**Per-target settings** — apply to **▶ Run Interaction Analysis (per target)**:

| Control | Default | What it does |
|---|---|---|
| **Group overall plots by immune line condition** | off | Splits the two per-target overall plots by the interacting immune cell's line condition (`im_{type}_line_condition`): the cumulative-contact curve draws one line per condition and the alive-vs-dead bar plot colours bars by condition. Only affects the per-target overall plots. |

**Interaction Overview settings** — apply to **▶ Run Interaction Overview** (the violin of cumulative contacts per organoid, the cumulative-to-death curves, and the active-killing dashboard):

| Control | Default | Applies to | What it does |
|---|---|---|---|
| **Before-death window** | 60 min | the **cumulative-to-death curve only** | Look-back window, in minutes before each target's Time of Death. Requires death classification. |
| **Temporal Range** — *All timepoints* / *Custom time range* → Start T – End T | All timepoints | the **cumulative-interactions-per-organoid** and **active-killing** plots | Restricts those plots to a timepoint window. *Custom time range* enables the Start T / End T spinboxes (0-indexed timepoints; defaults 0 and 100). |
| **Annotate by immune line condition** | off | the violin and the active-killing dashboard | Adds the immune line condition as a hue / line-style annotation within the chosen *Group by* grouping, instead of pooling all conditions. Does not affect the cumulative-to-death curve. |
| **Group by** | By target (organoid) type | violin, cumulative-to-death curve, active-killing dashboard | Sets the x-axis grouping. *By target (organoid) type* keeps each organoid type as its own group; *By treatment (immune cell)* pools all organoid types and groups by the interacting immune cell type instead. |

```{note}
The **Before-death window** and the **Temporal Range** are two separate controls for two separate sets of plots: the window only sets how far back before each death the **cumulative-to-death** curve looks, whereas the Temporal Range restricts the timepoints used by the **cumulative-interactions-per-organoid** and **active-killing** plots.
```

## Interaction outputs

| Run | Result PDF |
|---|---|
| Per target ↔ effector | `<output_dir>/analysis/<target>/interaction_analysis/interaction_analysis_<target>_vs_<effector>.pdf` |
| Combined | `<output_dir>/analysis/multi_organoid_comparison/multi_organoid_interaction_comparison.pdf` |

**The per-pair PDF contains:**

- **Cumulative target–effector contacts over time** — overall (all samples) and per sample. This shows how contact accumulates as the experiment runs.
- **If a dead column is present:** **alive-vs-dead comparisons** (overall and per sample) — i.e. do targets that end up dead accumulate more effector contact than those that survive?

**The combined PDF compares across targets** (grouped by target type or by treatment, per the Advanced setting):

- A **violin comparison** of interaction metrics across the groups.
- **If death data is present:** **cumulative-to-death curves** (contact accumulation aligned to each target's time of death).
- **If active-killing data is present:** an **active-killing dashboard** summarising killing-related metrics.

A per-pair statistics CSV is written next to each per-pair PDF in `interaction_analysis/`, and a summary CSV accompanies the combined PDF.

```{note}
The death-aware panels (alive-vs-dead, cumulative-to-death, active-killing dashboard) only appear when the death classification exists. With no dead channel, the PDF simply contains the contact-accumulation plots.
```
