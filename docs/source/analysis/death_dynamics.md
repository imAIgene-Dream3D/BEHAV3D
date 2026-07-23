# 💀 Death Dynamics, Interaction & Invasiveness

The first sub-tab of the **📊 Analysis** tab. It answers three related questions about your **target** cells (organoids, or any "other" cell type) and the **effector** cells that contact them (immune cells, or "other"):

1. **Death Dynamics** — how, and how fast, does each target population die over the course of the movie?
2. **Interaction Analysis** — how do contacts between targets and effectors relate to target death and effector behaviour?
3. **Invasiveness Analysis** — how deeply do the immune cells engage the organoid surface (from the *immune* cell's perspective)?

All three steps read the **filtered** per-timepoint feature tables produced by the [Filtering](filtering.md) tab, so you must run Feature Extraction and Filtering for every cell type involved before anything here will turn on.

![Death Dynamics sub-tab](../_static/screenshots/death_dynamics_tab.png)

```{note}
*Screenshot placeholder.*
```

## Before you start — prerequisites

The buttons in this sub-tab enable themselves only when their inputs already exist on disk. The two selector panels show small status badges so you can see at a glance what is missing.

| Requirement | How to satisfy it | Used by |
|---|---|---|
| A **filtered** features CSV for each target | Run [Filtering](filtering.md) for that cell type | Everything |
| A **dead channel** declared in the metadata | Set `dead_channel` in [Data Preparation](../preprocessing/data_preparation) | Death Dynamics (required); Interaction (optional, unlocks death-aware panels) |
| A **`dead` column** in the target's filtered CSV | Run [Feature Extraction](feature_extraction) with a **Dead mask % threshold > 0**, then re-run Filtering | Death Dynamics; death-aware Interaction panels |
| **Contact columns** for the chosen effectors in the target's CSV | Make sure Feature Extraction was run with the **Contact** family on, for both cell types | Interaction Analysis |

```{important}
Death Dynamics is **completely hidden** until a `dead_channel` is configured in the metadata. If you see only the Interaction step, your dataset has no dead channel — that is expected.
```

## Selecting cells

### 🎯 Target cells (organoid / other)

A checklist of every **organoid** and **other** cell type detected in the metadata. Each row carries a badge:

- `✅ filtered` — a filtered CSV exists (the checkbox is selectable).
- `⚠ no filtered CSV (run Filtering)` — the checkbox is disabled until you run Filtering for that cell type.
- `✅ dead-col` / `⚠ no dead column` — whether death classification is present in that CSV.

### 🤝 Interaction cells (immune / other)

A checklist of every **immune** and **other** cell type. These are the effectors whose contact with the selected targets defines an interaction. Each effector row shows, per selected target, whether that target's CSV actually contains a contact column for it (`✅` / `⚠`). An effector with no matching contact column for any selected target is disabled.

```{tip}
Select your target(s) **first**. The interaction-cell list only becomes meaningful once at least one target is selected, because availability is judged against the target's contact columns.
```

## Step 1 — Death Dynamics

Quantifies the death progression of the selected target population(s) over time, using the sticky `dead` flag and the dead-mask signal computed during Feature Extraction.

The headline curve is the **percentage of dead targets at each timepoint** — the number of targets flagged dead at that timepoint, out of the targets tracked at the first timepoint. The death-signal traces (`percentage_dead_mask` or `nr_dead_mask_pixels`) are **baseline-normalised** — shifted so each starts at 0 at the track's first timepoint — so an already-dead baseline doesn't skew the trend.

| Button | What it does | Enabled when |
|---|---|---|
| **▶ Run Death Dynamics (per target)** | Runs the analysis separately for **each** selected target. | ≥ 1 selected target has a `dead` column. |
| **▶ Run Combined Death Dynamics (≥2 targets)** | Produces a single cross-target comparison. | ≥ 2 selected targets have a `dead` column. |

Next to each button:

- **+🛒** — adds the step to the [Processing Queue](../plugin_essentials/processing_queue) to run later in a batch.
- **👁** — opens the resulting PDF in napari (enabled once it exists; if several targets are selected it offers a chooser menu).

### Death thresholds (read-only)

Below the buttons, a read-only panel lists the **Dead mask % threshold** currently configured for each target. This value is **owned by Feature Extraction** — there is nothing to tune here. To change it, go back to the [Feature Extraction](feature_extraction) tab, set a new Dead mask % threshold, and re-run Feature Extraction (and Filtering) for that cell type.

```{note}
If a selected target shows `⚠ no dead column`, the Death Dynamics buttons stay disabled and a disclaimer tells you to re-run Feature Extraction with a Dead mask % threshold > 0. The death classification is what Death Dynamics measures, so without it there is nothing to plot.
```

### Death Dynamics outputs

| Run | Result PDF |
|---|---|
| Per target | `<output_dir>/analysis/<target>/results/combined_general_<target>_dynamics_analysis.pdf` |
| Combined | `<output_dir>/analysis/multi_organoid_comparison/multi_organoid_death_dynamics_comparison.pdf` |

**The per-target PDF contains:**

- A **line plot of the percentage of dead targets over time**, one line per sample — the headline "how fast does this population die" curve.
- A **stacked bar of the end-of-experiment state** (alive / dead / disappeared) per sample.
- **Per-sample mean ± SEM pages** for the death signals over time (both the raw and a smoothed version), so you can see the average trajectory and its spread per sample.

**The combined PDF adds cross-target context:**

- Percentage-dead-over-time with **one line per sample × target type**.
- An **end-of-experiment stacked bar** (alive / dead / disappeared) for every target × sample.
- **Per-individual-target death-signal traces** (smoothed % dead-mask and smoothed dead-pixel count), and baseline-normalised versions, so single dying targets are visible rather than only averages.

The underlying numbers are written as CSVs alongside the PDFs (overall in `results/`, per-sample in `results/per_sample/`, and a combined table under `multi_organoid_comparison/`). After a successful run a pop-up offers to open the results folder.

## Step 2 — Interaction Analysis

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

### Interaction settings (two collapsible panels)

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

### Interaction outputs

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

## Step 3 — Invasiveness Analysis

Where Interaction Analysis looks at *contact*, Invasiveness measures how much of each **immune** cell's surface is engaging an organoid — the same **Organoid Invasiveness** feature computed during [Feature Extraction](feature_extraction) (percentage of the immune cell's surface within contact distance of an organoid; "invasive" once that reaches ≥ 50 %). This step has its **own** immune-cell picker and target checkboxes, independent of the target/effector selectors above.

```{important}
Invasiveness Analysis only works if you ran Feature Extraction for the immune cell type with the **Organoid Invasiveness** family enabled (which also requires **Contact**). Without the `{target}_invasiveness_perc` columns there is nothing to plot.
```

| Control | Default | Meaning |
|---|---|---|
| **Immune cell type(s)** | — | Which immune types to analyse (checkboxes). |
| **🎯 Targets to compare** | — | Which organoid targets to measure invasiveness against. |
| **Per-movie summary stat** | mean | How to collapse each movie's over-time curve to a single dot per movie: `mean`, `median`, `max`, or `AUC`. |
| **Separate by immune line condition** | OFF | Split results by the immune cell's `line_condition` instead of pooling all conditions. |
| **Timepoint range** | All timepoints | Restrict the analysis to a custom `Start T` – `End T` window. |

The analysis produces three views:

- **Fraction invasive over time** — the % of immune cells flagged invasive (≥ 50 % surface contact) at each timepoint.
- **Mean / median % over time** — the average / typical surface-contact percentage across *all* cells (including non-invasive 0 % ones) at each timepoint.
- **Per-movie summary** — the chosen over-time curve collapsed to one dot per movie using the summary stat.

```{note}
The **Per-movie summary stat** collapses each movie's over-time curve to a single dot per movie:
- **mean** — average of the per-timepoint values.
- **median** — the middle per-timepoint value.
- **max** — the highest per-timepoint value reached.
- **AUC** — the area under the curve normalized by its duration (a time-weighted average, on the same scale as the curve).
```

Buttons mirror the other steps: **▶ Run Invasiveness Analysis**, **+🛒** to queue, and **👁** to reopen the result PDF.

### Invasiveness outputs

Written under `<output_dir>/analysis/<immune_type>/invasiveness_analysis/`:

| File | Contents |
|---|---|
| `invasiveness_analysis_<immune_type>.pdf` | All the figures: the **fraction-invasive-over-time** curve, the **mean/median surface-contact-% over time** curve, the **per-movie summary** (one dot per movie, using the chosen stat), and — if organoid targets with death data are selected — **fate violins** comparing contact-% / invasive-cell counts on organoids that died vs. survived. |
| `invasiveness_fraction_over_time_<immune_type>.csv` | The per-timepoint fraction of invasive cells (the boolean-≥50 % curve), per sample / target. |
| `invasiveness_perc_over_time_<immune_type>.csv` | The per-timepoint mean/median surface-contact percentage. |
| `invasiveness_per_movie_summary_<immune_type>.csv` | One row per movie: the over-time curve collapsed with the chosen summary stat. |
| `invasiveness_by_fate_<target>_<immune_type>.csv` | Per-organoid contact-%/invasive counts split by whether that organoid died (only when organoid targets with death data are selected). |

## Run All Available

The big **▶▶ Run All Available** button runs every step that is currently possible for your selection — per-target and combined Death Dynamics (if a dead column is present) and per-target and combined Interaction Analysis (if effectors are selected) — in one go.

- The **button itself is strict**: it only enables when at least one step can run *right now* with the inputs already on disk.
- The **+🛒 next to it is lenient**: it queues every applicable step. Steps whose inputs are still missing at run time (for example because an earlier Filtering step in the queue hasn't finished) skip themselves with a log message instead of failing the whole queue.

## The Results panel

A shared **Results** panel sits below the sub-tabs (and is visible from Single Cell too). It scans the output directory for result PDFs and lets you re-open any of them in napari without hunting through folders. It refreshes automatically after each run and when you switch sub-tabs; the **DPI** spinner controls the rendering resolution when a PDF is opened in the viewer.

## Tips & best practices

- **Run [Filtering](filtering.md) for every cell type involved first.** Both steps read the *filtered* CSV; an unfiltered or stale CSV will either disable the buttons or feed downstream stats short / dead-at-start tracks.
- **Match the dead threshold to your data once, in Feature Extraction.** Because the threshold is shared from Feature Extraction, set it carefully there (using the live preview) rather than expecting to adjust it in this tab.
- **Use the combined runs for cross-condition figures.** Per-target runs are best for QC of a single population; the combined PDFs are what you want for comparing organoid types or treatments side by side.
- **Queue the heavy work.** For multi-sample experiments, use the +🛒 buttons to batch Death Dynamics and Interaction Analysis behind your segmentation / tracking / feature steps and let the whole pipeline run unattended.
- **Read the "disappeared" (grey) bar as QC.** A few targets disappearing (usually organoids fusing) is normal; a *large* grey bar is a red flag that segmentation or tracking went wrong for that sample, not a real biological signal.
- **Compare both % dead and absolute dead-pixel counts.** The binary dead threshold tends to over-represent death in *small* organoids (their % crosses sooner), while the absolute dead-pixel count leans toward *big* organoids (more pixels can die). Looking at both — ideally after a size filter so organoids are comparable in size — gives the cleanest picture. The normalised dead-dye increase (scaled to 0 at the start) prevents an already-dead baseline from skewing the trend.

## See also

- [Feature Extraction](feature_extraction) — computes the death and contact columns these analyses rely on, and owns the Dead mask % threshold.
- [Filtering](filtering.md) — produces the filtered CSV that both steps read.
- [Single Cell](single_cell/index) — per-cell behavioural-state and trajectory classification.
- [Output Directory & File Layout](../plugin_essentials/output_layout) — where the result PDFs and CSVs live.
