# 💀 Death Dynamics & Interaction

The first sub-tab of the **📊 Analysis** tab. It answers two related questions about your **target** cells (organoids, or any "other" cell type) and the **effector** cells that contact them (immune cells, or "other"):

1. **Death Dynamics** — how, and how fast, does each target population die over the course of the movie?
2. **Interaction Analysis** — how do contacts between targets and effectors relate to target death and effector behaviour?

Both steps read the **filtered** per-timepoint feature tables produced by the [Filtering](filtering.md) tab, so you must run Feature Extraction and Filtering for every cell type involved before anything here will turn on.

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
| **▶ Run Combined Interaction Comparison (≥2 targets)** | Compares interaction behaviour across targets. | ≥ 2 targets **and** ≥ 1 effector selected. |

```{note}
When **no dead channel** is configured, Interaction Analysis still runs, but the following are skipped automatically:

- Alive-vs-dead comparison panels (overall and per-sample)
- Fate-based statistics (how many targets survive vs die, and their contact counts)
- Cumulative-to-death curves (Combined run)
- The active-killing dashboard (Combined run)
```

### Advanced settings (Interaction Analysis)

A collapsible section with two controls that apply to the **Combined Interaction Comparison**:

| Control | Default | Meaning | When to change |
|---|---|---|---|
| **Time window before TOD** | 60 min | Look-back window before each target's **Time of Death** used to attribute interactions to that death. | Shorten it if you only want contacts immediately preceding death; lengthen it for slower killing kinetics. |
| **Combined group by** | By target (organoid) type | Whether combined plots are grouped **by target type** or **by treatment** (the effector / immune cell condition). | Switch to *By treatment* when your comparison of interest is across experimental conditions rather than across target types. |

These two settings are saved to `behav3d_parameters.yml` when you run the combined comparison, so they persist between sessions.

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
