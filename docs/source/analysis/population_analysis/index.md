# 💀 Population Analysis

**Population analysis** is the family of analyses that describe what happens to whole populations and to the contacts between them, as opposed to [Single Cell](../single_cell/index), which classifies the behaviour of individual cells. It answers four related questions about your **target** populations (organoid or "other" cell types) and the **effector** populations that contact them (immune or "other" cell types):

| Page | Question it answers | Where you run it |
|---|---|---|
| [💀 Death Dynamics](death_dynamics) | How, and how fast, does a signal rise across each target population over the movie? | 📊 Analysis tab, Step 1 |
| [🤝 Interaction Analysis](interaction_analysis.md) | How do contacts between targets and effectors relate to that signal, and to effector behaviour? | 📊 Analysis tab, Step 2 |
| [🫳 Invasiveness Analysis](invasiveness) | How deeply does an effector engage the target's surface, from the *effector's* perspective? | 📊 Analysis tab, Step 3 |
| [🎯 Active Killing](active_killing) | Which *individual* effectors are associated with a signal rise in the target they touched? | 🧪 Feature Extraction tab |

```{note}
**Why Active Killing is configured somewhere else.** Conceptually it belongs with the three analyses above — it is about targets, effectors and contact. Mechanically it has to run during [Feature Extraction](../feature_extraction), because it needs the per-timepoint contact and signal columns while they are being computed, and it writes extra columns back into the effector's own feature table. So you configure and run it there, and interpret it here.
```

The first three steps read the **filtered** per-timepoint feature tables produced by the [Filtering](../filtering.md) tab, so you must run Feature Extraction and Filtering for every cell type involved before anything here will turn on.

```{important}
The 📊 Analysis tab currently labels its first sub-tab **💀 Death Dynamics**. The documentation calls this group *Population Analysis*, because the same machinery covers signals other than death and because Active Killing belongs with it conceptually while living in another tab.

**These analyses are not restricted to cell death.** What they actually compute is the rise of a per-object signal — across a population over time (Death Dynamics), and in association with contact (Active Killing, configured in Feature Extraction). Cell death measured with a dye is the original and canonical use, which is why the controls are named after it.

To use them with **any reporter that switches on** — a differentiation marker, an activation reporter, a stress or damage signal — declare that reporter's channel as `dead_channel` in [Data Preparation](../../data_preparation). Everything downstream then works unchanged; read "dead" as "the signal has risen past your threshold".

This applies only to reporters that **go up and stay up**. A **fluctuating** reporter that switches on and off, such as calcium, is not suited to these analyses — use [Single Cell](../single_cell/index) analysis instead, where the reporter intensity becomes a behavioural feature.
```

![Death Dynamics sub-tab](../../_static/screenshots/death_dynamics_tab.png)

```{note}
*Screenshot placeholder.*
```

## Before you start — prerequisites

The buttons in this sub-tab enable themselves only when their inputs already exist on disk. The two selector panels show small status badges so you can see at a glance what is missing.

| Requirement | How to satisfy it | Used by |
|---|---|---|
| A **filtered** features CSV for each target | Run [Filtering](../filtering.md) for that cell type | Everything |
| A **dead channel** declared in the metadata — i.e. the channel carrying whichever rising signal you are measuring | Set `dead_channel` in [Data Preparation](../../data_preparation) | Death Dynamics (required); Interaction (optional, unlocks signal-aware panels) |
| A **`dead` column** in the target's filtered CSV | Run [Feature Extraction](../feature_extraction) with a **Dead mask % threshold > 0**, then re-run Filtering | Death Dynamics; death-aware Interaction panels |
| **Contact columns** for the chosen effectors in the target's CSV | Make sure Feature Extraction was run with the **Contact** family on, for both cell types | Interaction Analysis |

```{important}
Death Dynamics is **completely hidden** until a `dead_channel` is configured in the metadata. If you see only the Interaction step, your dataset has no such channel declared — that is expected. If you do have a rising reporter you want to analyse this way, declare its channel as `dead_channel` and the step will appear.
```

## Selecting cells

### 🎯 Target cells (organoid / other)

A checklist of every **organoid** and **other** cell type detected in the metadata. Each row carries a badge:

- `✅ filtered` — a filtered CSV exists (the checkbox is selectable).
- `⚠ no filtered CSV (run Filtering)` — the checkbox is disabled until you run Filtering for that cell type.
- `✅ dead-col` / `⚠ no dead column` — whether death classification is present in that CSV.

### 🤝 Interaction cells (immune / other)

A checklist of every **immune** and **other** cell type. These are the effectors whose contact with the selected targets defines an interaction — "effector" here simply means the population doing the contacting, whatever its biology. Each effector row shows, per selected target, whether that target's CSV actually contains a contact column for it (`✅` / `⚠`). An effector with no matching contact column for any selected target is disabled.

```{tip}
Select your target(s) **first**. The interaction-cell list only becomes meaningful once at least one target is selected, because availability is judged against the target's contact columns.
```

## Run All Available (Steps 1–3)

The big **▶▶ Run All Available** button runs every step that is currently possible for your selection — per-target and combined Death Dynamics (if a dead column is present) and per-target and combined Interaction Analysis (if effectors are selected) — in one go.

- The **button itself is strict**: it only enables when at least one step can run *right now* with the inputs already on disk.
- The **+🛒 next to it is lenient**: it queues every applicable step. Steps whose inputs are still missing at run time (for example because an earlier Filtering step in the queue hasn't finished) skip themselves with a log message instead of failing the whole queue.

## The Results panel

A shared **Results** panel sits below the sub-tabs (and is visible from Single Cell too). It scans the output directory for result PDFs and lets you re-open any of them in napari without hunting through folders. It refreshes automatically after each run and when you switch sub-tabs; the **DPI** spinner controls the rendering resolution when a PDF is opened in the viewer.

## Tips & best practices

- **Run [Filtering](../filtering.md) for every cell type involved first.** Both steps read the *filtered* CSV; an unfiltered or stale CSV will either disable the buttons or feed downstream stats short / dead-at-start tracks.
- **Match the dead threshold to your data once, in Feature Extraction.** Because the threshold is shared from Feature Extraction, set it carefully there (using the live preview) rather than expecting to adjust it in this tab.
- **Use the combined runs for cross-condition figures.** Per-target runs are best for QC of a single population; the combined PDFs are what you want for comparing organoid types or treatments side by side.
- **Queue the heavy work.** For multi-sample experiments, use the +🛒 buttons to batch Death Dynamics and Interaction Analysis behind your segmentation / tracking / feature steps and let the whole pipeline run unattended.
- **Read the "disappeared" (grey) bar as QC.** A few targets disappearing (usually organoids fusing) is normal; a *large* grey bar is a red flag that segmentation or tracking went wrong for that sample, not a real biological signal.
- **Compare both % dead and absolute dead-pixel counts.** The binary dead threshold tends to over-represent death in *small* organoids (their % crosses sooner), while the absolute dead-pixel count leans toward *big* organoids (more pixels can die). Looking at both — ideally after a size filter so organoids are comparable in size — gives the cleanest picture. The normalised dead-dye increase (scaled to 0 at the start) prevents an already-dead baseline from skewing the trend.

## See also

- [Feature Extraction](../feature_extraction) — computes the death and contact columns these analyses rely on, and owns the Dead mask % threshold.
- [Filtering](../filtering.md) — produces the filtered CSV that both steps read.
- [Single Cell](../single_cell/index) — per-cell behavioural-state and trajectory classification.
- [Output Directory & File Layout](../../plugin_essentials/output_layout) — where the result PDFs and CSVs live.

```{toctree}
:hidden:
:maxdepth: 1

death_dynamics
interaction_analysis
invasiveness
active_killing
```
