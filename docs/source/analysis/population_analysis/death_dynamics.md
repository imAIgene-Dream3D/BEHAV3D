# 💀 Death Dynamics

Quantifies how the measured signal progresses across the selected target population(s) over time, using the sticky `dead` flag and the dead-mask signal computed during Feature Extraction. If your `dead_channel` carries a reporter other than a death dye, read every "dead" below as "past the threshold you set".

The headline curve is the **percentage of dead targets at each timepoint** — the number of targets flagged dead at that timepoint, out of the targets tracked at the first timepoint. The death-signal traces (`percentage_dead_mask` or `nr_dead_mask_pixels`) are **baseline-normalised** — shifted so each starts at 0 at the track's first timepoint — so an already-dead baseline doesn't skew the trend.

| Button | What it does | Enabled when |
|---|---|---|
| **▶ Run Death Dynamics (per target)** | Runs the analysis separately for **each** selected target. | ≥ 1 selected target has a `dead` column. |
| **▶ Run Combined Death Dynamics (≥2 targets)** | Produces a single cross-target comparison. | ≥ 2 selected targets have a `dead` column. |

Next to each button:

- **+🛒** — adds the step to the [Processing Queue](../../plugin_essentials/processing_queue) to run later in a batch.
- **👁** — opens the resulting PDF in napari (enabled once it exists; if several targets are selected it offers a chooser menu).

## Death thresholds (read-only)

Below the buttons, a read-only panel lists the **Dead mask % threshold** currently configured for each target. This value is **owned by Feature Extraction** — there is nothing to tune here. To change it, go back to the [Feature Extraction](../feature_extraction) tab, set a new Dead mask % threshold, and re-run Feature Extraction (and Filtering) for that cell type.

```{note}
If a selected target shows `⚠ no dead column`, the Death Dynamics buttons stay disabled and a disclaimer tells you to re-run Feature Extraction with a Dead mask % threshold > 0. The death classification is what Death Dynamics measures, so without it there is nothing to plot.
```

## Death Dynamics outputs

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
