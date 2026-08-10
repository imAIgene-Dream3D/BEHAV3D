# 🫳 Invasiveness Analysis

Where Interaction Analysis looks at *contact*, Invasiveness measures how much of each **immune** cell's surface is engaging an organoid — the same **Organoid Invasiveness** feature computed during [Feature Extraction](../feature_extraction) (percentage of the immune cell's surface within contact distance of an organoid; "invasive" once that reaches ≥ 50 %). This step has its **own** immune-cell picker and target checkboxes, independent of the target/effector selectors above.

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

## Invasiveness outputs

Written under `<output_dir>/analysis/<immune_type>/invasiveness_analysis/`:

| File | Contents |
|---|---|
| `invasiveness_analysis_<immune_type>.pdf` | All the figures: the **fraction-invasive-over-time** curve, the **mean/median surface-contact-% over time** curve, the **per-movie summary** (one dot per movie, using the chosen stat), and — if organoid targets with death data are selected — **fate violins** comparing contact-% / invasive-cell counts on organoids that died vs. survived. |
| `invasiveness_fraction_over_time_<immune_type>.csv` | The per-timepoint fraction of invasive cells (the boolean-≥50 % curve), per sample / target. |
| `invasiveness_perc_over_time_<immune_type>.csv` | The per-timepoint mean/median surface-contact percentage. |
| `invasiveness_per_movie_summary_<immune_type>.csv` | One row per movie: the over-time curve collapsed with the chosen summary stat. |
| `invasiveness_by_fate_<target>_<immune_type>.csv` | Per-organoid contact-%/invasive counts split by whether that organoid died (only when organoid targets with death data are selected). |
