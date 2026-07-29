# 🧹 Filtering

Tab 6 of the BEHAV3D EXPLORER dock widget. The Filtering tab drops or trims tracks that won't survive downstream analysis — short tracks, tracks already dead at the start, tracks past the experiment duration, undersized starting cells — and produces:

1. A **filtered per-timepoint table** that every downstream analysis (Death Dynamics, Interaction, behavioural-state classification, Active Killing) actually reads.
2. A **summarised per-track table** with one row per surviving track. *Currently consumed only by the legacy DTW-based BEHAV3D analysis (see [below](#what-the-summarised-csv-is-for)).*
3. A set of **quality-control PDFs** that let you check at a glance how many tracks each filter removed and what the surviving tracks look like.

```{important}
**Filtering must be run even if you apply no filters at all.** Beyond dropping tracks, this step generates the filtered per-timepoint CSV (and interpolates missing timepoints) that every downstream analysis reads — the pipeline will not proceed to Death Dynamics, Interaction, Invasiveness or Single-Cell classification without it. The downstream state analysis is built to tolerate tracks of uneven length, so the length filters are genuinely optional; running Filtering with everything off still produces the CSV those steps need.
```

![Filtering tab](../_static/screenshots/filtering_tab.png)

```{note}
*Screenshot placeholder.*
```

## Per-cell-type sub-tabs

Filtering is organised as sub-tabs per cell type along the left side, with the usual colour code:

| Icon | Category |
|---|---|
| 🟣 | Organoid |
| 🔵 | Immune |
| 🟡 | Other |

Each sub-tab has its own filter parameters and its own batch run button.

## The five filters

All filters are **independent** — turn any combination on. When a filter is unchecked, its spinner is hidden and the filter is skipped. By default the three time-based filters (experiment duration, min length, max length) are **on**, while *min size at first timepoint* and *dead at first timepoint* are **off**.

### 1 · Trim experiment duration

> ☑ **Trim full time series to max timepoints** — *Max timepoints:* `350`

Drops every row past the chosen experiment duration. Use the **Unit for time-based filters** dropdown at the bottom of the tab to interpret this value as `frames` or `hours`:

- In `frames` mode, this is a 0-based timepoint cutoff applied to `position_t`.
- In `hours` mode, this is applied to the `time` column that Feature Extraction already computed from your metadata's frame interval.

Typical use: cut off the tail of a movie where the dye has bleached or the signal has degraded, so downstream stats don't include the unreliable late timepoints.

### 2 · Filter short tracks

> ☑ **Filter tracks shorter than minimal length** — *Min length:* `30`

Drops tracks whose **time span** (last timepoint − first timepoint) is below the threshold. The unit is the time-unit dropdown (frames or hours), same as filter 1.

Short tracks usually correspond to false detections, transient over-segmentations, or cells that briefly entered the field of view. They bias every per-track statistic, so it's almost always worth filtering them out.

Typical values:

- Immune-cell behaviour: ≥ 30 frames at 1-minute intervals (a 30-minute observation window is the practical minimum for movement / contact statistics).
- Organoids that should persist the whole experiment: a much larger threshold — close to the experiment duration — so you only keep organoids tracked from start to end.

### 3 · Trim long tracks to max length

> ☑ **Trim tracks to maximum length** — *Max length:* `30`
> &nbsp;&nbsp;☑ **Split long tracks into chunks** *(sub-option, on by default)*

For each track longer than the threshold, this controls what happens to the extra frames, via the **Split long tracks into chunks** sub-option:

- **Split on (default):** the track is cut into **consecutive full-length chunks** rather than cropped. The first chunk keeps the original `TrackID`; each following full-length chunk becomes a **new** track (a leftover remainder shorter than the max length is discarded). Every row keeps an `original_TrackID` so backprojection can still paint the real, unsplit track. E.g. a 100-timepoint track at max length 30 yields three 30-frame tracks.
- **Split off:** the track is simply **cropped** to the first *max-length* frames; the rest is dropped.

This is needed for track-level analyses that require uniform-length inputs. The [Track Classification](single_cell/track_classification.md) step offers an equivalent trajectory-size / split control, so this can be done here or there.

### 4 · Filter by minimum size at the first timepoint

> ☑ **Filter by minimal size at first timepoint** — *Min size:* `1000`

Drops whole tracks whose voxel count (`nr_pixels`) at their first appearance (`relative_time = 1`) is smaller than the threshold. If the feature table doesn't have an `nr_pixels` column, the physical **volume** (µm³) is used instead — in that fallback case, toggle the units switch to µm³ so the threshold you enter matches what's being compared.

Common use: removing dust, debris, or over-segmented fragments that briefly look like cells.

### 5 · Filter dead cells at the first timepoint

> ☑ **Filter dead cells at first timepoint**

Drops every track that is already flagged `dead` at its first appearance. This filter is only visible if your metadata has a `dead_channel` and requires that Feature Extraction has already produced a `dead` column.

Typical use: when you want to study **dying** rather than **already-dead** organoids.

## Time unit

> *Unit for time-based filters:* `frames` ▾ / `hours`

Applies to filters 1, 2 and 3 only (experiment duration, min length, max length). Pick whichever matches how you think about your experiment.

## Apply, Run, Queue

Per sub-tab:

- **Apply to all <Category>s** — copy this sub-tab's filter settings to every other sub-tab of the same category.
- **Apply to all** — copy to all cell types.
- **▶ Filter <CellType> Tracks & Summarize** — runs the filter and the summarisation for this cell type only. If a filtered CSV already exists, an **Overwrite Existing Filtered Data?** dialog asks before replacing it.

Main tab:

- **▶ Run Batch Filtering (All Cell Types)** — runs every cell type sequentially. The batch overwrite dialog offers Overwrite All / Skip Existing / Cancel.
- **+🛒** — adds a Filtering step to the [Processing Queue](../plugin_essentials/processing_queue).

Settings are persisted to `behav3d_parameters.yml` under `track_filtering.<cell_type>`.

## Outputs

For each cell type, filtering writes two CSVs and a folder of QC PDFs:

```text
<output_dir>/analysis/<cell_type>/
├── track_features/
│   ├── BEHAV3D_<cell_type>_combined_track_features_filtered.csv      # per-timepoint
│   └── BEHAV3D_<cell_type>_combined_track_features_summarized.csv    # one row per track
└── quality_control/
    ├── BEHAV3D_filter_counts.pdf                       # always
    ├── BEHAV3D_<target>_touching_distribution.pdf      # one per contact target
    ├── BEHAV3D_dead_dye_distribution.pdf               # if a dead dye column exists
    └── BEHAV3D_dead_dye_distribution_t0.pdf            # if a dead dye column exists
```

### The filtered CSV

Same columns as the input Feature Extraction CSV, with only the rows of surviving tracks (subject to experiment duration and max-length trims). One row per `(track, timepoint)`. **This is the file every downstream analysis reads** — Death Dynamics, Interaction, behavioural-state classification and Active Killing all start from here.

```{tip}
If a sample has an **Active Killing** advanced features CSV available, filtering will use *that* as input instead of the plain feature-extraction CSV — so the immune-cell killing columns are preserved through filtering.
```

### What the summarised CSV is for

> One row per `(sample_name, TrackID)`. Built immediately after the filtered CSV, by aggregating each track over time.

The columns are:

| Column | How it's aggregated |
|---|---|
| `sample_name`, `TrackID` | Identity |
| `track_length` | Number of timepoints the track has after filtering |
| `mean_<col>` | **Mean** over the track for every numeric column (e.g. `mean_speed`, `mean_volume`, `mean_directional_persistence`, `mean_<type>_contact`, `mean_percentage_dead_mask`, …) |
| `dies` | `True` if the track is dead at any timepoint (requires the `dead` column) |
| `displacement_from_origin` | **Maximum** displacement from origin reached by the track (i.e. how far it ever wandered from its start) |
| `cumulative_displacement` | Value at the **last** timepoint (total path length) |
| `active_<contact>` | For each contact column, the fraction of *contact timepoints* on which the cell was actively engaged |
| `well`, `exp_nr`, `*_line_condition`, … | First-row metadata, carried over so you can group / compare conditions later |

```{important}
**This file is currently used only by the legacy DTW-based T-cell analysis** (the `BehavioralAnalysisPanel` ipywidgets workflow that runs UMAP + K-Means clustering on per-track means and DTW-distance matrices). The current napari behavioural-state analysis pipeline does **not** use this file — it builds its own rolling-window features from the *filtered per-timepoint* CSV instead.

If you are not running the legacy T-cell analysis you can safely ignore the summarised CSV. It is still generated automatically because it is cheap to compute and a few exploratory notebooks rely on it.
```

### Quality-control PDFs

Always written, one folder per cell type:

| File | Pages |
|---|---|
| `BEHAV3D_filter_counts.pdf` | **(1)** Per-sample histogram of track length **before** filtering. **(2)** Per-sample bar chart showing how many unique tracks survived **each successive filter stage** (the stages that appear depend on which filters you enabled — experiment-duration cut, min size, min length, dead-at-t0). **(3)** Per-sample histogram of track length **after** filtering. |
| `BEHAV3D_<target>_touching_distribution.pdf` | One PDF per `*_contact` column in your features (e.g. `BEHAV3D_organoid_touching_distribution.pdf` for immune-cell tracks). Per-sample bar showing the fraction of timepoints with vs without contact. |
| `BEHAV3D_dead_dye_distribution.pdf` | Violin + strip plot of `mean_dead_dye` across all surviving timepoints. Only generated if the column exists. |
| `BEHAV3D_dead_dye_distribution_t0.pdf` | Same plot but only at the first timepoint of each track — quick check that the dead-at-t0 filter is doing the right thing. |

All PDFs are written at 600 dpi.

```{note}
The QC plots are always generated by the napari Filtering tab — there is no UI toggle to turn them off.
```

## Tips & best practices

- **Always look at `BEHAV3D_filter_counts.pdf` after a run.** The per-stage bar chart is the fastest way to spot a filter that is removing far too many or far too few tracks.
- **Re-run Filtering whenever Feature Extraction parameters change.** The filtered and summarised CSVs become stale otherwise, and downstream analyses will keep reading the old version.
- **The minimum-track-length filter is the single most impactful one** for behavioural analyses. Short tracks add noise to every per-track statistic; ≥ 30 frames at 1-minute intervals is a reasonable minimum for immune cells.
- **For organoids that should persist the whole experiment**, set a *large* min length (close to the experiment duration) rather than a small one — this drops organoids that disappeared mid-run.
- **The time-unit dropdown is per cell type**, but only one unit applies to the three time-based filters of a given sub-tab. Pick the one that matches how you think about that cell type's biology.

```{tip}
**What we filter in practice:**

- **Organoids** are usually tracked accurately and rarely need filtering — at most a small minimum-length filter to drop fragments that fuse at the start of the run.
- **T cells / immune cells:** filter out tracks shorter than roughly **10–30 timepoints**; these short fragments are typically junk.
- **Avoid the legacy fixed-length trim.** Old BEHAV3D set both a minimum *and* a maximum length around 100 so every track was exactly 100 frames long — a workaround for fragmented tracking. With today's accurate tracking that just throws data away; keep long tracks so per-timepoint classification can use them. Reach for *Trim tracks to maximum length* only when you deliberately need every track aligned to the same window.
- **Trim experiment duration** is the right tool when the dead dye bleaches or signal degrades late in the movie — cut off the unreliable tail instead of trimming individual tracks.
```

## See also

- [Feature Extraction](feature_extraction) — produces the input CSV that Filtering reads.
- [Death Dynamics & Interaction](death_dynamics) — consumes the filtered per-timepoint CSV.
- [Single Cell](single_cell/index) — classifies per-cell behavioural states from the filtered CSV.
- [Output Directory & File Layout](../plugin_essentials/output_layout) — where the CSVs and QC PDFs live.
