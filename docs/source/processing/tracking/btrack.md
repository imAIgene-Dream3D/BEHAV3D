# 🧠 btrack (Bayesian)

[btrack](https://btrack.readthedocs.io/) is a Bayesian multi-object tracker. Its
**motion model** is *"used to make forward predictions about the location of objects
using historical information and estimates of the error in the measurements"*, and
its **hypothesis model** is *"used by the global optimizer to build the final set of
tracks"* ([btrack documentation](https://btrack.readthedocs.io/en/stable/user_guide/configuration.html)).

In BEHAV3D EXPLORER, btrack runs on a cell type's per-timepoint segmentation in two
phases: the Kalman-filter tracker (always) and the global hypothesis optimiser
(optional). Per the reference workflow it is the tracker used in routine practice
for immune / motile cells (see [Tracking overview](index)).

## Two-step workflow

btrack runs in two phases inside BEHAV3D EXPLORER:

1. **Step 1 — Kalman filter tracking**. Always runs.
2. **Step 2 — Global hypothesis optimizer**. Optional, off by default. Resolves track merges, splits and false positives over the whole track set.

Tune Step 1 first and **validate visually** before turning on Step 2.

### Step 1 — Kalman Filter Tracking *(always runs)*

| Parameter | Default | Range | Effect |
|---|---|---|---|
| **Config preset** | Cell (bundled) | combo | Picks the motion / hypothesis model. `Cell` and `Particle` are bundled JSON configs shipped with BEHAV3D EXPLORER. `Custom JSON…` enables the **Custom JSON path** field below for your own file. |
| **Custom JSON path** | — | path | Only active when the preset is "Custom JSON…". See `models/README_btrack.md` in the repo for the expected format. |
| **Max search radius (px)** | 100 | 1–9999 | Maximum isotropic search distance for linking objects between frames. Increase for fast cells; decrease to prevent long-range false links. |
| **Update method** | EXACT | EXACT / APPROXIMATE | EXACT = full Bayesian belief matrix (accurate, slower). APPROXIMATE = local spatial search (faster, recommended once you exceed ~1000 cells/frame). |
| **Step size (frames)** | 100 | 1–10000 | Frames processed per tracking iteration — e.g. `100` processes 100 frames' worth of data at a time rather than loading a whole 1000-frame movie into memory at once. Only lower it (e.g. to 50 or 20) if you hit an **out-of-memory error** while tracking a long movie: each chunk then uses less RAM, at the cost of a little extra overhead stitching more chunks together. |
| **Use visual features** | OFF | checkbox | When enabled, raw image intensity statistics (**mean and standard deviation per channel**) of each cell are computed alongside the centroids and fed to the Kalman filter as an extra cue. Requires `raw_image_path` in the metadata. |
| **Workers** | 1 | 1 to (CPU cores − 1) | CPU cores for parallel feature extraction and zarr writing during btrack. The tracker itself is serial. |

### Step 2 — Global Hypothesis Optimizer *(opt-in)*

| Parameter | Default | Effect |
|---|---|---|
| **Enable global track optimization** | OFF | Master switch for Step 2. |
| **Hypotheses** | `P_FP, P_init, P_term, P_link` | Which hypotheses the optimizer considers (definitions below). `P_FP` is always required (the checkbox is disabled in the UI). `P_branch`, `P_dead` and `P_merge` are opt-in. |
| **Distance threshold** (`dist_thresh`) | 60 | BEHAV3D's help text: *"Maximum distance for generating link/branch hypotheses in the optimizer"* (in the unit selected by the µm/px toggle). btrack's own documentation does not define this field in detail. |
| **Time threshold** (`time_thresh`) | 3 | BEHAV3D's help text: *"Maximum frame gap for generating link hypotheses in the optimizer."* btrack's own documentation does not define this field in detail. |

**Hypothesis definitions** (quoted from the [btrack documentation](https://btrack.readthedocs.io/en/stable/user_guide/configuration.html)):

| Hypothesis | btrack definition |
|---|---|
| `P_FP` | "Hypothesis that a tracklet is a false positive detection." |
| `P_init` | "Hypothesis that a tracklet starts at the beginning of the movie or edge of the FOV." |
| `P_term` | "Hypothesis that a tracklet ends at the end of the movie or edge of the FOV." |
| `P_link` | "Hypothesis that two tracklets should be linked together." |
| `P_branch` | "Hypothesis that a tracklet can split onto two daughter tracklets." |
| `P_dead` | "Hypothesis that a tracklet terminates without leaving the FOV." |
| `P_merge` | "Hypothesis that two tracklets merge into one tracklet." |

```{tip}
**Recommended workflow**: tune Step 1 only → preview → fix segmentation issues that surface → tune Step 1 again → only *then* enable Step 2 for refinement. Step 2 makes Step 1's mistakes much harder to diagnose, so debug Step 1 first.
```

## Step-by-step in the napari plugin

1. **Open the Tracking tab** and pick the sub-tab for the cell type you want to track.
2. From the **Method** dropdown, pick **btrack (Bayesian)**.
3. **Configure Step 1 — Kalman Filter Tracking**:
   - **Config preset**: start with **Cell**. Switch to **Particle** only for very small objects (bacteria, organelles). `Custom JSON…` is for advanced users with their own motion model.
   - **Max search radius**: set to roughly the fastest cell's single-frame displacement in your data, in the same units as your metadata's `pixel_distance_xy`. Default 100 µm is generous; lower it for slower cells.
   - **Update method**: leave at **EXACT** for now (more accurate). Switch to **APPROXIMATE** later only if Step 1 is too slow on your data (typically when cells-per-frame exceed a few hundred).
   - **Step size**: default 100 frames is fine for most datasets. Lower it if Step 1 runs out of RAM on long movies.
   - **Use visual features**: tick this if two cells of the same type cross paths and pure motion isn't enough to disambiguate them. Requires `raw_image_path` in your metadata.
4. **Click *Run … Tracking*** (the green button at the bottom of the sub-tab, named after the cell type) with Step 2 still disabled, then **inspect the result** in the Visualization tab.
5. **Iterate Step 1** until the tracks look correct under pure motion-based linking. See [Tuning the parameters](#tuning-the-parameters).
6. **Only now consider enabling Step 2 (Global Hypothesis Optimizer)**:
   - Tick **Enable global track optimization**.
   - Tick the additional hypotheses you need (`P_branch` for cell division, `P_dead` for cell death, `P_merge` for genuine cell-cell fusion). `P_FP` is always on. Don't tick hypotheses you don't biologically expect — false positives are very disruptive.
   - **Distance threshold** and **Time threshold** govern the candidate links / branches the optimizer considers. Defaults (60 / 3) are conservative.
7. **Re-run** and re-inspect. If Step 2 makes things worse, turn it back off and go back to step 5.
8. When happy, queue the work with the tab-level **+🛒** button below the sub-tabs. The queued step is global — it runs every cell type's currently-configured method at run time, using your saved btrack settings (see [Processing Queue](../../plugin_essentials/processing_queue)).

```{tip}
**Always tune Step 1 first, in isolation.** A poor Step 1 + Step 2 produces failures that look like Step 2 problems but are really Step 1 problems. Once Step 1 is correct, Step 2 becomes a refinement step rather than a black box.
```

## Tuning the parameters

After inspecting the result in the Visualization tab, identify the failure mode and tune accordingly.

**Step 1 — Kalman Filter Tracking**

*Tracks are fragmented*
- Cause: `Max search radius` is smaller than the fastest cell's single-frame displacement.
- Fix: **raise `Max search radius`** until the fragmentation goes away.

*Adding an intensity cue*
- **`Use visual features`** computes each cell's per-channel intensity mean and standard deviation and feeds them to the Kalman filter alongside the centroids. Requires `raw_image_path` in the metadata.

*Step 1 is very slow / runs out of memory*
- Cause: large cell counts per frame combined with EXACT update.
- Fix: switch **Update method** to **APPROXIMATE**, and/or lower **Step size** (e.g. 100 → 50) so each iteration uses less RAM.

*The `Cell` preset doesn't suit your data*
- Cause: your cells are unusually small (organelles) or have very specific motion patterns.
- Fix: try `Particle`, or build a custom JSON config — the bundled JSON files in the repo are a useful starting template.

**Step 2 — Global Hypothesis Optimizer**

*Cell divisions are tracked as two unrelated tracks*
- Cause: `P_branch` is not enabled.
- Fix: **tick `P_branch`** in the hypothesis list. Confirm divisions are now tracked as parent → two daughters in the output.

*Spurious branches / merges appear where there shouldn't be any*
- Cause: too many opt-in hypotheses are enabled.
- Fix: **untick the hypotheses you didn't biologically expect**. The defaults (`P_FP, P_init, P_term, P_link`) are deliberately conservative.

*Step 2 doesn't bridge gaps you expected it to*
- Cause: `Time threshold` is too short, or `Distance threshold` is too small.
- Fix: **raise both gradually**. Default `60 / 3` is conservative; doubling either is a reasonable next attempt.

*Step 2 produces too many bridges across unrelated cells*
- Cause: `Time threshold` / `Distance threshold` are too generous.
- Fix: **lower one of them**, usually `Time threshold` first.

*The optimisation step is very slow*
- Cause: a large hypothesis set combined with large `Time threshold` / `Distance threshold`.
- Fix: tighten the thresholds, or run Step 1 alone first to confirm the issue isn't in the input tracks.

## Outputs

Standard BEHAV3D EXPLORER tracking outputs:

- `<output_dir>/images/<sample>/<sample>_<cell_type>_tracked.zarr`
- `<output_dir>/trackdata/<sample>/<cell_type>/<sample>_<cell_type>_tracks.csv`

CSV columns: `TrackID, SegmentID, position_t, position_x/y/z, pixel_position_x/y/z` — see [Tracking overview](index).

## Good practices & tips

- **Start with the `Cell` preset.** It is the bundled default and the one the lab tunes against. Switch to `Particle` only for genuinely small objects (bacteria, organelles).
- **Use `EXACT` update unless tracking is too slow.** EXACT is the safer default; switch to `APPROXIMATE` once your cells-per-frame counts exceed a few hundred and runtime becomes a problem.
- **`Use visual features`** feeds each cell's per-channel intensity mean and standard deviation to the Kalman filter alongside the centroids. Requires `raw_image_path` in the metadata.
- **Enable Step 2 only after Step 1 looks correct.** A poor Step 1 + Step 2 = mysterious failures that are very hard to diagnose. Build confidence in Step 1 first.
- **`P_branch` is the right hypothesis for cell division.** Without it, mitosis events will be tracked as two unrelated tracks starting where the parent ended.
- **Workers are capped at one less than your CPU core count.** btrack itself is serial; the workers parallelise per-timepoint feature extraction and zarr writing, not the tracker.

## Things this page does **not** claim

- Specific per-sample runtimes — depends on cells-per-frame, frames, and Step 2 hypothesis set.

## See also

- [btrack documentation](https://btrack.readthedocs.io/) — upstream library, including precise definitions of every hypothesis and motion-model parameter.
- `models/README_btrack.md` in the BEHAV3D repo — JSON config format and example files.
- [Propagation](propagation) — the overlap-based method for slow, non-dividing organoids.
