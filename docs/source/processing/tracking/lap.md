# 🔗 LAP (laptrack)

LAP — Linear Assignment Problem tracking via the [laptrack](https://github.com/yfukai/laptrack) library — links per-timepoint detections by solving a global optimization on the centroid-distance cost matrix (Hungarian-style assignment). It also supports gap closing across short disappearances and optional merge/split events.

## When to use it

- You have moderately dense, actively migrating cells.
- Segmentation occasionally drops a cell for a few frames and you want gap closing to recover the track.
- You may need to model cell **division / fusion** events (off by default — see Merging / Splitting cost below).

## What runs under the hood

For each sample BEHAV3D EXPLORER extracts per-timepoint segment centroids, converts them from pixel coordinates to **physical units** (µm if the metadata's `pixel_distance_xy/z` are in µm), and hands the centroid table to `laptrack`. `laptrack` solves a global assignment problem to link detections across frames.

## Parameters

| Field | Default | Range | Effect |
|---|---|---|---|
| **Track cost (px)** | 45 | 1–999 | Maximum distance a cell can move between two **consecutive** frames to be linked as the same track. |
| **Gap close cost (px)** | 60 | 1–999 | Maximum distance allowed when reconnecting a track across **one or more missing** frames. Should be ≥ Track cost. |
| **Gap close max frames** | 3 | 0–100 | Maximum number of *consecutive* frames a cell can be missing before the gap is considered too large to close. |
| **Merging cost (px)** | 0 | 0–999 | Maximum distance for detecting merge events (two tracks → one). `0` disables merging detection. |
| **Splitting cost (px)** | 0 | 0–999 | Maximum distance for detecting split events (one track → two). `0` disables splitting detection. |

## Step-by-step in the napari plugin

1. **Open the Tracking tab** and pick the **sub-tab** for the cell type you want to track (one sub-tab per cell type — organoids 🟣, immune 🔵, other 🟡).
2. From the **Method** dropdown, pick **LAP (laptrack)**.
3. **Estimate the fastest single-frame displacement in your data**, in the same units as your metadata's `pixel_distance_xy` (usually µm). The easiest way: open the Visualization tab, load one sample, scrub through two consecutive frames, and visually measure how far the fastest cell moves between them.
4. **Set `Track cost` ≈ that displacement**, with a small margin (e.g. fastest cell moves ~30 µm → set Track cost to 35–45). The default 45 µm covers a wide range of typical immune-cell motility.
5. **Set `Gap close cost` ≈ Track cost × 1.2 – 2** (default 60 µm). A cell that disappears for a few frames may have moved further than a single-frame neighbour.
6. **Set `Gap close max frames`** = how many consecutive frames you're willing to bridge. Default 3 is fine for most timelapses.
7. **Leave `Merging cost` and `Splitting cost` at 0** unless you actually expect cell fusion or division — both are very disruptive when triggered by mistake.
8. (Optional) **Use the Apply buttons** to propagate these settings to other immune cell types if they have similar dynamics: **Apply to all Immune** or **Apply to all**.
9. **Click *Run Immune Tracking*** (the green button at the bottom of the sub-tab, named after the cell type) for an immediate run across all samples. To run every cell type at once, or to queue the work, use the tab-level **Run Batch Tracking (All Cell Types)** button or the **+🛒** queue button below the sub-tabs — see [Generic workflow](./index.md#generic-workflow-per-sub-tab).
10. **Inspect the result** in the Visualization tab. If you see fragmented tracks, ID swaps, or cells that should have been bridged across a gap but weren't, see [Tuning the parameters](#tuning-the-parameters).

## Tuning the parameters

After looking at the result in the Visualization tab, identify which failure mode you're seeing and tune accordingly:

**Tracks are fragmented (one biological cell has several short TrackIDs over time)**

- Cause: the cell moves further between consecutive frames than `Track cost` allows.
- Fix: **raise `Track cost`** (e.g. 45 → 60). Re-run and re-inspect.

**Track IDs swap between neighbouring cells**

- Cause: `Track cost` is large enough that two nearby cells' candidate links cross — LAP picks the wrong assignment.
- Fix: **lower `Track cost`**. If lowering it then fragments tracks, the dataset is genuinely ambiguous and you should consider btrack with `Use visual features` for better disambiguation.

**A cell disappears for a few frames and then comes back as a new TrackID**

- Cause: gap closing didn't bridge the disappearance.
- Fix: **raise `Gap close cost`** so the cost of bridging the gap is high enough that LAP picks it over starting a new track. Also check `Gap close max frames` is at least the length of the disappearance.

**False bridges: two unrelated cells are linked across a gap**

- Cause: `Gap close cost` is too generous, or `Gap close max frames` is too high.
- Fix: **lower one of them**. Typically the right move is to lower `Gap close max frames` first (it's easier to interpret).

**One cell becomes two (split detected where there shouldn't be one)** / **Two unrelated cells are merged**

- Cause: `Splitting cost` or `Merging cost` is > 0 and is being triggered by false positives.
- Fix: **set the offending cost to 0** (disabled). Only re-enable if you actually see biological splits / merges and the cost was previously too low to catch them.

**LAP takes too long on dense data**

- Cause: LAP solves a global assignment problem; cost is roughly quadratic in the number of detections per frame.
- Fix: switch to [TrackPy](trackpy) (adaptive search radius) or [btrack](btrack) with `Update method = APPROXIMATE` for very dense data.

## Outputs

Standard BEHAV3D EXPLORER tracking outputs:

- `<output_dir>/images/<sample>/<sample>_<cell_type>_tracked.zarr`
- `<output_dir>/trackdata/<sample>/<cell_type>/<sample>_<cell_type>_tracks.csv`

CSV columns: `TrackID, SegmentID, position_t, position_x/y/z, pixel_position_x/y/z` — see [Tracking overview](index) for the full schema.

## Good practices & tips

- **Start with the defaults and look at the result before tuning.** The defaults (Track cost 45, Gap close 60, Gap frames 3) cover a wide range of typical immune-cell motility. Tune only if you see specific failure modes.
- **Set Track cost to roughly the maximum displacement of the fastest cell between two frames** in your data. Inspect a few frames in napari and measure visually — too small = fragmented tracks; too large = false links across cells.
- **Gap close cost should be ≥ Track cost.** A cell that disappears for *n* frames may have moved further than a single-frame neighbour, so allow more room.
- **Keep merging / splitting at 0 unless you have a biological reason to enable them.** They are off by default because false-positive merges/splits are very disruptive downstream. Turn them on only when you actually see fusion / division events in your data.
- **LAP gets expensive on very dense data** (thousands of cells per timepoint) because it solves a global assignment problem. Consider TrackPy (adaptive search radius) or btrack (approximate update mode) for that case.

## Things this page does **not** claim

- Specific per-sample runtimes — they depend on cell count per frame and total frame count. Profile your own data.

## See also

- [laptrack documentation](https://laptrack.readthedocs.io/) — upstream library, including a precise definition of the cost function.
- [TrackPy](trackpy) — simpler nearest-neighbour alternative.
- [btrack](btrack) — Bayesian alternative with motion model + global hypothesis solver.
