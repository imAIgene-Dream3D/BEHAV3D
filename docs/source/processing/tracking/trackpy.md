# 🎯 TrackPy

[TrackPy](http://soft-matter.github.io/trackpy/) is a particle-tracking library originally written for colloidal-particle tracking (Crocker-Grier algorithm). In BEHAV3D EXPLORER it provides a **simpler and faster alternative to LAP** for sparser data — when cells are well-separated and you don't need built-in merge / split modelling.

## When to use it

- Cells are **sparse** — most candidate links are unambiguous nearest neighbours.
- Movement is stochastic or moderately directional — without strong persistence that would require a motion model.
- LAP is slow or over-linking on your dataset.

## What runs under the hood

For each sample BEHAV3D EXPLORER extracts per-timepoint segment centroids in **physical units** (µm if the metadata uses µm) and hands them to TrackPy's linker.

## Parameters

| Field | Default | Range | Effect |
|---|---|---|---|
| **Search range (px)** | 31 | 1–999 | Maximum centroid distance between consecutive frames considered for linking. Set roughly to the displacement of the fastest cell between frames. |
| **Memory (frames)** | 2 | 0–100 | Number of frames a cell can disappear and still be reconnected to its previous track — analogous to LAP's gap-closing max frames. |
| **Adaptive stop** | 10.0 | 0–100 | Minimum search range allowed during *adaptive* search. When TrackPy hits a subnet that is too large, it shrinks the search radius iteratively; it stops shrinking once the radius reaches this value. Higher = TrackPy stops shrinking sooner, so it gives up less search radius before falling back. |
| **Adaptive step** | 0.95 | 0.01–1.0 | Multiplier applied to the search range at each adaptive-search iteration. Closer to 1.0 → slow, gentle reduction; closer to 0.5 → aggressive. |

## Step-by-step in the napari plugin

1. **Open the Tracking tab** and pick the sub-tab for the cell type you want to track.
2. From the **Method** dropdown, pick **TrackPy**.
3. **Estimate the fastest single-frame displacement in your data** (same logic as LAP — load a sample in Visualization, scrub two consecutive frames, eyeball the fastest cell's step).
4. **Set `Search range`** ≈ that displacement, in the same units as your metadata's `pixel_distance_xy` (default 31 px).
5. **Set `Memory`** = how many consecutive missing frames a cell can have and still be reconnected (default 2). This is the analogue of LAP's gap closing.
6. **Leave `Adaptive stop` and `Adaptive step` at their defaults** (10.0 and 0.95) unless you specifically see TrackPy taking very long on a dense sample.
7. (Optional) **Use Apply buttons** to copy these settings to other sub-tabs.
8. **Click *Run … Tracking*** (the green button at the bottom of the sub-tab, named after the cell type) for an immediate run across all samples. To run every cell type at once, or to queue the work, use the tab-level **Run Batch Tracking (All Cell Types)** or **+🛒** controls below the sub-tabs — see [Generic workflow](./index.md#generic-workflow-per-sub-tab).
9. **Inspect the result** in the Visualization tab; if you see fragmented tracks, ID swaps, or "SubnetOversizeException" warnings, see [Tuning the parameters](#tuning-the-parameters).

## Tuning the parameters

**Tracks are fragmented (one cell has several short TrackIDs over time)**

- Cause: the cell moves further between frames than `Search range` allows.
- Fix: **raise `Search range`** to roughly the fastest cell's single-frame displacement.

**Track IDs swap between neighbouring cells**

- Cause: `Search range` is so large that two close cells become indistinguishable nearest neighbours.
- Fix: **lower `Search range`**. If that fragments other tracks, the dataset is genuinely ambiguous — consider LAP (global optimisation) or btrack (motion model).

**A cell flashes off for one or two frames and comes back as a new TrackID**

- Cause: `Memory` is too low.
- Fix: **raise `Memory`** to at least the longest expected disappearance (e.g. 2 → 5). Don't set it super high — large memory values invite spurious bridges between unrelated cells.

**"SubnetOversizeException" warnings in the log**

- Cause: a frame had too many candidate links inside a single search radius — TrackPy is forced to fall back to its slower recursive linking strategy.
- Tracking still finishes correctly, but each occurrence slows the run down.
- Fix: **lower `Search range`** to the smallest value that still produces continuous tracks for the fastest cells. If you can't lower it without fragmenting tracks, consider switching to LAP or btrack.

**TrackPy is extremely slow on a particular sample**

- Cause: adaptive search keeps shrinking iteratively because subnets are too large.
- Fix: **lower `Adaptive step`** (e.g. 0.95 → 0.7 or 0.5) so each iteration shrinks the search radius more aggressively. The trade-off is occasional slightly worse links — TrackPy's adaptive logic is conservative by default.
- **Raising `Adaptive stop`** also helps: TrackPy stops shrinking sooner once the radius hits this floor, so it spends less time recursing.

**Cells fuse or divide in your dataset**

- TrackPy has no native merge / split detection. Switch to LAP (with `Merging cost` / `Splitting cost` > 0) or btrack (with `P_branch` / `P_merge` hypotheses).

## Outputs

Standard BEHAV3D EXPLORER tracking outputs:

- `<output_dir>/images/<sample>/<sample>_<cell_type>_tracked.zarr`
- `<output_dir>/trackdata/<sample>/<cell_type>/<sample>_<cell_type>_tracks.csv`

CSV columns: `TrackID, SegmentID, position_t, position_x/y/z, pixel_position_x/y/z` — see [Tracking overview](index).

## Good practices & tips

- **Adaptive search is what keeps TrackPy fast on dense data.** Leave `Adaptive stop` and `Adaptive step` at their defaults unless tracking is failing badly. Lowering `Adaptive step` (e.g. 0.5) is the right move only when you see TrackPy taking very long on a single sample.
- **Set Search range = the maximum reasonable single-frame displacement** in your data. Inspect the raw movie in napari to estimate — too small = fragmented tracks; too large = ID swaps between neighbouring cells.
- **`Memory` recovers brief segmentation gaps**, similar to LAP's gap closing. Start at 2; raise it if you see TrackID flipping when a cell briefly disappears.
- **No native merge / split detection.** If cells fuse or divide in your dataset, prefer [LAP](lap) (with merging / splitting cost > 0) or [btrack](btrack).
- **If you see "SubnetOversizeException" warnings in the log**, BEHAV3D EXPLORER has already retried with the slower recursive strategy — tracking still finishes but the run is slower. Consider lowering Search range if this happens often.

## Things this page does **not** claim

- That TrackPy is faster than LAP in every case — for very sparse data it generally is, but for moderately dense data, LAP's global optimisation can produce better tracks at comparable speed. Measure on your own data.

## See also

- [TrackPy documentation](http://soft-matter.github.io/trackpy/) — upstream library, including a detailed explanation of the adaptive-search algorithm.
- [LAP](lap) for dense data and merge/split modelling.
- [btrack](btrack) for a Bayesian motion model.
