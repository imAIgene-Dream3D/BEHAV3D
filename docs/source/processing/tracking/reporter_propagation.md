# 🔆 Reporter Propagation

Reporter Propagation is a tracking method for **near-static objects whose
segmentation is unreliable or intermittent over time** — objects that flicker on
and off between timepoints because their fluorescent signal is not always visible
(for example a calcium or other reporter that switches on and off while the cell
itself stays in place).

Instead of linking detections frame to frame, it works across the whole movie at
once:

1. It pools **every** segment detected at **every** timepoint.
2. It groups together any segments that **spatially overlap** each other, regardless
   of how far apart in time they were detected.
3. For each group it takes the **single largest** detected instance and uses that
   one shape as the object's mask, **stamped onto every timepoint** unchanged.

Because the same static mask is reused at every frame, an intensity trace can be
read for each object even at the timepoints where it was not segmentable on its own.

```{important}
Reporter Propagation reuses **one static shape** for every timepoint, so it assumes
each object's shape and position are genuinely static. Any real movement or
morphology change is masked out, because the same shape is pasted into every frame.
```

## Parameters

The Tracking sub-tab exposes two parameters for this method:

| Control | Default | Range | What it does |
|---|---|---|---|
| **Min overlap fraction** | 0.1 | 0.0 – 1.0 | Minimum spatial overlap — as a fraction of the **smaller** segment's area/volume — required for two segments (from any two timepoints, however far apart) to be grouped as the same object. `0` means any shared pixel counts. |
| **Min segment size** | 100 | ≥ 1 (voxels) | Minimum segment size in voxels. Segments smaller than this are ignored entirely (treated as noise), both as candidates for the "largest instance" and when grouping. |

Both descriptions are the definitions shown in the widget's own help text. Raising
**Min overlap fraction** is the control for the case where unrelated nearby objects
are being grouped together as one.

## Step-by-step in the napari plugin

1. **Open the Tracking tab** and pick the sub-tab for the cell type you want to track.
2. From the **Method** dropdown, pick **Reporter Propagation**.
3. Set **Min overlap fraction** and **Min segment size** if the defaults need adjusting.
4. Click **Run *Cell type* Tracking** (the green button at the bottom of the sub-tab),
   or queue it with the tab-level **+🛒** button — see
   [Generic workflow](./index.md#generic-workflow-per-sub-tab).
5. **Inspect the result** in the [Visualization](../../plugin_essentials/visualization) tab.

## Outputs

Standard BEHAV3D EXPLORER tracking outputs:

- `<output_dir>/images/<sample>/<sample>_<cell_type>_tracked.zarr`
- `<output_dir>/trackdata/<sample>/<cell_type>/<sample>_<cell_type>_tracks.csv`

## See also

- [Tracking overview](index).
- [Fragmentation tracking](fragmentation_tracking) — overlap-based tracking for objects
  that are segmented reliably at every timepoint.
