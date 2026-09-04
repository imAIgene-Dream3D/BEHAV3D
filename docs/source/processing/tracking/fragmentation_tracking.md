# 🧩 Fragmentation propagation

Fragmentation propagation is **not** a centroid-distance method. For each consecutive pair of frames `t → t+1` it propagates the labels of frame `t` into frame `t+1` segmentation mask (watershed-based overlap tracking).

## When it applies

Because linking is by spatial overlap, the method relies on **each object substantially overlapping itself between consecutive frames** — a 2-pixel dilation absorbs small drifts. That single condition, measured in your own data, is what decides whether this method fits. It is not a property of a cell type: the same object overlaps itself at a short frame interval and does not at a long one.

It suits populations that:

- **stay put**, or drift only slightly, so consecutive frames overlap;
- are a **closed set** — no new objects appear part-way through the movie, whether by entering the field or by dividing;
- change mainly by **fragmenting**, for example as a structure breaks up on death.

Multicellular structures such as organoids and spheroids are the usual case (see [Tracking overview](index)). If new objects appear during the movie, or objects no longer overlap themselves between frames, use [btrack](btrack) instead.

## Parameters

The napari plugin **does not expose any tunable parameters** for fragmentation propagation.

## "Track all organoids together" checkbox (organoids only)

By **default**, all organoid cell types are tracked as one cohort. The Tracking tab shows a single combined **🟣 All Organoids** sub-tab (locked to fragmentation propagation) with the checkbox already enabled:

> ☑ Track all organoids together

```{admonition} When to set it
:class: important

- **Leave it ON (default)** to track all organoid types in **one merged watershed with a single, shared set of track IDs** — the combined `all_organoids` output. Per-type tracked files are still written from that same run (split out of the combined result), so each type can still be analysed on its own.
- **Turn it OFF** to track each organoid type **separately** — its own pass, its own independent IDs, and each sub-tab free to pick a **different tracking method**. (The combined tab is locked to Fragmentation Propagation; separate runs share neither a watershed nor an ID space.)
```

When checked, fragmentation propagation runs **once on merged organoid masks** so track IDs are assigned in a single shared scene. This is useful when:

- Multiple organoid types coexist in the same well (e.g. `or_organoidA`, `or_organoidB`) and your main analysis is the combined **`all_organoids`** output (one cohort, one set of paths).
- You want a unified "any organoid" view without running separate tracking jobs per type.

```{tip}
Track structure populations together only when they coexist in the same movie and
should share one overlap-based tracking scene. Keep independent runs when populations
must not influence one another or need different tracking methods. Signal bleed-through
should be resolved or modelled during segmentation; it is not by itself a tracking rule.
```

You **can still analyse each organoid type separately** after a combined run: BEHAV3D also writes the usual per-type tracked zarr and CSV (`<sample>_<organoidA>_tracked.zarr`, etc.), filtered from that same run. Use those files for type-specific feature extraction — use `all_organoids` only when you want everything in one cohort. For **fully independent** tracking (a separate pass per type, no shared ID space), untick the checkbox instead.

**Unticking** rebuilds one independent sub-tab per organoid type. Each sub-tab has its own "Track all organoids together" checkbox (unchecked there) so you can switch back; the choice is saved in `behav3d_parameters.yml`.

## Step-by-step in the napari plugin

1. **Open the Tracking tab.** By default organoids appear as a single **🟣 All Organoids** sub-tab; if you previously unticked "Track all organoids together", pick the individual organoid sub-tab (🟣) instead.
2. **Method is Fragmentation Propagation.** On the combined All Organoids sub-tab the method is fixed; on an individual sub-tab pick it from the **Method** dropdown. Either way the parameter panel shows just one message and the checkbox — there is nothing else to tune from the GUI.
3. **Decide whether to keep `Track all organoids together`:**
   - **On** (default) → one combined run; primary combined output is `all_organoids`, plus per-type tracked files from the same run. Use when you want the unified cohort or will analyse each type from its own `*_tracked.zarr` paths.
   - **Off** → each organoid type is tracked in **its own** pass (isolated masks, no shared ID space). Use when types must not influence each other during tracking (e.g. touching regions in the merged scene).
4. **Run it.** On the combined tab click **▶ Run All Organoids Tracking** (lists the types it covers); on an individual sub-tab click **Run *Organoid* Tracking**. To queue instead, use the tab-level **+🛒** button — see [Generic workflow](./index.md#generic-workflow-per-sub-tab).
5. **Inspect the result** in the Visualization tab. Because the method is overlap-based, its failure modes are different from centroid methods — see [Tuning the result](#tuning-the-result).

## Tuning the result

There are no parameters to change. Tuning it means **improving the segmentation** or **switching to a different tracking method**.

**A track ID switches mid-movie (one organoid suddenly becomes a different ID)**

- Cause: the same organoid was fragmented in the next frame's segmentation — two pieces of foreground instead of one. The watershed step could not decide which fragment inherits the parent ID, so it relabelled.
- Fix: improve the segmentation. Use stricter post-processing (higher `Min size`, higher `EDT` if appropriate, more `Opening`) so the organoid stays as one connected component frame-to-frame, or correct with [Manual Editing](manual_editing) after tracking.

**Tracks are lost as soon as the organoid moves a bit**

- Cause: the method assumes the foreground mask of frame `t+1` overlaps the labels of frame `t` (after a 2-pixel dilation). If your organoids actually drift more than ~2 pixels between frames, the overlap fails.
- Fix: switch to **btrack** (with a small `Max search radius`), which handles organoid-scale motion via a centroid motion model instead of pixel overlap.

**Two objects that briefly touch end up sharing a track ID**

- Cause: when two distinct foreground regions merge in one frame, the watershed assigns the merged blob to one parent ID and that ID propagates from then on.
- Fix: try **Bounded Propagation**, which restricts each ID to one connected region before watershed. If the masks should never have joined, also improve segmentation and verify the selected strategy's EDT direction before changing its threshold.

**You used `all_organoids` but needed fully separate tracking logic**

- Together mode still writes per-type tracked files, but IDs come from one merged run. If types should not share that step (or you need a different method per type), untick **Track all organoids together**, re-run, and track each organoid sub-tab separately.

## Outputs

Standard BEHAV3D EXPLORER tracking outputs:

- `<output_dir>/images/<sample>/<sample>_<cell_type>_tracked.zarr`
- `<output_dir>/trackdata/<sample>/<cell_type>/<sample>_<cell_type>_tracks.csv`

When *Track all organoids together* is on, you get **both**:

- **Combined:** `<sample>_all_organoids_tracked.zarr` and `.../all_organoids/<sample>_all_organoids_tracks.csv` — all organoid types in one cohort.
- **Per type:** the usual paths above for each organoid cell type (same run, split from the combined tracking).

Metadata columns for each organoid type point at the **per-type** tracked files, not `all_organoids`.

## Good practices & tips

- **It assumes frame-to-frame overlap.** If a label in frame `t` no longer overlaps itself in frame `t+1`, the track is lost. Measure that overlap rather than assuming it from the cell type, and use [btrack](btrack) when it fails.
- **Because it is overlap-based, segmentation quality dominates.** If the same organoid is over-segmented into two pieces in one frame, those pieces get different track IDs from then on. Fix the segmentation first, or correct with [Manual Editing](manual_editing) after tracking.
- **The 2-pixel dilation handles small motion drifts** — even slowly-drifting organoids that move a few pixels per frame still propagate cleanly.
- **Choose *Track all organoids together* based on how tracking should run**, not whether you can analyse types separately afterward. Per-type outputs are still written; the difference is **one merged run** (checkbox on) vs **one run per type** (checkbox off). Prefer **off** if organoid types must not affect each other during tracking.

## See also

- [Tracking overview](index).
- [Bounded Propagation](bounded_propagation).
- [Manual Editing](manual_editing).
