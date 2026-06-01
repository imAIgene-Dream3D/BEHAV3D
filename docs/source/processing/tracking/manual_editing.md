# ✏️ Manual editing of tracked segments

After automated segmentation + tracking has run, if you find mistakes — split labels, wrong merges, missed cells, junk objects. The Manual Edition workflow lets you fix these interactively inside napari, with undo / redo, time-range scoping and lazy disk persistence.

```{important}
Manual editing operates on **tracked segments** (`*_tracked.zarr`), not on the raw segmentation output (`*_segments.zarr`). Run [Tracking](index) first so each cell has a stable ID over time — then come back here to clean up.
```

## Where to find it

The editor is reached from the **Visualization tab**:

1. Load a sample with **Load Dataset**.
2. The **Manual Edition** group at the bottom of the Visualization tab becomes visible (it only appears when there is at least one `*_tracked.zarr` layer in the viewer).
3. Pick the cell type from the dropdown.
4. Click **✏️ Edit tracked segments**.

The editor panel appears inside the Visualization tab's scroll area.

## The six tools

The toolbar exposes six operations:

| Tool | What it does | Key control |
|---|---|---|
| **Split** | Split one selected label into two or more using seed points placed inside it. | Seed points layer |
| **Merge** | Merge two or more selected labels into one. You pick which ID survives from a **Target ID** dropdown; every other selected label is rewritten to that ID. | Target ID dropdown |
| **Erode** | Erode the selected label with a 3D isotropic ball. | `Radius XY (px)` spinbox (0–50, default 1) |
| **Dilate** | Dilate the selected label with a 3D isotropic ball; expansion never crosses into another label. | `Radius XY (px)` spinbox (0–50, default 1) |
| **Delete** | Erase the selected label(s) across the chosen time range. | Confirmation checkbox "Yes, erase the selected label(s)" |
| **Create** | Drop ≥ 1 seed point(s) on unlabelled background pixels; each seed becomes a new label, propagated forward and backward in time. | Seed points layer |

Every tool has a **👁 Preview … (single timepoint)** button next to the Apply button — it runs the operation on the currently displayed frame only, so you can tune the radius / confirm the selection before committing to the full lifetime.

```{note}
The **Create** tool uses Otsu thresholding + image-guided watershed to grow a new label from each seed point.
```

## Selecting labels

Three ways to pick which label(s) the next operation applies to:

| Input | Effect |
|---|---|
| **Left-click on the layer** | Adds the clicked label to the selection. |
| **Shift-click** | Adds another clicked label to the selection. |
| **Right-click** | Removes the clicked label from the selection. |
| **Type a TrackID in "Add by ID" and press Enter** | Add by numeric ID — handy when you already know it. |
| **Clear** button | Clear the selection. |

The current selection is shown next to **Selected labels:** above the toolbar (the in-UI hint reads: *"left-click adds a label, Shift-click adds another, right-click removes it"*).

## Time range

Each operation runs over a time range:

- **Full lifetime of the selected label(s)** (recommended) — the editor looks up the label's first/last appearance in the tracks CSV and applies the edit across that whole range.
- **Custom range** — pick `Start T` and `End T` manually. Useful for "delete this object only after t=30" type edits.

## Undo / redo

- **↶ Undo** — revert the latest edit (in memory).
- **↷ Redo** — re-apply.

Operations are kept in memory and are not written to disk until you click **Save**. You can stack many edits and then undo back to any point.

## Saving

- **💾 Save tracked segments** — flushes every in-memory edit to the underlying `*_tracked.zarr` and updates the tracks CSV (e.g. shortened lifetimes after a Delete).
- **✗ Discard changes** — drops the in-memory edits.
- **Stop editing** — closes the editor; if there are unsaved edits it first pops a dialog offering **Save and close**, **Discard and close**, or Cancel.

The same Save / Discard / Cancel dialog also appears if you try to switch to another tab with unsaved edits.

## Tips & best practices

- **Track first, then edit**. Editing pre-tracking segmentation gives you back the same broken cell IDs after tracking — pointless.
- **Use full lifetime by default**. Custom ranges introduce errors if you forget that a track might extend further than you remember.
- **Save often** if you've made big edits — the in-memory buffer survives plugin reloads only via Save.
- **Junk objects**: use Delete. If you have many junk objects, it's usually cheaper to redo the segmentation with stricter filtering (e.g. raise the minimum object size in the segmentation method).

## See also

- [Visualization tab](../../plugin_essentials/visualization) — where the editor is launched from.
- [Tracking overview](index) — get tracked segments before editing.
