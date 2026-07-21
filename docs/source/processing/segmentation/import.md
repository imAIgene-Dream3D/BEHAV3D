# 📥 Import existing segmentation

The **Import** option in the Segmentation tab lets you bring in pre-existing segmentation masks — typically produced by Imaris, Ilastik, an external Python script — and integrate them into the canonical BEHAV3D EXPLORER output layout so downstream steps (tracking, feature extraction) can consume them.

```{important}
"Import" is **not** a method of *creating* segmentations; it's a method of *validating and copying* existing ones into the expected zarr layout. The actual segmentation must have been done elsewhere.
```

## How it works

The widget reads the **`{prefix}_{ct}_segments_image_path`** columns from your metadata (`or_` for organoids, `im_` for immune cells, `ot_` for other types) and shows one row per sample × cell type, grouped under a `📁 <sample>` header. Each row shows either a status label or an action button:

| Status shown | Meaning | Action |
|---|---|---|
| **No segmentation available** | The metadata column for this sample × cell type is empty. | Fill in the path in Data Preparation → Metadata Builder. |
| **⚠️ File not found** | A path is set but the file does not exist. | Fix the path in the Metadata Builder. |
| **✅ Ready for tracking** | A valid, already-imported `.zarr` (single root array, first chunk dim `1`) whose spatial (`Z,Y,X`) and time (`T`) sizes match the raw image, **and** whose path is already saved in the metadata. | Nothing — already importable as-is. |
| **💾 Save path to metadata** *(button)* | A structurally valid, dimension-matching `.zarr` that the metadata CSV doesn't point at yet (e.g. you just typed/browsed to a new path). | Click to persist the path to the metadata CSV — no file conversion needed. |
| **⚠️ Dimension mismatch** | A structurally valid `.zarr` whose spatial or time size does **not** match the raw image (hover the label for the exact reason). | Re-export with matching dimensions. |
| **🔄 Fix zarr format** *(button)* | A `.zarr` whose first chunk dim is not `1`, or that is a `zarr.Group` instead of a root array. | Click to re-chunk. |
| **🔄 Convert TIFF → zarr** *(button)* | A `.tif` / `.tiff` file (including `.ome.tif`). | Click to convert + re-chunk to zarr. |
| **⚠️ Format not supported (…)** | Any other extension (e.g. `.h5`). | Re-export the segmentation as TIFF first. |

When you click a per-row button, **⚡ Convert / Save all for `<sample>`** (per sample), or **⚡ Convert / Save All Samples** (global), the widget:

1. Reads the source file (TIFF, or an existing zarr that only needs re-chunking).
2. Writes a zarr at the canonical path `<output_dir>/images/<sample>/<sample>_<cell_type>_segments.zarr` chunked one timepoint per chunk, matching every other output in BEHAV3D EXPLORER, and updates the metadata column to point at it.

If every file is already correct, the widget shows **✨ All available segmentation files are in the correct format** instead of a Convert button.

```{important}
The validator checks the **zarr structure** *and* that the segmentation's spatial (`Z,Y,X`) and time (`T`) dimensions match the raw image. It does **not** verify that pixel values are sensible integer labels, or that label IDs are unique per timepoint. Inspect the imported segmentation in the [Visualization tab](../../plugin_essentials/visualization) before running tracking on it.
```

## What "valid" looks like

BEHAV3D EXPLORER expects each segmentation zarr to be:

- A **single root array** (not a zarr Group with sub-arrays).
- A **4-D `(T, Z, Y, X)` integer array** with `0` = background, positive integers = unique cell IDs per timepoint.
- Chunked **one timepoint per chunk** — the first chunk dimension is `1`.

The validator reports the exact reason on failure (chunk mismatch, missing root array, zarr open error, …).

## Importing from common sources

### Imaris

Imaris saves segmentations inside the `.ims` file. The recommended workflow is:

1. Export the segmentation as a `.tiff` time-series from Imaris.
2. Set the `{prefix}_{ct}_segments_image_path` column in `metadata.csv` to the exported `.tiff`.
3. Open the Segmentation tab → Import method → click **⚡ Convert / Fix All Samples**.

### Ilastik / scikit-image / external scripts

Export the result as a TIFF time-series (or a zarr) and point the metadata column at it, then click **Convert**. The import UI only handles TIFF and zarr — for any other format (`.h5`, …) export to TIFF first. The label values are preserved as-is; they only have to be integer-valued (`0` = background, positive integers = cell IDs per timepoint).

## After importing

Once each row shows ✅ Ready for tracking, the Visualization tab can load these segments as Labels layers, and the Tracking tab can run on them just like segmentations produced by APOC / ConvPaint / Pixel Classifier / Cellpose.

## Tips

- **Always preview a converted sample in the Visualization tab before committing the rest.** Channel-order swaps or T ↔ Z swaps in the source file are silent failure modes: the zarr writes fine, but downstream tracking will produce nonsense.
- **Absolute paths are the most robust choice.** Use absolute paths when you set {prefix}_{ct}_segments_image_path for external segmentations (Metadata Builder). They are the most robust choice; relative paths work only if they resolve from the metadata CSV’s directory or the working directory.

## See also

- [Segmentation overview](./index.md) — method comparison and shared concepts.
- [Data Preparation](../../preprocessing/data_preparation) — where the `*_segments_image_path` columns are filled in.
- [Output Directory & File Layout](../../plugin_essentials/output_layout) — the canonical zarr layout.
- [Import existing tracks](../tracking/import) — analogous workflow for tracking data.
