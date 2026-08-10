# 📥 Import existing tracks

The **Import tracking** option lets you bring in pre-existing tracking results outside BEHAV3D EXPLORER and integrate them into the canonical BEHAV3D EXPLORER output layout.

```{important}
"Import" is **not** a tracking method per se; it's a method of *validating, converting, and copying* existing track data into the expected zarr + CSV layout so downstream Feature Extraction can consume it.
```

The Import page is shown when you pick **Import tracking** in any cell-type sub-tab of the Tracking tab. It uses the same zarr validation logic as the [segmentation importer](../segmentation/import): a single root array, chunked one timepoint per chunk.

## How it works

The Import page reads the **`{prefix}_{ct}_tracks_image_path`** column from the metadata. For each sample × cell type it shows:

| Status shown | Meaning | Button offered |
|---|---|---|
| *No source path set in metadata* | The metadata column is empty | — (edit metadata first: Data Preparation → Metadata Builder) |
| ⚠️ *File not found* | Path is set but the file doesn't exist | — (fix the path) |
| TIFF source, not yet imported | A `.tif`/`.tiff` source needs conversion + tracks-CSV generation | **🔄 Convert TIFF → zarr** |
| Zarr source, not yet imported | A `.zarr` (or zarr directory) source needs importing + tracks-CSV generation | **📄 Import zarr** |
| ✅ *Already processed* | Both `_tracked.zarr` and `_tracks.csv` already exist at the canonical paths | **🔄 Re-process** (optional, overwrites) |
| ⚠️ *Unsupported format* | The source is neither TIFF nor zarr | — |

When you click the per-row button (**Convert TIFF → zarr**, **Import zarr**, or **Re-process**), or the global **⚡ Process All Samples** button, the widget:

1. Reads the source (TIFF or zarr).
2. For a zarr source it checks the structure (single root array, chunked one timepoint per chunk) and re-chunks it if needed; a TIFF source is converted directly.
3. Writes a zarr at `<output_dir>/images/<sample>/<sample>_<ct>_tracked.zarr` chunked one timepoint per chunk.
4. Computes centroids per timepoint and writes the tracks CSV at `<output_dir>/trackdata/<sample>/<ct>/<sample>_<ct>_tracks.csv` in the canonical schema.

The **⚡ Process All Samples** button only appears once at least one sample has a source path that still needs importing.

If you also have a pre-existing tracks CSV in `{prefix}_{ct}_tracks_csv_path`, the importer will use it as-is rather than recomputing centroids from the labels image — **provided the CSV already follows BEHAV3D EXPLORER's column schema** (`TrackID, SegmentID, position_t, position_x/y/z, pixel_position_x/y/z` — see [Tracking overview](index)). If the schema does not match, recompute by leaving that meadata column empty.

## What "valid" means

A tracked-segments zarr is valid if it satisfies:

- A **single root array** (not a `zarr.Group`).
- A **4-D `(T, Z, Y, X)` integer array** where the pixel value is the **TrackID** (stable across timepoints — same physical cell keeps the same ID).
- Chunked **one timepoint per chunk** (first chunk dim = `1`).

For TIFF input, the structure is assumed to be `(T, Z, Y, X)` with track IDs.

```{important}
As with segmentation import, the validator only checks **zarr structure** — it does not verify the semantics of the labels. Always preview a converted sample in the Visualization tab before running Feature Extraction.
```

## Importing from common sources

### Imaris

Imaris stores tracks within the `.ims` file. Export the tracked labels as a TIFF time-series (with the same integer ID per track across timepoints), set `{prefix}_{ct}_tracks_image_path` in the metadata to that file, and click **Convert TIFF → zarr** in the Import page.

### btrack run outside BEHAV3D EXPLORER

If you ran btrack standalone and wrote a labelled zarr, set the metadata path to it and Import will validate + re-chunk.

### A previous BEHAV3D EXPLORER run on a different output directory

If you tracked once and want to copy the result into a new output directory, set the metadata path to the previous run's `*_tracked.zarr` and Import will copy it into the new output dir's canonical layout.

## Tips

- Always **inspect the converted result** in the Visualization tab before running Feature Extraction.
- **Absolute paths are the most robust choice.**Use absolute paths when you set {prefix}_{ct}_tracks_image_path for external segmentations (Metadata Builder). They are the most robust choice; relative paths work only if they resolve from the metadata CSV’s directory or the working directory.

## See also

- [Data Preparation](../../data_preparation) — where the `*_tracks_image_path` metadata columns live.
- [Output Directory & File Layout](../../plugin_essentials/output_layout) — the canonical paths.
- [Import existing segmentation](../segmentation/import) — analogous workflow at the segmentation stage.
