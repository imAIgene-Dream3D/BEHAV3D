# 📋 Data Preparation

Tab 1 of the BEHAV3D EXPLORER dock widget. This is where every experiment starts. It owns four artefacts that the rest of the plugin reads:

1. The **output directory**, everything BEHAV3D EXPLORER writes goes under it.
2. The **metadata.csv**, the single source of truth for samples, channels, units and cell types.
3. The **dimension order**, of each raw image (e.g. `TCZYX`).
4. The **zarr-converted raw images**, that downstream tabs read for fast access.

```{image} ../_static/screenshots/data_preparation_tab.png
:alt: Data Preparation tab
:align: center
```

## 1 · Output Directory

```{important}
You **must** set the output directory before any other tab can run.
```

Click **Browse** and pick a folder. Everything BEHAV3D EXPLORER writes goes under it, following the canonical layout documented in [Output Directory & File Layout](../plugin_essentials/output_layout). We recommend short paths and no special characters.

## 2 · Metadata Builder

The Metadata Builder is collapsed by default. Tick the checkbox in the section header to expand it and build a metadata.csv from scratch.

```{note}
**Supported raw image formats**: `.czi`, `.lif` / `.liff`, `.tif` / `.tiff`, `.ims`, and already-converted `.zarr`.
```

Workflow:

1. **Number of samples** — number of samples in the experiment you are going to process/analyze.
2. **Population counts** — how many *organoid* / *immune* / *other* cell types your experiment has. Tick **Include dead channel** if you have it in your samples.
3. Click **Configure Cell Types**. Per-cell-type name fields appear; defaults are `organoid1`, `tcell` (for the first immune), `other1`. Edit to change.
   - For each *immune* type you can tick **Multicolor** + set a channel count to declare a per-channel-split immune cell type (see *Multicolor cell types* below).
4. Click **Create Sample Forms**. One form per sample appears.
5. Fill in the form for sample 1, then click **Fill All from Sample 1** if shared fields (channels, units…) are the same across samples.
6. Click **Save Metadata CSV** and pick a location. The tab then loads it back via the same code path used by section 3.


## 3 · Metadata Loader

```{important}
The loader validates a fixed list of required columns on load: 9 basic columns + one `_line_condition` column per declared cell type. See [Metadata schema reference](#metadata-schema-reference) below for the full list.
```

Click **Browse** and select a `.csv`: the one you have already created, or if you have created it in a previous session, you can skip the Metadata Builder and load the CSV directly here. The metadata loader validates that every required column is present and non-empty for every sample. If validation passes, every other tab (Visualization, Segmentation, Tracking, Feature Extraction, Filtering) receives the loaded metadata.

## 4 · Dimension Order

After metadata loads, the tab populates a table with one row per sample showing the detected dimension and dimension order of the raw image (e.g. `TCZYX`). If a sample reads in the wrong order edit it inline; the corrected order is used by the zarr conversion.


## 5 · Convert to Zarr

Click **Convert to Zarr**:
- The converted raw images are written here:

```
<output_dir>/images/<sample_name>/<sample_name>.zarr
```
- The original raw file is left untouched.

- Optional time-range trimming:

- `t_start`, `t_end`: if your acquisition has dead frames at the start or end, trim them here before conversion.

- The conversion step skips samples whose zarr already exists, unless you tick the "Overwrite existing" checkbox.

---

## Metadata schema reference

BEHAV3D EXPLORER uses a **dynamic, prefixed** schema: instead of hard-coded `tcell_*` / `organoid_*` columns, it discovers cell types from the metadata at load time. The prefix marks the **category** of the cell type:

| Prefix | Category |
|---|---|
| `or_` | Organoid-like (large, roughly spherical, slow-moving) |
| `im_` | Immune (small, fast-moving) |
| `ot_` | Other (any cell that doesn't fit the above) |

Within a category you can have **multiple** named cell types — `or_organoid1`, `or_organoid2`, `im_tcell`, `im_macrophage`, `ot_endothelial`, etc.

### Required columns (every row, every sample)

These 9 columns must be present and non-empty for every sample. The Metadata Loader refuses to proceed otherwise.

| Column | Description | Example |
|---|---|---|
| `sample_name` | Unique sample identifier, used for output paths. | `Sample_A1` |
| `exp_nr` | Experiment number. | `1` |
| `well` | Plate well identifier. | `A1`, `well01` |
| `raw_image_path` | Absolute path to the raw image file (`.czi`, `.lif` / `.liff`, `.tif` / `.tiff`, `.ims`, or already-converted `.zarr`). Forward slashes recommended even on Windows. | `D:/data/exp1/sample_A1.czi` |
| `pixel_distance_xy` | Pixel size in the XY plane. | `0.5` |
| `pixel_distance_z` | Z-step size, same unit as `pixel_distance_xy`. | `2.0` |
| `distance_unit` | Unit of the two above. | `μm` |
| `time_interval` | Time between frames. | `1.0` |
| `time_unit` | `s`, `m`, or `h`. | `m` |

### Optional columns (auto-detected or conditional)

| Column | When you need it | Notes |
|---|---|---|
| `dimension_order` | Optional | Axis order of the raw image (e.g. `TCZYX`, `ZCTYX`). If absent or invalid, the converter auto-detects it from the file format. You normally don't fill it in the Builder, you edit it in Section 4 (Dimension Order) of this tab, only if the auto-detected order is wrong.
| `dead_channel` | Only if you have a dead channel | If the **Include dead channel** option is ticked, this field is mandatory, fill in the channel index (starting from 0) and make sure it matches the correct channel in your raw image. |
| `dead_mask_path` | Only if `dead_channel` is set | Path to a per-sample dead mask. Soft-required: the loader **warns** if missing (asking you to run the dead-mask step), but does not block. |

### Per-cell-type columns

For every cell type you declare, one of the following columns must be present (the prefix is `or_` / `im_` / `ot_` depending on category; `{ct}` is the cell type name):

| Column | Required? | Description |
|---|---|---|
| `{prefix}_{ct}_line_condition` | **Yes** | A free-form label (e.g. `10T_unstim`) used to split visualisations by condition. You enter **Line** and **Condition** as separate fields in the Metadata Builder; the tool merges them into this single column when it saves the CSV. |
| `{prefix}_{ct}_segments_image_path` | Optional | Path to a per-sample segmentation mask. |
| `{prefix}_{ct}_tracks_image_path` | Optional | Path to a per-sample tracked-labels image. |
| `{prefix}_{ct}_tracks_csv_path` | Optional | Path to a per-sample tracks CSV (paired with `_tracks_image_path`). |

```{note}
The three per-cell-type path columns above, together with `dead_mask_path` from the Optional columns table, are **blank on a first run**, and you do not have to fill them yourself: when the Segmentation, Tracking or dead-mask step finishes, the plugin writes the produced file paths back into the metadata CSV for you, so they appear automatically next time you reload.

To plug in pre-existing results from outside BEHAV3D EXPLORER, first write the path in the corresponding field in the Metadata Builder, then use the dedicated **Import** workflow in the [Segmentation](../processing/segmentation/import) or [Tracking](../processing/tracking/import) tab. The external file must be in **`.tif` / `.tiff` format**, the Import widget then validates it, converts it into BEHAV3D EXPLORER's canonical zarr layout, and updates the metadata CSV for you. **Typing a path manually in the CSV is not enough**: every downstream step reads the canonical zarr, not the original file.
```

### Which channel is used for segmentation?

Your raw image often has **several fluorescence channels** (e.g. organoid, T cell, dead). The metadata tells BEHAV3D EXPLORER *which cell types exist* and where the files live, but it does **not** record which channel index to use when you train a pixel classifier, paint labels, or run Cellpose.

That choice is made later, in the **Segmentation** tab of the same plugin, and it depends on the method.

The **only** channel number that belongs in the metadata CSV is `dead_channel` (when you tick **Include dead channel** in the Metadata Builder): the **0-based index** of the dead-stain channel in the raw stack (0 = first channel).

### Multicolor cell types

Use multicolor when the **same biological cell type appears in more than one fluorescent channel** of the same image, for example, T-cell sub-clones each tagged with a different fluorophore (red + green TEGs, RFP + GFP + YFP, etc.). It is also useful when you want to track different fluorescent labels of the *same* population separately and then pool the results.

**How to declare it in the Metadata Builder:**

1. In the cell-type naming section, tick **Multicolor** next to the immune name.
2. Set `N` = number of fluorescent channels for that cell type (typically 2 or 3).
3. Click *Configure Cell Types*. The Builder now creates **`N` separate cell-type entries** in the metadata, named `{base}_1_multicolor`, `{base}_2_multicolor`, …, `{base}_N_multicolor`. Each one gets its own per-cell-type columns (`im_{base}_n_multicolor_line_condition`, `…_segments_image_path`, `…_tracks_image_path`, `…_tracks_csv_path`).
4. In each sample form, fill in **Line** and **Condition** for every multicolor channel, they can be the same (just different colours of the same biology) or different (genuinely different sub-populations).

**What happens downstream:**

- Segmentation and tracking treat each `_multicolor` entry as a normal cell type, one segmentation + one set of tracks per channel.
- A merge step combines all `N` per-channel results into a single output named `{base}_merged` (no `im_`/`or_`/`ot_` prefix). This merged output appears in the CSV after the merge step runs.
- Feature Extraction, Filtering and the rest of the analysis pipeline operate on the merged `{base}_merged` output, not on the individual `_multicolor` entries.

```{tip}
Multicolor is only meaningful when the channels really belong to the *same* cell type. If you have two genuinely different cell types (e.g. T cells + macrophages), declare them as two separate immune cell types instead — that way they keep independent track IDs and features all the way through.
```

### Examples

**A minimal example (no multicolor)** — one organoid type and one T-cell type, each on a single channel.

| sample_name | exp_nr | well | raw_image_path | pixel_distance_xy | pixel_distance_z | distance_unit | time_interval | time_unit | or_organoid_line_condition | im_tcell_line_condition |
|---|---|---|---|---|---|---|---|---|---|---|
| Sample_A1 | 1 | A1 | D:/data/exp1/A1.czi | 0.5 | 2.0 | μm | 1.0 | m | 10T_unstim | TEG_actA |
| Sample_A2 | 1 | A2 | D:/data/exp1/A2.czi | 0.5 | 2.0 | μm | 1.0 | m | 10T_unstim | TEG_actB |

(Optional `or_organoid_segments_image_path` etc. columns omitted because we'll segment from scratch.)

**A multicolor example** — same organoid type, plus T cells split across 2 fluorescent channels (red and green TEGs) with a dead-stain channel. Only the columns that are different from the example above are shown:

| sample_name | … | dead_channel | or_organoid_line_condition | im_tcell_1_multicolor_line_condition | im_tcell_2_multicolor_line_condition |
|---|---|---|---|---|---|
| Sample_A1 | … | 3 | 10T_unstim | TEG_RFP | TEG_GFP |
| Sample_A2 | … | 3 | 10T_unstim | TEG_RFP | TEG_GFP |

- The Builder created **two** T-cell entries: `tcell_1_multicolor` and `tcell_2_multicolor`. Each one is segmented and tracked independently.
- After the merge step, an additional `tcell_merged_segments_image_path` (and `_tracks_image_path`, `_tracks_csv_path`) column will be added automatically. These are the columns the rest of the pipeline consumes.

## Combining multiple experiments

There is no one-click "combine experiments" feature, but you can merge several runs manually and have BEHAV3D treat them as a single batch:

1. **Segment and track each experiment separately first**, with whatever channels and settings each one needs.
2. **Copy each experiment's per-sample `images/<sample>` folders into one shared output directory** so all samples live side by side under the same `images/` folder.
3. **Merge the experiments' `metadata.csv` files into one** (one row per sample, sample names unique) in that shared directory, and **Load Metadata** on it.

BEHAV3D then sees every sample as one experiment, so you can run **Feature Extraction**, **Filtering** and **Analysis** jointly across all of them. Keep sample names unique across experiments to avoid path collisions.

## Tips & troubleshooting

- **Spaces in paths**: avoid them where possible.
- **"File Not Found"** errors during zarr conversion usually mean a typo in `raw_image_path`, open `metadata.csv` and double-check.
- **Channel index off-by-one**: indices are 0-based.
- **`AssertionError: Not all required columns are present`**: the message lists exactly which columns the validator wants. 
## See also

- [Output Directory & File Layout](../plugin_essentials/output_layout) — what gets written where.
- [Visualization tab](../plugin_essentials/visualization) — open one sample after zarr conversion to verify channels.
- [Segmentation overview](../processing/segmentation/index.md) — next step in the pipeline.