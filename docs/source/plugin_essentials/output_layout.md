# Output Directory & File Layout

Every BEHAV3D EXPLORER step writes into a single, canonical folder structure under the **output directory** you pick in the Data Preparation tab. This page is the reference for where things go.

## The output directory

You pick it once, in **Data Preparation → Output Directory**. Every other tab is constructed with a reference to the Data Preparation tab, so they all read and write to the same place.

## Top-level tree

```text
<output_dir>/
├── metadata.csv                          # the metadata table for this output dir
├── behav3d_parameters.yml                # GUI parameter snapshot (auto-saved)
│
├── images/                               # raw images, masks and trained pixel classifiers
│   ├── PixelClassification/              # trained Pixel Classifier / APOC / ConvPaint models
│   └── <sample_name>/                    # one folder per sample: raw zarr, segments, tracked segments, dead mask
│
├── trackdata/                            # per-sample, per-cell-type track tables (CSV)
│   └── <sample_name>/
│       └── <cell_type>/
│
└── analysis/                             # all analytical outputs
    ├── <cell_type>/                                # per-cell-type results
    │   ├── track_features/                         # feature CSVs (raw, filtered, summarised)
    │   ├── quality_control/                        # filtering QC plots
    │   ├── active_killing/                         # immune cell types only — killing events
    │   ├── invasiveness_analysis/                  # immune cell types only — surface-engagement analyses
    │   ├── results/                                # organoid-type dynamics analyses
    │   ├── interaction_analysis/                   # immune ↔ organoid contact analyses
    │   ├── behavioral_states/                      # state classification (.h5ad) + backprojection/ subfolder (state-labelled zarrs)
    │   └── behavorial_trajectories/                 # track classification (DTW) results
    └── multi_organoid_comparison/                  # multi-organoid death-dynamics comparisons
    
```

```{note}
File names use the cell type name from your metadata. If your metadata defines `or_organoid1` and `im_tcell1`, the cell-type names used in paths are `organoid1` and `tcell1` (the `or_` / `im_` / `ot_` prefix is metadata-side only).
```

## Per-step writes

### Data Preparation → Convert to Zarr

| Writes | Path |
|---|---|
| Zarr-converted raw image | `images/{sample}/{sample}.zarr` |

### Segmentation

| Writes | Path |
|---|---|
| Segments | `images/{sample}/{sample}_{cell_type}_segments.zarr` |
| Auxiliary masks  | `images/{sample}/{sample}_*mask*.zarr` |
| Trained pixel-classifier | `images/PixelClassification/` |

### Tracking

| Writes | Path |
|---|---|
| Tracked-segments image | `images/{sample}/{sample}_{cell_type}_tracked.zarr` |
| Tracks table | `trackdata/{sample}/{cell_type}/{sample}_{cell_type}_tracks.csv` |

### Feature Extraction

| Writes | Path |
|---|---|
| Per-timepoint features | `analysis/{cell_type}/track_features/` |
| Extended Analysis → Active Killing (immune only) | `analysis/{immune_type}/active_killing/` |

### Filtering

| Writes | Path |
|---|---|
| Filtered & summarised features | `analysis/{cell_type}/track_features/` |
| Quality-control plots | `analysis/{cell_type}/quality_control/` |

### Analysis & Backprojection

| Writes | Path |
|---|---|
| Death Dynamics results (per target) | `analysis/{cell_type}/results/` |
| Interaction Analysis results (per target) | `analysis/{cell_type}/interaction_analysis/` |
| Behavioural-state classification | `analysis/{cell_type}/behavioral_states/` |
| Cross-cell-type comparisons (combined Death Dynamics / Interaction) | `analysis/multi_organoid_comparison/` |

```{note}
Backprojection exports and Track Classification outputs are also written under each `analysis/{cell_type}/` folder. The easiest way to reopen them is the per-step **👁** buttons or the shared **📄 Results** panel at the bottom of the **Analysis**, **Feature Extraction**, and **Filtering** tabs — a collapsible tree of everything under `analysis/` with **Open in napari**, **Open externally**, and **Reveal in folder** — rather than browsing these paths by hand.
```

## How tabs find prior outputs

The Visualization tab and downstream steps **never re-derive** paths from raw metadata, they always look under `output_dir/...` using the conventions above. As a result:

- If you delete a `*_tracked.zarr`, the Visualization tab will not show tracked segments for that sample, even though the raw metadata is unchanged.
- If you move the output directory between runs, change it in Data Preparation; tabs will then look in the new place.
- The Processing Queue's dependency checker uses this exact path scheme to decide whether a step's inputs are missing.

## Single source of truth

The Data Preparation tab owns:

- The loaded `metadata.csv` (DataFrame).
- The output directory path.
- The combined config (`behav3d_parameters` dict, persisted as YAML in `output_dir/behav3d_parameters.yml`).

## Importing pre-existing data

If you bring data in via the **Import existing** option in the Segmentation or Tracking tab, BEHAV3D EXPLORER copies / validates it into the same canonical paths above. See:

- [Import segmentation](../processing/segmentation/import)
- [Import tracking](../processing/tracking/import)

## Cleaning up

Safe to delete (will just be re-computed next run):

- `analysis/`, analysis outputs only.
- Specific `*_tracked.zarr` if you want to re-track.
- Specific `*_segments.zarr` if you want to re-segment.

Be careful with:

- `images/{sample}/{sample}.zarr` — this is the raw zarr conversion. Deleting it means redoing zarr conversion in Data Preparation.
- `models/*` — trained classifiers / Cellpose models. Keep these unless you intend to retrain.
- `behav3d_parameters.yml` — the saved tab parameters. Deleting resets all tabs to defaults the next time you load metadata.
