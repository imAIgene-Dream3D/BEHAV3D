# Recommended Settings & Practical Tips

This page collects **hands-on recommendations** for the BEHAV3D EXPLORER napari
workflow — typical parameter values and practical techniques distilled from
training sessions. It complements the per-step method pages (which explain *what*
each setting does) by suggesting *where to start* and *what tends to work in
practice*.

```{note}
These are **recommended starting points, not hard rules**. They may differ from a
form's built-in default, and the best value always depends on your data
(magnification, frame interval, cell type, signal quality). Preview / inspect the
result and adjust. Treat the numbers below as a sensible first guess.
```

## Pipeline order

Work through the dock tabs left to right:

**Data Preparation → Visualization → Segmentation → Tracking → Feature Extraction
(+ Active Killing) → Filtering → Analysis (Death Dynamics & Single-Cell Behaviour).**

All other tabs require metadata to be loaded first, so always start on **Data
Preparation**.

## Data Preparation & metadata

- **Output directory** — point it at the folder that holds your metadata CSV and
  your images directory (your "run" directory). Avoid network/cloud-synced
  locations for large data; processing writes a lot.
- **Number of samples** = number of unique fields of view (ideally one image file
  per FOV).
- **Cell types** are split into **organoid types** and **immune types** — e.g.
  enter `CD4` and `CD8` as two separate immune types. A **multicolor** option lets
  you segment a cell type that is split across several fluorophores separately and
  then merge the channels (e.g. T cells stained with multiple dyes).
- After setting counts, click **Configure cell types**, name them, then **Create
  sample forms**.

Per-sample fields:

| Field | What to enter |
|---|---|
| Name | Experiment / image name (e.g. `Exp1`) |
| Well | e.g. `A4` — used to group FOVs from the same well (can be custom) |
| Image path | Path to the image file for this sample |
| Dim order | Typically `TCZYX` |
| Pixel size | Read it from Zen / IMARIS / ImageJ for your acquisition |
| Time interval | Often ~120 seconds (use your real acquisition interval) |
| Dead channel # | **Zero-based** index (so `0` = the first channel) |

- Use **Fill All from Sample 1** to copy shared values (pixel size, time interval)
  across samples, then manually correct the per-sample fields that actually
  differ (name, image path, etc.).
- **Save** the metadata CSV, then click **Load Metadata** to make the latest
  version active.
- **Convert to Zarr** for fast access (e.g. CZI → Zarr): a worker count around 7
  is a reasonable default; expect roughly 7–10 minutes per image.

## Visualization

- Select a sample to load it into napari. The **Layers** panel shows each channel
  (e.g. T-cell, organoid, dead); toggle visibility with the eye icon and adjust
  colormap / contrast limits per layer.
- Toggle **3D view** with the small square button at the bottom-left of the canvas.
- When **not** in 3D, the bottom slider is **time** and the top slider is **Z**.

## Segmentation (APOC pixel classifier)

Recommended method: **APOC (GPU)** with a probability-map + watershed strategy.

```{note}
APOC is the practical first choice on a normal workstation, especially when you
have many similar datasets (train once, reuse everywhere). If you have a **strong
GPU / HPC**, clean high-resolution images, and accuracy is the priority — and you'd
rather not paint labels — reach for **[Cellpose-SAM](processing/segmentation/cellpose_sam)**
(zero-shot, no training). It is the most accurate but by far the slowest method,
so it fits a modest number of movies. Set up its one-time sidecar environment
first.
```

### Generating training data

- **Examples per sample** defaults to ~3 (first / middle / last time point). Add
  more to capture cells dimming over time or specific events (e.g. touching
  T cells). Counts are zero-indexed (so `0–17` means 18 images).
- After **Generate training data**, work on one channel at a time: deselect the
  layers you're not drawing on, and make sure the target channel/label layer is
  selected (highlighted). Drawing only works in **2D** view.

### Drawing ground-truth labels

- Use **Label 1 = background** and **Label 2 = foreground**. Always draw
  **background first**, then overlay foreground; focus on borders and uncertain
  regions.
- **Organoid technique:** draw background overlapping the background *and the edge*
  of the organoid (defining the border), then draw foreground over the interior.
  Layering it this way gives a very clean organoid border.
- **Touching objects:** draw a thin background line between them, bordered by
  foreground, so the classifier learns to split them into separate objects. For
  organoids, focus this on **time point 1** of each image — propagation tracking
  copies the segmentation forward, so the classifier only has to separate objects
  correctly on the first time point.
- **Elongated T cells:** draw the label slightly **thicker than the signal
  (~3 px)**. Segmentation connects adjacent pixels, so a one-pixel line leaves
  gaps; a little extra thickness yields one clean segment.
- **Bleed-through:** where another channel bleeds into this one (e.g. organoid
  signal in the T-cell channel at higher Z), paint those regions explicitly as
  **background** so the model ignores them — even when they're much brighter.
- **Dead mask:** mostly draw the bright, high-intensity dead "blobs"; including the
  dimmer signal is optional and tends to be noisier.

```{tip}
Whenever you change anything in the training tab, **train the classifier before
moving on** — drawn labels are not necessarily kept if you leave or reload
training data without training first.
```

### Reading the probability map

- Yellow = confident foreground; dark/blue = low confidence. Where an area should
  be foreground but isn't, draw on it with the background brush on the probability
  map and retrain — usually only small tweaks are needed.
- Aim for **background clearly below ~0.3–0.4** and **foreground above ~0.7**, so
  borderline (~0.5) pixels don't flicker over the threshold between slices.
  (Anything with background probability above 0.5 is *not* treated as background.)

### Training parameters (per cell type)

- **Organoids** — input channels: organoid + dead only (leave out the T-cell
  channel; it adds little and slows things down). Feature preset: `medium_preset`
  (sigmas `1, 2, 5, 15`), or `large_preset` (sigmas `1, 2, 5, 10, 25`) for
  multi-scale data.
- **T cells** — pick the channels/features relevant to T cells. Feature preset:
  `small_preset` (sigmas `1, 2, 5`).

Shared random-forest / instance settings:

| Setting | Recommended | Note |
|---|---|---|
| Max depth | 5 | Too low can't learn; too high overfits |
| n_trees | 100 | Number of classifiers that vote on foreground/background |
| Foreground / mask threshold | 0.50 | |
| Seed threshold | 0.80 | Used for splitting touching objects |
| Opening (px) | 0 | Can smooth the outside of organoids |
| Minimum size filter | 0 during training; ~500–1000 px for organoids once trained | |

- Use **Apply Config to All Tabs** after configuring one cell type, then customize
  channels/features for the others.
- Train one cell type (e.g. organoids first) or **Train all classifiers**, then
  review the segmentation and probability map and refine.
- **Show classifier statistics** to see which features mattered; removing
  high-sigma features that rank low speeds up processing.

### Batch segmentation

- Process **all time points**; set workers to what your machine can handle;
  overwrite existing results if re-running.
- A full time series can take on the order of ~1.5 hours per image — run long jobs
  overnight, plugged in.
- Trained classifiers are saved in the `pixel_classification` subfolder and can be
  reused across experiments by copying that folder.

## Tracking

- **Organoids — propagation:** select the organoid layer; the method auto-sets to
  propagation. It copies segmentation frame-to-frame, so there's nothing to tune.
  It's the slower of the two (~5–10 min/sample).
- **T cells — btrack:** select btrack, enable **global track optimization**, and
  set workers around 4 on most laptops (keeps the machine usable). Typical
  ~1–2 min/sample. Useful starting parameters:
  - **Max search radius:** ~100 px (how far ahead to look for the next-frame
    candidate along the predicted trajectory).
  - **Distance threshold:** ~60 px (used to link tracklets during global
    optimization).
- **Run batch tracking (all cell types)** tracks organoids first, then T cells
  (choose "skip existing" if already tracked).
- **Verify it worked:** tracked segments keep the **same colour/ID across time**;
  untracked segments get a new colour/ID every frame. Hover a label to confirm its
  ID stays constant.

## Feature extraction (+ active killing)

- Intensity, contact, and depth are always extracted (fast). **Movement** (fast)
  and **Morphology** (slower) are optional — select everything if you're willing to
  wait; you can prune later.
- **Death threshold** = the fraction of dead-mask pixels inside a cell at which it
  counts as dead. Use **Preview** to see alive (green) vs dead (red) against a time
  point you trust:
  - **Organoids:** ~5 % (0.05) — sounds low, but organoids are large 3-D objects,
    so a few dead cells are a small volume fraction.
  - **T cells:** ~10 % (0.1). A binary is-dead flag is stored alongside the actual
    intensity and percentage, so you can re-tune later.
- **Run Batch Extraction** covers all cell types in one pass (~5–10 min/sample).
  Output lands in `track_data/` as per-track, per-time-point feature tables; files
  can reach ~1 GB on dense data.

**Active killing** — run *after* feature extraction (it relies on extracted
features). It links organoid death to contacting T cells: for a T cell touching an
organoid, it looks ahead by the observation window to see whether that organoid's
dead mask increases.

| Setting | Recommended |
|---|---|
| Immune cell type | `tcell` |
| Observation window | 5 time points (≈10 min) |
| Measure | number of dead-mask pixels |
| Threshold mode | absolute (not multiplier) |
| Absolute threshold | ~30 (dead-mask pixels rising by ~20–30 within the window) |
| Min contact duration | 1 |

Result: a binary active-killing column (1/0) usable in clustering.

```{note}
Caveat: if many T cells contact one organoid at once, all of them may be labelled
active killers even if only one actually kills.
```

## Track filtering

Decide whether to drop or trim tracks before analysis, then run batch filtering.

- **Filter tracks shorter than a minimum length** — removes short, fragmented,
  noisy tracks. In practice: **T cells** ~10–30 time points; **organoids** rarely
  need length filtering (occasionally a small minimum for start-of-run fusing).
- **Trim to a maximum time point** — useful when a dye bleaches late (discard data
  after, say, time point ~600).
- **Maximum-length trimming** is legacy behaviour and generally **not recommended**
  now: with accurate tracking, long tracks should be kept for per-time-point
  classification, so trimming just discards data.
- **Filter dead-at-start cells** — optionally exclude cells already dead at the
  first time point.
- Time-based filters can be set in time points/frames or hours.

## Analysis — death dynamics

- **Requires** organoid segmentation, tracking and feature extraction only — no
  T-cell steps needed.
- Select the organoid type(s) (you can combine multiple) and run. It reuses the
  death threshold from feature extraction.
- **Outputs per sample:** % dead organoids over time, % alive, % disappeared, and
  the normalized dead-dye increase (scaled to 0 at the start so an already-dead
  baseline doesn't skew it). Raw CSVs are saved next to the PDFs for re-plotting.
- **Reading the "disappeared" (grey) bar:** a few organoids disappearing (fusing)
  is normal; a large grey bar suggests something went wrong.
- **Why both % and absolute dead-pixel counts:** the binary % over-represents death
  in small organoids (their threshold is reached sooner), while absolute counts
  lean toward big organoids. Using both — plus a size filter so organoids are
  comparable — gives the cleanest picture.

## Analysis — single-cell behavioural classification

The key idea: behaviour can be classified **per time point**, not just one
super-behaviour per whole track. There are two tabs — state classification and
track classification.

**State classification (per time point)**

- **Feature selection:** cluster on dynamic/morphology features (e.g. elongation,
  extent, solidity, speed, net displacement). Keep **binary** features (dead /
  organoid-contacting / active-killing) *separate* — clustering on them creates
  artificial clusters.
- **Smoothing:** averaged over ~5 time points by default to ride out segmentation
  errors.
- Run it to cluster cells into primary dynamic states, then **rename** clusters by
  meaning (e.g. high speed → scanner; round + static → static; elongated → plastic;
  reusing a name merges clusters). Combine these with the binary groups to get a
  full per-time-point behaviour, then collapse the long list of combinations into
  the behaviour types you care about.
- **Outputs:** state composition over time, per-track behaviour bars, and
  transition matrices (the diagonal dominates because behaviour persists; removing
  self-transitions reveals true cross-class flow).

**Track classification (super-behaviours)**

- Uses dynamic time warping on the sequence of per-time-point states to cluster
  whole tracks into super-behaviours. It clusters on state *ordering*, not raw
  speed.
- A legacy "original BEHAVE" mode (under Advanced Configuration) runs the classic
  DTW clustering on raw features without state classification.
- State and track classifications can be saved as a **`.pkl` pipeline** and
  re-applied to new, similar data for reproducible auto-clustering and renaming.

## Combining multiple experiments

To analyse several runs together: copy the per-sample image folders from each
output folder into one directory and merge their metadata CSVs, then load the
merged folder. You can segment/track each experiment with its own channels/settings
first, then merge and run feature extraction and analysis jointly.

## Performance & CPU workers

- Check your core count (e.g. via system information).
- Use roughly **half to three-quarters** of available cores so the machine stays
  usable while it processes.
- Keep the laptop plugged in and **disable sleep** for long runs — if the machine
  sleeps, processing pauses (locking the screen is fine).

## Settings quick reference

| Step | Setting | Recommended value |
|---|---|---|
| Segmentation | Draw elongated T cells | ~3 px thick |
| Segmentation | Background probability | clearly < 0.3–0.4 (>0.5 = NOT background) |
| Segmentation | Foreground probability | > 0.7 |
| Segmentation | Example images / sample | ~3 |
| Segmentation | Max depth / n_trees | 5 / 100 |
| Segmentation | Foreground / seed threshold | 0.50 / 0.80 |
| Segmentation | Min size (trained) | ~500–1000 px (organoids) |
| Tracking | Organoid method | propagation (~5–10 min/sample) |
| Tracking | T-cell method | btrack (~1–2 min/sample) |
| Tracking | Global track optimization | Enabled |
| Tracking | Workers | ~4 |
| Tracking | Max search radius / distance threshold | ~100 px / ~60 px |
| Feature extraction | Organoid death threshold | ~5 % (0.05) |
| Feature extraction | T-cell death threshold | ~10 % (0.1) |
| Active killing | Immune cell type | `tcell` |
| Active killing | Observation window | 5 time points (≈10 min) |
| Active killing | Threshold mode / value | absolute / ~30 |
| Active killing | Min contact duration | 1 |
| Filtering | T-cell min track length | ~10–30 time points |
| Filtering | Organoid filtering | rarely needed |
| Single-cell | Smoothing window | ~5 time points |
| Data prep | Dim order | `TCZYX` |
| Data prep | Time interval | ~120 s (use your real value) |
| Data prep | Zarr conversion workers | ~7 |
