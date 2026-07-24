# 🦠 Segmentation

Tab 3 of the BEHAV3D EXPLORER dock widget. Segmentation turns each timepoint of the raw zarr into **instance labels**, one integer label per cell, written to:

```
<output_dir>/images/<sample>/<sample>_<cell_type>_segments.zarr
```

The zarr is a 4-D `(T, Z, Y, X)` `uint16` array; pixel value `0` is background; positive integers are individual cell IDs *within a timepoint* (track IDs are assigned later by the [Tracking](../tracking/index)).

![Segmentation tab](../../_static/screenshots/segmentation_tab_overview.png)

```{tip}
In a hurry? [Troubleshooting](troubleshooting) is a one-page, symptom → fix quick reference for method choice, classifier issues, and instance-processing parameters.
```

## The six methods

A top dropdown lets you pick the segmentation method. The body of the tab swaps to show that method's parameters.

| Method (dropdown label) | What it is | Hardware |
|---|---|---|
| **APOC (GPU)** | GPU pixel classifier — random forest trained on classical image features (Gaussian, edges, Laplacian, …) | OpenCL GPU |
| **ConvPaint (DL pixel classifier)** | Pixel classifier where features come from a pretrained deep network (VGG / DINOv2 / …) and a gradient-boosted classifier predicts every class at once | PyTorch (CUDA GPU recommended; CPU possible) |
| **Pixel Classifier (Random Forest)** | CPU equivalent of APOC — same family of classifier, classical features computed on the CPU | CPU |
| **Cellpose (Deep Learning)** | Pretrained / retrainable deep-learning instance segmenter (cellpose v3) — outputs cell-level masks directly | PyTorch (CUDA GPU strongly recommended) |
| **Cellpose-SAM (zero-shot)** | Foundation-model instance segmenter (cellpose v4) — segments many cell shapes with no training or annotation. Most accurate, most expensive. | PyTorch (good CUDA GPU or HPC) |
| **Import segmentation** | Validates and copies segmentations produced outside BEHAV3D EXPLORER into the canonical layout | n/a |

```{note}
**Two Cellpose entries.** **Cellpose (Deep Learning)** is the classic cellpose v3 backbone you load a pretrained/retrained model into (see [Cellpose](cellpose)). **Cellpose-SAM (zero-shot)** is the newer cellpose v4 foundation model that needs no model file and no training (see [Cellpose-SAM](cellpose_sam)). They run in separate Python environments — Cellpose-SAM sets up a one-time sidecar env on first use.
```

## How to pick a method

Method choice depends mostly on **what your objects look like**, **whether you have a GPU**, and **whether someone has already trained a model for your sample type**:

```{mermaid}
flowchart TD
    A["Your situation?"] --> D["Already segmented elsewhere<br/>(Imaris, Ilastik, script, …)"]
    D --> J["Import segmentation"]

    A --> B["Roundish, well-separated<br/>(organoids, T cells, …)"]
    A --> C["Subtle texture or weak boundaries"]
    A --> E["Other shapes<br/>(elongated, irregular, clustered, …)"]

    B --> G{"Strong GPU / HPC and<br/>accuracy is the priority?"}
    G -- "Yes, zero-shot" --> SAM["Cellpose-SAM<br/>(no training, no labels)"]
    G -- "Pretrained model fits" --> Cellpose["Cellpose v3<br/>(load a model)"]
    G -- "Neither" --> L

    C --> L
    E --> L

    L{"Can you paint a few<br/>foreground/background pixels?"}
    L -- "Yes" --> PC["Train a pixel classifier<br/>(all need painted labels):<br/>APOC · ConvPaint · Pixel Classifier<br/>— pick by GPU; ConvPaint often best<br/>for subtle texture"]
    L -- "No" --> N["Import existing masks, or<br/>Cellpose training notebook"]
```

```{note}
The **Cellpose** methods appear on the *roundish, well-separated* branch because their models target roughly convex cells. **Cellpose-SAM** is zero-shot (works out of the box); classic **Cellpose v3** needs a pretrained or retrained model.
```

Practical rules:


- **APOC** is the right pixel classifier when you have a GPU with working OpenCL drivers and you can paint a few foreground/background pixels per cell type (same labeling step as Pixel Classifier; ConvPaint paints all classes on one layer instead — see [Labeling foreground and background](#labeling-foreground-and-background)). Training is quick; inference is fast on a modern GPU.
- **ConvPaint** uses pretrained deep features (VGG-16 or DINOv2), useful when APOC's classical filter bank cannot separate cells with subtle intensity differences, ambiguous boundaries, or texture-based contrast. Requires a PyTorch-compatible device (CUDA GPU strongly preferred; the widget exposes a CPU fallback).
- **Pixel Classifier (Random Forest)** is the CPU-only fallback, same family of classifier as APOC, but features are computed with `scikit-image` instead of `pyclesperanto`. Use it when no GPU is available.
- **Cellpose** and **Cellpose-SAM** perform **true instance segmentation in one shot** — they directly output separated cell labels for roughly-spherical / convex cells, no pixel labeling and no watershed post-processing. The three pixel classifiers (APOC, ConvPaint, Pixel Classifier) instead output foreground/background or per-class probability, which BEHAV3D EXPLORER turns into instance labels by watershed post-processing.
  - Use **Cellpose (v3)** when the lab already has a pretrained model for your sample type (see [Cellpose page](cellpose) and [Zenodo](https://zenodo.org/records/18872978)), or when you can retrain one on your own masks.
  - Use **Cellpose-SAM (v4, zero-shot)** when you have a strong GPU/HPC and want maximum accuracy with **no training and no annotation** — at the cost of being the slowest, most compute-heavy method (see [Cellpose-SAM page](cellpose_sam)).
- **Import segmentation** is *not* a way to create segmentations, it validates existing masks (.tiff/ .tif) into BEHAV3D EXPLORER's canonical zarr layout so the rest of the pipeline can consume them.

## Common structure of every method page

The three pixel-classifier methods (APOC, ConvPaint, Pixel Classifier) all follow the same three-section UI inside the BEHAV3D EXPLORER dock widget:

1. **Training / configuration** — pick a strategy, generate training data (loads a few timepoints into the napari viewer with empty Labels layers to paint into), train the classifier.
2. **Preview / fine-tune** — run the trained model on the currently displayed timepoint to check quality before committing to a batch run.
3. **Batch segmentation** — process either all timepoints or a custom range across every sample in the metadata. See [Batch segmentation controls](#batch-segmentation-controls) below.

The per-method pages below document the exact widgets, parameters, and tips for each one. The labeling step, the post-processing parameters, and the batch controls are shared across the pixel-classifier methods and are described once in the three sections below.

(labeling-foreground-and-background)=
## Labeling foreground and background

The three pixel-classifier methods (APOC, Pixel Classifier, ConvPaint) are trained by painting a few labels into the napari viewer on the Labels layer(s) created by **Generate Training Data**. How the labels are organised differs:

- **APOC and Pixel Classifier are binary, one Labels layer per cell type.** You paint two classes with these brush values (the legend is shown in the widget):

| Brush value | Meaning | Colour in legend |
|---|---|---|
| `0` | Eraser | — |
| `1` | Background | red |
| `2` | Foreground (the cell type you are training) | cyan |

- **ConvPaint is multi-class on a single layer.** All classes are painted on one `User Provided Labels` layer using the index shown in its **Legend** tab (`1` = background, `2` = first cell type, `3` = next, …).

**Labeling tips (apply to all three methods):**

- You only need **enough pixels to disambiguate the classes** — typically dozens to a few hundred per class per timepoint. You do *not* need to paint whole cells.
- **Bias your strokes toward boundaries** and any visually-confusing region (touching cells, low-contrast areas) — that is where pixel classifiers fail most.
- Label examples from **different timepoints and Z-slices** (e.g. start / middle / end and a few Z levels) so the classifier sees the variability in your data.
- **Balance examples across all classes** — don't only label the easiest cells.
- **Avoid contradictory labels**; only label regions you are confident about.
- **Save your labels often** (the *Save Labels* / *Save User Labels* button) — unsaved labels are lost when the viewer closes.

### Drawing technique

These practical conventions come from the training sessions and give noticeably cleaner segmentations:

- **Draw background first, then overlay foreground.** Lay down the background class along and just over the object edge, then paint the foreground over the part that sits on the object. Layering it this way produces a very clean border.
- **Split touching objects with a thin background line.** Where two objects touch, draw a thin background line between them, bordered by foreground, so the classifier learns to separate them. For organoids tracked by **propagation**, you only need to get this right on **timepoint 1** of each sample — propagation copies the segmentation forward.
- **Draw thin / elongated cells thicker than the signal (~3 px).** Segmentation connects adjacent pixels, so a one-pixel-wide stroke on an elongated T cell leaves gaps; a slightly over-thick line yields one clean segment.
- **Paint channel bleed-through as background.** Where another channel's signal bleeds into the one you are training (e.g. organoid signal appearing in the T-cell channel at higher Z), explicitly paint it as background — even when it is much brighter — so the model learns to ignore it.
- **For the dead mask, label the bright dead blobs.** Concentrate on the clear high-intensity dead signal; including very dim red is optional and tends to be noisier.

### Reading the probability map

After training, the probability map shows the classifier's confidence (e.g. yellow = confident foreground, dark/blue = low confidence). Use it to target corrections:

- Aim for **background clearly below ~0.3–0.4** and **foreground above ~0.7**. Anything above **0.5 is not treated as background**, so keep borderline pixels away from that line to stop them flickering across slices.
- Where an area *should* be foreground but reads low, paint a little more **foreground** there (or **background** where it wrongly reads high) and **retrain** — only small tweaks are usually needed.

(instance-post-processing-parameters)=
## Instance post-processing parameters

APOC, ConvPaint and the Pixel Classifier are *pixel* classifiers: their raw output is a foreground/background mask or a per-class probability map. BEHAV3D EXPLORER turns that into **instance labels** with a watershed-based post-processing pipeline. The same parameters appear (with method-specific defaults) wherever a watershed strategy is selected; their meaning is identical everywhere:

| Parameter | Meaning |
|---|---|
| **EDT threshold** | Erosion-style threshold on the Euclidean Distance Transform of the mask: voxels with EDT ≥ this value become watershed seeds. **Higher → seeds shrink toward each object's core, breaking the thin neck between touching objects → more splitting**; **lower → neighbouring cores stay connected → objects stay merged**. `0` disables splitting. |
| **Opening px** | Morphological opening applied to the mask before watershed. Smooths boundaries and removes small protrusions / speckles. `0` disables. |
| **Min size** | Segments smaller than this many voxels are discarded after watershed (noise / debris filter). |
| **Fill holes** | Fill internal gaps inside objects before watershed (recommended for solid 3-D objects). |
| **Mask threshold** *(probability strategies)* | Foreground cutoff applied to the probability map. |
| **Seed threshold** *(probability strategies)* | Higher cutoff used to place watershed seeds; should be ≥ Mask threshold. **Higher → keeps only each object's confident core as a seed, splitting more** touching objects; lower merges neighbouring cores together, splitting fewer. Same direction as EDT threshold above. |

### Instance post-processing tuning

- Touching cells **merged** into one label → **raise EDT threshold** (Mask + EDT/Watershed strategy), or **raise Seed threshold** (Probability Map strategy).
- One cell **split** into several labels → **lower EDT threshold**, or **lower Seed threshold**.
- Background labelled as cells → **raise Mask threshold**, or add more background labels.
- Cells lost at the edges → **lower Mask threshold**.
- Small speckles appear as tiny labels → **raise Min size**; real small cells dropped → **lower Min size**.
- Fuzzy / jagged boundaries or single stray voxels → **raise Opening**; thin structures eroded → **lower Opening** (or `0`).

(advice-large-vs-small-objects)=
### Advice: large vs small objects

For **APOC**, **ConvPaint**, and the **Pixel Classifier**: tune **each cell type on its own tab**. Settings for **big** objects (e.g. organoids) are usually wrong for **small, crowded** ones (e.g. immune cells / T cells).

- Work **one cell type at a time**: train (or paint), run a **preview**, look at the labels, adjust sliders, repeat. Settings on one tab do not apply to another.
- **Big** (e.g. organoids): use a larger minimum size than for single cells. Stuck together → **raise** EDT for Mask + EDT (the direction differs for Peak EDT). Fuzzy outline → try a **probability-map** strategy if your method has one.
- **Small / crowded** (e.g. immune, T cells): **small** minimum size; **lower** EDT to start. One cell split into many → **lower** EDT further.
- **Calculate before copying:** estimate a single cell's diameter from the image and convert it with the current XY/Z spacing. For Minimum size, a tolerant first preview is about half the estimated object volume. For EDT, use the object's pixel diameter and the active strategy, then adjust by preview. Historical experiment values are useful only as named examples from acquisition-matched data.
- **Several cell types in one experiment?** In APOC or ConvPaint, **Advanced (per cell type)** lets each type use a different strategy; still tune post-processing separately on each tab.

```{note}
Which sliders you see depends on the strategy you pick. Default numbers differ by method — see each method page. Treat built-in values as UI defaults, not biological recommendations. Always adjust for the current resolution and the **nature of your data** (big vs small, crowded vs sparse, sharp vs fuzzy boundaries).
```

(batch-segmentation-controls)=
## Batch segmentation controls

The pixel-classifier methods (APOC, ConvPaint, Pixel Classifier) share the same batch-run controls:

| Control | Default | Effect |
|---|---|---|
| **Process ALL timepoints** | ON | Segment the full movie of every sample. Uncheck to enable the From–To range. |
| **From t / to t** (range) | — | Inclusive timepoint range, active only when *Process ALL* is off. |
| **Workers** | 1 | Parallel CPU workers for zarr I/O, capped at one less than your CPU core count. For the GPU methods (APOC/ConvPaint) the GPU is the actual segmenter, so 1–2 is usually optimal. |
| **Overwrite existing results** | OFF | Re-segment samples that already have an output zarr. Otherwise they are skipped. |
| **▶ Run …** | — | Runs immediately and blocks the GUI until done. |
| **+🛒** | — | Queues the step in the [Processing Queue](../../plugin_essentials/processing_queue) to run later, sequentially with the rest of the pipeline. |

```{note}
Cellpose has its own Run / **+🛒** controls in its Section 3, and Import uses Convert buttons instead — see their pages.
```

## Outputs

Every method writes its instance labels to the same canonical path:

```
<output_dir>/images/<sample>/<sample>_<cell_type>_segments.zarr
```

Trained pixel-classifier models (APOC, ConvPaint, Pixel Classifier) are written under `<output_dir>/images/PixelClassification/`. Cellpose models are *not* stored under the output directory, they are loaded from wherever you point the Browse field.

See the [Output Directory & File Layout](../../plugin_essentials/output_layout) page for the full canonical layout.

```{toctree}
:hidden:
:maxdepth: 1

apoc
convpaint
pixel_classifier
cellpose
cellpose_sam
import
troubleshooting
```
