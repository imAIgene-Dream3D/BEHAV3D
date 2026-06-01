# Cellpose Training Notebook

When the lab's [pretrained Cellpose models on Zenodo](https://zenodo.org/records/18872978) don't give clean segmentations for your data, the next step is to train your own Cellpose model. BEHAV3D EXPLORER ships a dedicated Jupyter notebook for this:

```
notebooks/train_behav3d_cellpose.ipynb
```

```{important}
This is the **only** notebook in the BEHAV3D EXPLORER workflow that is meant for routine use. Everything else — segmentation, tracking, feature extraction, filtering, analysis — happens inside the napari plugin. The other notebook (`run_behav3d.ipynb` in `notebooks/`) is for **advanced users** who want to run the pipeline from Python or explore the code directly; it is not covered in this wiki.
```

## What it does

The notebook walks you through six steps to produce a custom Cellpose model file that you can later point the [Segmentation → Cellpose](processing/segmentation/cellpose) tab at.

```{mermaid}
flowchart LR
    A["1. Load packages"] --> B["2. Pick model type<br/>(organoid / immune / other)"]
    B --> C["3. Set train + val directories"]
    C --> D["4. Configure channels"]
    D --> E["5. Train Cellpose"]
    E --> F["6. Validate predictions<br/>in napari"]
```

## Requirements

| Requirement | Detail |
|---|---|
| Training set | 3D annotated images, ~70 % of your data, in `training_data/train/`. |
| Validation set | 3D annotated images, ~30 % of your data, in `training_data/val/`. Must be images the model has **never** seen during training. |
| Image format | `png`, `jpg`, `jpeg`, `tif`, `tiff`. |
| Mask pairing | Each image must have a paired mask file in the same folder, same basename but different suffix (e.g. `sample01_img.tif` + `sample01_masks.tif`). |

Expected folder layout:

```text
training_data/
├── train/
│   ├── sample01_img.tif
│   ├── sample01_masks.tif
│   ├── sample02_img.tif
│   └── sample02_masks.tif
└── val/
    ├── sample03_img.tif
    └── sample03_masks.tif
```

```{tip}
You don't need fully annotated 3D volumes. The notebook's `generate_2D_planes` step extracts 2D planes from your 3D volumes for training. In BEHAV3D EXPLORER, the [Cellpose tab](processing/segmentation/cellpose) runs the trained model in **3-D mode** (`do_3D=True`) on the preprocessed zarr volumes.
```

## Step-by-step

### Step 1 — Load packages

Run the first code cell. It imports:

- `behav3d.preprocessing.segmentation.cellpose_training` (the BEHAV3D EXPLORER-side helpers).
- `cellpose.io`, `cellpose.models` (Cellpose itself).
- `napari` (for the visual validation step).

### Step 2 — Pick model type

Use the widget to pick one of:

- **Organoid** — 3D, slowly-changing, often large objects.
- **Immune / Other** — smaller, fast-moving, often round.

This sets the **base Cellpose backbone** for fine-tuning (`cyto3` for organoids, `nuclei` for immune / other) and a default **mask filename** suffix (`_organoids` vs `_immune_other`). Step 6 validation and the [Cellpose tab](processing/segmentation/cellpose) at inference both drop small masks using **minimum object size** — **1000** voxels for organoids, **10** for immune / other (not editable in the napari UI).

### Step 3 — Set training & validation directories

Use the `PathPicker` widgets to point at your `training_data/train/` and `training_data/val/` folders.

### Step 4 — Configure channel labels

For each channel in your **training** images, assign a short label (e.g. `tcell`, `organoid`). These names are only used inside the notebook so you can pick which channels to train on in Step 5. They do **not** have to match metadata cell-type names in EXPLORER.

When you later use the model in the [Cellpose tab](processing/segmentation/cellpose), you separately map each **zarr** channel index to a cell type and choose which cell type to segment. The model must see the same it was trained on — wrong mapping gives bad masks, but the plugin does not enforce label text or channel index parity with the notebook.

```{warning}
If you train on **more than one** input channel, **order matters**: model channel 0 is the first channel you checked in Step 5, channel 1 the second, and so on. At inference, BEHAV3D EXPLORER currently passes **one** mapped channel per run (the cell type you select). Multi-channel models are best trained for a **single** input channel unless you change the training or inference code.
```

### Step 5 — Train

This is the slow step. In the **Train Cellpose** tab you only choose:

- **Which channels** to use as Cellpose inputs (checkboxes).
- **Model basename** (prepended to the saved folder name).

Click **Train the Cellpose model**. The notebook then extracts 512×512 2D crops (`generate_2D_planes`) and calls `train_cellpose` with **fixed** settings (not editable in the UI):

| Setting | Value |
|---|---|
| Epochs | **400** |
| Learning rate | **0.1** |
| Optimizer | SGD, `weight_decay=1e-4` |
| Patch size | **512×512** (from plane extraction) |

To change epochs, learning rate, or Cellpose's internal batch size, edit the notebook cell or `behav3d/preprocessing/segmentation/cellpose_training.py`.

Training needs a **CUDA GPU**; runtime depends on dataset size and GPU (often on the order of tens of minutes to a few hours for the default 400 epochs).

### Step 6 — Validate predictions in napari

The final cell:

1. Loads the trained model.
2. Runs inference on the validation set.
3. Opens napari and overlays predictions on the validation images.

Scrub through samples. Look for:

- Under-segmentation (one mask covers two cells).
- Over-segmentation (one cell split into multiple masks).
- Missing cells.
- Boundary fidelity.

If the validation looks bad, improve the training labels or retrain (edit the epoch count in code if you need more than 400 epochs).

## Output: your custom model

The trained weights are saved under `models/` next to the notebook. The folder name encodes the channels you trained on, for example:

```text
models/cellpose_organoids__channel0-organoid
models/cellpose_immune_other__channel0-tcell-channel1-organoid
```

When you **Load Model** in the Cellpose tab, EXPLORER parses this name and shows a short channel summary as a reminder — it does **not** auto-configure Section 1 or pass multiple channels on **Run**.

## After training: use the model

1. Open napari with the BEHAV3D EXPLORER plugin.
2. Go to **Segmentation → Cellpose**.
3. **Browse** to the trained model directory and **Load Model**.
4. In **Channel Configuration**, map each zarr channel index to the correct cell type (same biological channel you trained on).
5. Select the **cell type to segment** and hit **▶ Run Cellpose** (or queue it).

## Tips

- **Annotate fewer images, well, rather than many images, badly.** A handful of clean masks beats a hundred shaky ones.
- **Validate in Step 6** before running a full dataset in EXPLORER. The default training run is always **400 epochs** unless you change the code.
- **Save intermediate models.** Cellpose writes checkpoints during training — you can keep a better checkpoint instead of the last epoch if needed.
- Training is **2D** (extracted planes); EXPLORER applies the model in **3-D** on zarr volumes. There is no separate “3-D training” path in this notebook.

## See also

- [Cellpose segmentation tab](processing/segmentation/cellpose) — where you'll use the model.
- [Cellpose documentation](https://cellpose.readthedocs.io/) — full reference for training hyperparameters.
- [Pretrained Cellpose models on Zenodo](https://zenodo.org/records/18872978) — try these before training from scratch.
