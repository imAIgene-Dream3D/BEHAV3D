# ✨ Cellpose-SAM (zero-shot)

Cellpose-SAM is the **foundation-model instance segmenter** in BEHAV3D EXPLORER. Like classic [Cellpose](cellpose) it outputs cell-level instance labels directly from raw intensities (no pixel labeling, no watershed post-processing), but it is built on the newer **Cellpose-SAM** backbone ([cellpose v4](https://github.com/MouseLand/cellpose)) and is **zero-shot** — it segments a wide range of cell shapes with *no training and no annotation*.

For method choice across all six segmentation options (including when Cellpose-SAM is the recommended choice versus APOC / ConvPaint / classic Cellpose), see the [segmentation overview](./index.md#how-to-pick-a-method).

```{important}
Cellpose-SAM segments **one cell type per run** (the channel(s) for that cell type must be selected). T cells and organoids cannot be predicted in a single pass — either run the tool once per cell type, or tick **Run all cell types in one batch** to loop through every declared type automatically.
```

## Section 1 · Cellpose-SAM environment (one-time setup)

Cellpose-SAM needs **cellpose ≥ 4**, which cannot share a Python environment with the pinned **cellpose 3** that the rest of BEHAV3D EXPLORER uses (your existing fine-tuned classic-Cellpose models depend on v3). To keep both, BEHAV3D EXPLORER creates a small **side ("sidecar") environment** containing cellpose v4.

| Control | Effect |
|---|---|
| **Status label** | Reports the current state. Before setup it prompts you to install; after setup it shows a green **Ready** status with the detected cellpose / torch version and GPU. |
| **Set up Cellpose-SAM environment** | Creates the sidecar env. This is a **one-time step per computer** — click it once before first use. |

The sidecar **reuses this environment's existing PyTorch install**, so it only adds a few MB. Your classic cellpose v3 install and trained models are left completely untouched.

## Section 2 · What to segment

| Control | Default | Effect |
|---|---|---|
| **Cell type** | first detected | The cell type to segment. Determines the output filename and which raw channel(s) feed the model. |
| **Channel Selection** | — | A per-cell-type checkbox picker choosing which raw channel(s) feed Cellpose-SAM for the chosen cell type. This is Cellpose-SAM's **own** channel config (`cellpose_sam.channel_selection`), separate from classic Cellpose's per-channel panel. |
| **Run all cell types in one batch** | OFF | Segment every detected cell type in one go instead of only the selected one. When ticked, the Cell type dropdown is disabled. |
| **Process all timepoints** | ON | Segment the whole movie. Uncheck to expose the From/To range. |
| **From / To** | 0 / 0 | 0-indexed, **inclusive** timepoint range, active only when *Process all timepoints* is off. Timepoints outside the range keep whatever they already contain — a partial run never clears other timepoints. |

## Section 3 · Device

| Control | Default | Effect |
|---|---|---|
| **GPU Device** | auto | Which CUDA device to run on. Leave on the default/auto GPU. |
| **Force CPU-only processing** | OFF | Hides all GPUs from Cellpose-SAM. Much slower — only force CPU if the GPU is busy with something else. |

```{note}
A full GPU run already saturates the GPU **and** every CPU core (see **CPU threads / Cooldown** in the Advanced block). CPU-only mode is dramatically slower for 3-D volumes.
```

## Section 4 · Segmentation settings

A shared **px / µm unit toggle** at the top of this section governs every size- or radius-like field below (diameter, sharpen/smooth radius, tile-norm blocksize, flow3D smooth, and the Section 5 size filter). Values are stored internally in pixels/voxels regardless of what unit is displayed; µm values are converted using each sample's own pixel size from the metadata.

| Parameter | Default | Range | What it does |
|---|---|---|---|
| **Model** | `cpsam` | `cpsam`, `cpsam_v2` | Which shipped Cellpose-SAM model to use. `cpsam` is the original (April 2025); `cpsam_v2` is newer (June 2026). Both are zero-shot, so you can switch freely to compare. |
| **Diameter (0 = auto)** | 0 (auto) | ≥ 0 | Expected object diameter, used to rescale the image before segmentation. Leave at 0 — Cellpose-SAM was trained on diameters **7.5–120 px (mean 30 px)** and is robust to size, so only set it manually if your objects are far outside that range. |
| **Flow threshold** | 0.4 | 0.0–3.0 | Maximum allowed flow error for a predicted mask to be kept. **Raise** to detect more (possibly poorly-shaped) objects; **lower** to keep only cleanly-shaped ones. **Ignored (and greyed out) when 3D segmentation is on.** |
| **Cell prob threshold** | 0.0 | −6 to 6 | Pixels with predicted probability above this value are used to build masks. **Lower** to detect more / larger objects; **raise** to drop dim / uncertain detections. |
| **3D segmentation (do_3D)** | ON | checkbox | **On:** cellpose averages flows over the YX / ZY / XZ views and runs dynamics in true 3D — most accurate, but roughly **10–15× slower** (measured on the reference dataset). **Off:** each Z slice is segmented in 2D and stitched across Z by IoU overlap (see *Stitch threshold*). |
| **Stitch threshold** | 0.0 | 0.0–1.0 | IoU overlap required between a mask on one Z slice and the next for them to be joined into one 3D object. **Only used (and only editable) when 3D segmentation is off.** |
| **Batch size** | 8 | 1–128 | Number of 256×256 patches processed together per GPU pass. **Reduce this first if you hit a CUDA out-of-memory error** — it is the main lever for GPU memory use, and has no effect on segmentation quality. |
| **Remove flat (single-slice) segments** | ON | checkbox | Drops any output object that occupies only a single Z, Y or X slice. This is **BEHAV3D's own post-processing**, not a native Cellpose parameter: cellpose leaves flat fragments on the first/last slice (~100–350 voxels vs. ~2400–4000 for real objects on the reference dataset), most common in 2D+stitch mode. Turn off only if your objects genuinely occupy a single slice. |

```{tip}
**Do you want 3D or 2D+stitch?** Keep **do_3D on** for the most accurate result. Switch it **off** for speed, or when the ZY/XZ views are unusable (very anisotropic data) — in that case tune the **Stitch threshold** instead. The two are deliberately kept together as an on/off workflow choice rather than a fine-tuning knob.
```

### Advanced (normalization, filters)

Collapsed by default. Every setting inside defaults to Cellpose's own recommended value or "off", so you can **leave it untouched for routine runs**. The split follows Cellpose's own documentation: routine tuning knobs (diameter, flow/cellprob, do_3D/stitch, batch size) stay in the main panel above; edge-case overrides live here.

| Parameter | Default | What it does |
|---|---|---|
| **Norm percentiles (low / high)** | 1 / 99 | Intensity percentiles rescaled to 0/1 before segmentation. Defaults are fine for normal data. |
| **Normalize across whole Z stack (norm3D)** | ON | On: normalization is computed once for the whole Z stack. Off: each Z slice is normalized independently — turn off if brightness drifts a lot slice-to-slice within one stack. |
| **niter dynamics** | 0 (auto) | Iterations used to run the flow dynamics. 0 = automatic, proportional to the diameter. Long, thin objects may need a higher value (e.g. ~2000). |
| **flow3D smooth** | 0 (off) | Standard deviation of a Gaussian filter applied to the 3D flows before computing dynamics. Increase if you see ring artifacts. Shown in the selected unit. |
| **Sharpen radius** | 0 (off) | High-pass sharpening filter applied before segmentation. Use ~1/4–1/8 of the expected object diameter for dim, blurry images. |
| **Smooth radius** | 0 (off) | Low-pass smoothing filter applied before segmentation. Use for noisy images. |
| **Tile norm blocksize** | 0 (off) | Computes normalization in tiles of this size across the image rather than once per frame, brightening dark areas in unevenly-illuminated images. |
| **Tile norm smooth3D** | 0 | Smooths the tiled normalization across Z. Only matters once *Tile norm blocksize* > 0. |
| **Max size fraction** | 1.0 | Cellpose deletes any mask larger than this fraction of the image area. Cellpose's own upstream default is 0.4, which silently removes large organoids — BEHAV3D defaults it to 1.0 (keep everything) and relies on the size filter in Section 5 instead. |
| **CPU threads (0 = all)** | 0 | Caps CPU threads used for mask post-processing and compression. |
| **Cooldown between frames (s)** | 0 | Idle time inserted after each timepoint. |

```{warning}
**Power-delivery crashes on laptops.** A full Cellpose-SAM run saturates the GPU *and* every CPU core at once, which on some laptops draws more power than the machine can sustain, causing an abrupt shutdown mid-run. If that happens, **cap CPU threads** (a value around one quarter of your core count is a reasonable start) and **add a few seconds of cooldown**. A crashed run always **resumes automatically** from where it left off — progress is tracked in a per-sample sidecar journal — so these are recovery knobs, not an auto-apply toggle.
```

## Section 5 · Size filter (preview before batch)

This is a **post-processing** step, exactly like APOC's [Minimum size](apoc) filter: it **excludes** objects, it does not merge them — a small object touching a larger one is dropped, not absorbed.

| Control | Default | Effect |
|---|---|---|
| **Sample / t** | first sample / 0 | The single sample and timepoint that **Preview segmentation** runs on. Pick one representative of your data so the size histogram reflects real object/artifact sizes. |
| **Min size** | 0 (no bound) | Minimum object volume to keep; smaller objects are removed. Set to remove debris and flat fragments. |
| **Max size** | 0 (no bound) | Maximum object volume to keep; larger objects are removed. |
| **👁 Preview segmentation** | — | Runs Cellpose-SAM once, **unfiltered**, on the chosen sample/timepoint and caches the result as a napari Labels layer. The Min/Max size fields then **re-filter that cached result instantly** — no Cellpose-SAM re-run while you tune the thresholds. |

## Run

| Control | Effect |
|---|---|
| **▶ Run Cellpose-SAM** | Runs Cellpose-SAM for every sample in the loaded metadata, using all the settings in Sections 1–5 (saved automatically on Run). Runs in the background so the GUI stays usable; progress appears in the log and the shared progress bar. |
| **+🛒** | Queues a Cellpose-SAM segmentation step in the [Processing Queue](../../plugin_essentials/processing_queue) to run later, sequentially with the rest of the pipeline. |

## Outputs

| File | Path |
|---|---|
| Per-cell-type instance labels | `<output_dir>/images/<sample>/<sample>_<cell_type>_segments.zarr` |

Cellpose-SAM writes to the **same canonical segments path** as every other method, so tracking and everything downstream consume it identically. There is no model file to store — the `cpsam` / `cpsam_v2` weights ship with the sidecar environment.

## Tuning the result

Because Cellpose-SAM is zero-shot, tuning is mostly picking **do_3D vs 2D+stitch** and the two thresholds, then filtering by size:

- **Touching cells fused into one label** → Cellpose under-segments by design (no watershed splitting). Try raising **Cell prob threshold** slightly, or lower the **Diameter** if auto-diameter overestimated it. If it persists, this is a genuine limitation — split manually after tracking via [Manual Editing](../tracking/manual_editing), or use a watershed-based method ([APOC](apoc) / [ConvPaint](convpaint)) for that cell type.
- **Many spurious small / flat fragments** → keep **Remove flat segments** on and set a **Min size** in Section 5 (preview first to read real object sizes).
- **Large organoids silently disappearing** → check **Max size fraction** is 1.0 (BEHAV3D's default), not Cellpose's upstream 0.4.
- **Out-of-memory (CUDA)** → lower **Batch size** first.
- **Objects over-split in Z** → confirm your metadata pixel sizes are correct (they set the anisotropy), and consider turning **norm3D** on.

## Things this page does **not** claim

- Specific per-timepoint runtimes — they depend strongly on your GPU, volume size, channel count, and do_3D vs 2D+stitch. Profile your own data.
- That `cpsam_v2` always beats `cpsam` — both are zero-shot; compare them on your data.

## See also

- [Segmentation overview](./index.md) — how Cellpose-SAM fits with the other methods.
- [Cellpose (classic)](cellpose) — the trainable/retrainable v3 backbone, and the **dead-mask (Otsu)** step.
- [Cellpose project on GitHub](https://github.com/MouseLand/cellpose) — upstream documentation for the Cellpose-SAM model and every inference parameter.
- [Processing Queue](../../plugin_essentials/processing_queue) — batching segmentation steps.
