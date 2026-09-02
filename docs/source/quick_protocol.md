# ⚡ Quick protocol

The whole BEHAV3D EXPLORER pipeline on one page — the minimal path and the few real choices at each tab. Work the dock tabs **left to right**. Each section links to its full page when you want the detail.

```{note}
Method names and defaults here are **what the pipeline does**, not recommended values. Derive every threshold from your own images (spacing, frame interval, object size) and **preview before batching**.
```

## 📋 Data Preparation

*Full detail: [Data Preparation](data_preparation).*

1. **Set the Output Directory** — nothing else runs until this is set.
2. **Build metadata:** number of samples (one per field of view) → population counts (+ dead channel) → **Configure Cell Types** → **Create Sample Forms** → fill sample 1 → **Fill All from Sample 1** → **Save CSV**.
3. **Load the CSV** — must pass validation to unlock the other tabs.
4. **Check the dimension order** (e.g. `TCZYX`); fix any wrong row inline.
5. **Convert to Zarr** *(optional: trim dead start/end frames)*.

**Population vs Line/Condition:**

| If… | Declare it as… |
|---|---|
| You can tell it apart **in the image** (own or overlapping channel) | a **population** — own segmentation model + track IDs |
| Same biology, deliberately **split across colours** | **one** population, tick **Multicolor** |
| **Identical on screen**, differs only by treatment / donor / subclone | **one** population; use **Line** / **Condition** |

**Required per sample:** name · well · image path · XY & Z size + unit · time interval + unit *(dead-channel index only if ticked)*.

## 🦠 Segmentation

*Full detail: [Segmentation](processing/segmentation/index). Tune **each cell type on its own tab**.*

**Pick one method per cell type:**

| Your situation | Method |
|---|---|
| Already segmented elsewhere (Imaris, Ilastik, script) | **Import** |
| You have a trained Cellpose model for this data | **Cellpose (v3)** |
| No model, but GPU + clean channels + can wait + roughly convex cells | **Cellpose-SAM** (zero-shot) |
| Otherwise: GPU, good signal, simple blob-like objects | **APOC** |
| Otherwise: harder features / touching objects | **ConvPaint** |
| No GPU available | **Pixel Classifier** |

**Pixel-classifier flow (APOC / ConvPaint / Pixel Classifier):** Generate Training Data → paint a few foreground/background pixels → Train → **Run instance segmentation** (preview) → tune instance post-processing → Batch.

**Cellpose / Cellpose-SAM:** set the channel mapping → *(v3)* load the model → Run.

**Key knobs** — instance post-processing ([full detail](processing/segmentation/index.md#instance-post-processing-parameters)):

| Knob | What it does |
|---|---|
| **EDT / Seed threshold** | splits touching objects — **raise** to split more, **lower** to keep them merged |
| **Mask threshold** | foreground cutoff (probability strategies) — **lower** to keep dim cells |
| **Min size** | drops segments below this many voxels (noise filter) |
| **Opening** | smooths boundaries / removes speckles |

## 📍 Tracking

*Full detail: [Tracking](processing/tracking/index). One sub-tab per cell type.*

**Pick the method from overlap, not cell type** — between two consecutive frames, does the object still overlap where it was?

| Between frames… | Method |
|---|---|
| Large overlap, mainly **fragments** | **Fragmentation Tracking** |
| Overlap but masks **touch / join** | **Bounded Propagation** |
| **Near-static**, detection flickers on/off | **Reporter Propagation** |
| **Little / no overlap** (moves ~1 cell diameter or more) | **btrack** |
| Already tracked elsewhere | **Import** |

*(LAP and TrackPy are centroid-distance linkers, also available.)*

1. Per sub-tab: pick method → set params → **Run <cell type> Tracking**.
2. **Inspect in Visualization** — a correctly tracked cell keeps the **same ID / colour across frames**.
3. **Several organoid types in one movie?** *Track all organoids together* **ON** = one shared ID space; **OFF** = each type tracked independently.

**Key knobs** — the propagation methods have few or no parameters; **btrack** is the tunable one:

| Knob | What it does |
|---|---|
| **Max search radius** | furthest an object may move between frames and still link — set from the fastest one-frame displacement |
| **Distance / Time threshold** | (optimizer) largest spatial gap / frame gap to bridge when reconnecting broken tracks |

## 🧪 Feature Extraction

*Full detail: [Feature Extraction](analysis/feature_extraction). One sub-tab per cell type.*

Computes, for every cell at every timepoint, the features every later analysis reads from. Per cell type, choose which feature families to compute — some are forced on:

| Family | Computed for |
|---|---|
| **Intensity**, **Contact** | all cell types (always on) |
| **Movement** | immune / other cell types (always on) |
| **Death** | any cell type, when a dead channel exists (always on) |
| **Morphology** | optional, any cell type |
| **Organoid Invasiveness** | optional, immune only (needs Contact) |

1. Per sub-tab: confirm the feature families → set **Contact Threshold** (µm) → set **Dead mask % threshold** (preview it with **👁 Preview Dead Threshold**).
2. **Run** per cell type, or batch all cell types.
3. **Active Killing** (optional, immune only): an **extra feature-extraction step** run from this tab that adds killing columns to the immune feature table (`is_active_killing`, `killing_efficiency`, `targeted_track_id`, …). Its meaning and parameters are documented with the analyses → [Active Killing](analysis/population_analysis/active_killing).

**Key knobs:**

| Knob | What it does |
|---|---|
| **Contact Threshold** (µm) | max distance between two cells to count as "in contact" (`0` = must physically touch) |
| **Dead mask % threshold** | fraction of a cell's voxels overlapping the dead mask before it is flagged **dead** |

## 📊 Analysis

*Full detail: [Analysis](analysis/index).*

1. **Run [Filtering](analysis/filtering)** first (Tab 6) — the analyses read the filtered tracks, not the raw ones.
2. Run whichever analyses you need. The 📊 Analysis tab has two groups.

**Population Analysis** (📊 → 💀 Death Dynamics) — population- and interaction-level readouts:

| Analysis | What it gives you |
|---|---|
| **[Death Dynamics](analysis/population_analysis/death_dynamics)** | how a target signal rises across each population over time |
| **[Interaction Analysis](analysis/population_analysis/interaction_analysis)** | effector–target contact and interaction over the movie |
| **[Invasiveness](analysis/population_analysis/invasiveness)** | how much of an immune cell's surface is embedded in a structure |
| **[Active Killing](analysis/population_analysis/active_killing)** | contact-associated signal rise attributed to individual effectors *(configured in Feature Extraction)* |

**Single Cell** (📊 → 🧬) — per-cell behaviour:

| Analysis | What it gives you |
|---|---|
| **[State Classification](analysis/single_cell/state_classification)** | per-timepoint behavioural states (HMM) |
| **[Track Classification](analysis/single_cell/track_classification)** | whole-trajectory clustering into behaviour types |

Both Single Cell workflows end in a **Backprojection** step that paints the resulting labels back onto the raw images.

**Key knobs:**

- **Filtering:** minimum track length, plus experiment-duration coverage, dead-at-start, and minimum-size cutoffs.
- **Active Killing:** *observation window* (frames after each contact timepoint to look for the signal rise), *threshold* (a multiplier of the target's own signal, or a fixed absolute rise), and *minimum contact duration*.
- **State / Track Classification:** the number of behavioural **states** / trajectory **clusters** to fit.
