# 🧪 Feature Extraction

Tab 5 of the BEHAV3D EXPLORER dock widget. Feature Extraction takes your **tracked segments** plus the **raw image** and computes, for every cell at every timepoint, a panel of biological measurements (movement, intensity, morphology, contact with neighbours, death). The result is the table that every downstream analysis — Filtering, Death Dynamics, Interaction, and behavioural-state classification — reads from.

![Feature Extraction tab](../_static/screenshots/feature_extraction_tab.png)

## What this tab produces

For each cell type, Feature Extraction writes one combined table:

```
<output_dir>/analysis/<cell_type>/track_features/
    BEHAV3D_<cell_type>_combined_track_features.csv
```

Each row in that CSV = one cell (`TrackID`) at one timepoint (`position_t`). The columns are grouped into six families (see [below](#the-six-feature-families)). Per-sample intermediate CSVs are also written under `trackdata/<sample>/<cell_type>/`.

Optionally, an **Active Killing** detector can run on top of the immune-cell features to flag the timepoints at which an immune cell appears to kill a target. That produces several extra files in `analysis/<immune_type>/active_killing/`.

```{important}
The main **▶ Run Feature Extraction** step writes CSVs only — no plot PDFs. Other tabs produce most QC and behavioural figures ([Filtering](filtering.md), [State Classification](single_cell/state_classification.md), etc.).

Graphics from **this** tab are limited to:

- **Preview Dead Threshold in Viewer** — live napari overlay only (not saved).
- **Active Killing** (immune panel, after baseline extraction) — when you run the analysis or **Display Existing Results**:
  - Per sample: `plots/<sample>/killing_kinetics_summary_<sample>.png` (four panels: efficiency, kinetics, cumulative events, events per cell).
  - All samples: `plots/combined_killing_efficiency_distribution.png`.
  - Top killers: `gallery/<sample>/killing_event_*.gif` (Event GIFs of cropped 3D views (maximum-intensity projections) around each top killing event; count set by **Top-N killers to display**).

See [Active Killing](population_analysis/active_killing.md) for file names and how to read the plots.
```

## Per-cell-type sub-tabs

Like the Tracking tab, Feature Extraction is organised as **sub-tabs per cell type** along the left side. Each cell type declared in the metadata gets its own sub-tab, colour-coded by category:

| Icon | Category |
|---|---|
| 🟣 | Organoid |
| 🔵 | Immune |
| 🟡 | Other |

Multicolor channel splits (e.g. `tcell_1_multicolor`) are **not** shown as separate sub-tabs — the merged track set is used instead.

A collapsible **Active Killing** panel sits at the bottom of the tab, hidden until at least one immune cell type is present.

## The six feature families

Each sub-tab exposes the same six feature-family checkboxes. Some are **mandatory** for certain cell types (the checkbox is forced ON and greyed out, because downstream analyses depend on them):

| Family | What it measures | Mandatory for |
|---|---|---|
| **Movement** | Step length, cumulative distance, mean square displacement, speed, directional persistence | Immune, Other |
| **Intensity** | Mean intensity of each raw channel inside each cell's mask | All cell types |
| **Morphology** | Volume, shape descriptors, surface area, principal axes, orientation | Optional |
| **Contact** | Which neighbours (of any cell type) each cell is touching or close to | All cell types |
| **Organoid Invasiveness** | How much of the immune cell’s surface is in contact with an organoid (0–100%; “invasive” when ≥ about half) | Optional, **immune only**, requires Contact |
| **Death** | Fraction of the cell's mask that overlaps the dead mask, plus a sticky alive/dead flag | Whenever a dead channel exists |

### Per-category summary

| Cell type | Always on (greyed) | Optional |
|---|---|---|
| Organoid | Intensity, Contact (+ Death if dead channel) | Movement, Morphology |
| Immune | Intensity, Contact, Movement (+ Death if dead channel) | Morphology, Organoid Invasiveness |
| Other | Intensity, Contact, Movement (+ Death if dead channel) | Morphology |

## Parameters

Each sub-tab exposes the following spinboxes:

| Control | Default | Range | Units | Meaning |
|---|---|---|---|---|
| **Contact Threshold** | 0.0 | 0.0 – 500.0 | µm | Maximum distance between two cells to be counted as "in contact" (used for distance-based contact columns and `touching_*`). Pixel contact is fixed at one-voxel-diagonal (≈ 1.73 px). Visible only when Contact is checked. |
| **Dead mask % threshold** | 10 % (stored as 0.1) | 0 – 100 % | percent shown, fraction stored | Minimum share of the cell's voxels that must overlap the dead mask before the cell is flagged as dead. The spinbox is shown as a **percentage** (default 10 %), but the value is saved to `behav3d_parameters.yml` and compared as a **0–1 fraction** — 10 % is stored and compared as `0.1`, directly against `percentage_dead_mask` (itself a 0–1 fraction). Set to **0** to skip dead classification entirely. |
| **Preview sample** | first sample | dropdown | — | Sample used by the dead-threshold preview button. |
| **Workers** | balanced default | 1 – (CPU cores − 1) | cores | Parallel workers for image-based computations (morphology / intensity / contact / dead mask). |

```{important}
The **Dead mask % threshold** is **shared across all organoid sub-tabs**. If you change it on one organoid sub-tab, every other organoid sub-tab updates to the same value in real time. Immune and Other cell types each have their own independent threshold.
```

## Preview Dead Threshold in Viewer

The **👁 Preview Dead Threshold in Viewer** button (only visible when Death is enabled — death channel exists ) is the quickest way to set a sensible Dead threshold for a dataset. It:

1. Clears the napari viewer.
2. Loads, for the selected preview sample, the raw channels, the dead mask, and the tracked segments of every cell type.
3. Overlays each cell of the current sub-tab's cell type with a colour based on its dead-fraction:

| Colour | Meaning |
|---|---|
| 🟢 Green | **Alive** — dead-mask overlap **below** threshold |
| 🔴 Red | **Dead** — dead-mask overlap **at or above** threshold |

The overlay updates **live** as you drag the Dead mask % threshold spinner. Hover over any cell to see its track ID and exact dead percentage in the napari status bar.

For organoid sub-tabs, all organoid types are merged into a single overlay so you can see them together.

```{tip}
Picking the dead threshold is dataset-specific, but use these starting points (the value is the fraction of the mask that must overlap the dead signal):

- **Organoids ≈ 0.05** (about 5 % of the mask). It sounds low, but organoids are large 3-D objects, so a few dead cells are only a small volume fraction.
- **T cells ≈ 0.1** (about 10 %, the default). For a clean dead-dye signal there is a clear central red dot.

Adjust from there until the green/red split visually matches the dead cells in your raw image at a timepoint you trust. The binary alive/dead flag is stored alongside the actual intensity and percentage, so you can re-tune later with **▶ Re-run death**.
```

If you already ran extraction with one threshold and want another, change Dead mask % threshold and run ▶ Re-run death

## Active Killing detection

A collapsible **▶ Extended Analysis — Active Killing (Immune Cells)** section sits at the bottom of this tab. It detects, for each effector cell, the timepoints at which the target it touched shows a signal rise large enough to count as killing.

**It is configured and run here, but explained with the population analyses.** It lives in this tab because it needs the per-timepoint contact and signal columns while they are being computed, and it writes its results back into the effector's own feature table. Conceptually it belongs with Death Dynamics, Interaction and Invasiveness, because it is about targets, effectors and contact.

The panel only appears when the metadata contains at least one **immune** cell type, and it requires that cell type's combined feature CSV to exist already — so run baseline Feature Extraction for it first.

**Full explanation, parameters, calibration and outputs: [Active Killing](population_analysis/active_killing.md).**

## Apply, Run, Queue

Each sub-tab has:

- **Apply to all <Category>s** — copy this sub-tab's settings (feature checkboxes, contact threshold, workers, and the non-organoid dead threshold) to every other sub-tab of the same category.
- **Apply to all** — same but across all cell types.
- **▶ Run <CellType> Feature Extraction** — runs feature extraction for this cell type across every sample. If the combined CSV already exists you'll be asked to **Overwrite** or **Skip**.

The main tab (above the sub-tabs) has:

- **▶ Run Batch Feature Extraction (All Cell Types)** — runs every cell type sequentially.
- **+🛒** — queues a Feature Extraction step in the [Processing Queue](../plugin_essentials/processing_queue). The queue snapshots the current GUI state and applies it at run time.

Progress is reported in the **Log** at the bottom of the tab and in the console (progress bars per per-sample loop).

## Column reference

The combined CSV always contains a base set of identity columns, plus whatever the enabled feature families add.

### Identity / time columns (always present)

| Column | Meaning |
|---|---|
| `sample_name` | Sample identifier |
| `TrackID` | Stable track identifier across timepoints |
| `SegmentID` | Per-timepoint segment label (NaN on interpolated rows) |
| `position_t` | Timepoint index |
| `position_x` / `_y` / `_z` | Centroid in µm |
| `relative_time` | Time within this track, starting at 1 |
| `time` | Real time in hours |
| `distance_unit`, `time_unit` | Always `"um"` / `"h"` |
| `interpolated` | `True` if this row was filled in to bridge a missing frame |

Any extra metadata columns you have in the tracks CSV (e.g. `well`, `condition`) are carried over too.

### Movement (when checked)

All measurements are in µm and hours. Movement is computed only after missing timepoints are interpolated, so all columns are continuous along each track.

| Column | Meaning |
|---|---|
| `displacement` | Step length from the previous timepoint (µm) |
| `cumulative_displacement` | Running total of step lengths up to this timepoint |
| `displacement_from_origin` | Straight-line distance from the cell's first position |
| `mean_square_displacement` | Mean squared distance between the current position and every earlier position on the track (a diffusion measure) |
| `directional_persistence` | Cosine of the angle between consecutive step vectors (−1 = U-turn, 0 = right-angle turn, +1 = straight ahead) |
| `speed` | Step length divided by the time between frames (µm/h) |
| `mean_speed` | Speed averaged over the last 10 frames (rolling window); also used to decide `active_*_contact` |

#### How the movement features are computed

Let $\mathbf{p}(t) = (x, y, z)$ be a cell's centroid at timepoint $t$ (in µm), $\Delta t$ the frame interval (in hours), and $\mathbf{s}(t) = \mathbf{p}(t) - \mathbf{p}(t-1)$ the step vector between consecutive frames.

$$
\text{displacement}(t) = \lVert \mathbf{s}(t) \rVert
\qquad
\text{speed}(t) = \frac{\lVert \mathbf{s}(t) \rVert}{\Delta t}
$$

$$
\begin{aligned}
\text{cumulative\_displacement}(t) &= \sum_{i \le t} \lVert \mathbf{s}(i) \rVert \\[4pt]
\text{displacement\_from\_origin}(t) &= \lVert \mathbf{p}(t) - \mathbf{p}(0) \rVert
\end{aligned}
$$

The **mean square displacement** at timepoint $t$ averages the squared distance from the current position to *every* position up to and including $t$:

$$
\text{MSD}(t) = \frac{1}{t+1} \sum_{j=0}^{t} \lVert \mathbf{p}(j) - \mathbf{p}(t) \rVert^{2}
$$

MSD grows roughly linearly in time for random (diffusive) motion and quadratically for directed motion, so its shape distinguishes wandering cells from cells under a directional cue.

**Directional persistence** is the cosine similarity of successive step vectors (clipped to $[-1, 1]$, and 0 wherever a step has zero length):

$$
\text{persistence}(t) = \frac{\mathbf{s}(t-1) \cdot \mathbf{s}(t)}{\lVert \mathbf{s}(t-1) \rVert \, \lVert \mathbf{s}(t) \rVert}
$$

A value near +1 means the cell keeps heading the same way (persistent, directed migration); near 0 means each step's direction is independent of the last (random walk); near −1 means it reverses.

```{note}
The Feature Extraction tab does not produce the broader rolling-window motility statistics used by the state classifier (e.g. straightness, net displacement over a window, median turning angle, run-length summaries). Those are computed later from this table by [State Classification](single_cell/state_classification) when you run the HMM — they are not extra columns in the Feature Extraction CSV.
```

### Intensity (when checked)

One column per raw channel:

| Column pattern | Meaning |
|---|---|
| `mean_intensity_ch<N>` | Mean intensity of channel *N* inside the cell's mask |
| `mean_dead_dye` | Convenience alias for the dead channel's mean intensity (only when the metadata declares a `dead_channel`) |

Only **mean** intensity is computed per channel.

### Morphology (when checked)

When Morphology is on, the following are added. Units: volumes in µm³, lengths in µm, surface area in µm².

| Column | Meaning |
|---|---|
| `nr_pixels` | Voxel count of the mask |
| `volume` | Physical volume = voxel count × voxel volume (µm³) |
| `bbox_volume` | Volume of the axis-aligned bounding box |
| `equivalent_diameter` | Diameter of a sphere with the same volume |
| `axis1_length` / `axis2_length` / `axis3_length` | Longest / middle / shortest principal-axis length (from PCA of the mask's voxel coordinates, spacing-aware) |
| `major_axis_length`, `minor_axis_length` | The longest (`axis1`) and shortest (`axis3`) principal axis |
| `elongation` | `axis1_length / axis2_length` — longest over middle axis |
| `extent` | Volume divided by the **oriented** bounding-box volume (`axis1 × axis2 × axis3`) |
| `oblateness` | Flatness: `1 − axis3/axis1` (0 = not flattened, →1 = disc-like) |
| `prolateness` | Elongation: `1 − axis2/axis1` (0 = not elongated, →1 = cigar-like) |
| `surface_area` | Physical surface area from exposed voxel faces (µm²) |
| `sphericity` | How sphere-like the object is (1 = perfect sphere); from volume and surface area |
| `surface_to_volume_ratio` | `surface_area / volume` |
| `convex_volume` | Volume of the convex hull (`volume / solidity`) |
| `solidity` | `volume / convex_volume` — how much the object fills its convex hull (1 = convex) |
| `orientation_vector` | Unit vector of the longest principal axis (the cell's 3-D orientation) |
| `border_touching_segment` | `True` when the mask touches the edge of the image volume (a QC flag — border objects are partially cropped, so their morphology is unreliable) |

When Morphology is **unchecked**, only the cheap `nr_pixels` and `volume` columns are kept (the others are skipped to save time).

#### How the morphology features are computed

Let $V$ be the object's physical volume and $A$ its surface area. **Surface area** is measured directly from the voxelised mask by counting the exposed faces (each internal face between a foreground and a background voxel), weighted by the physical area of that face — not a smoothed mesh. **Sphericity** compares the object's surface area to that of a sphere of equal volume:

$$
\begin{aligned}
\Psi &= \frac{\pi^{1/3}\,(6V)^{2/3}}{A} \\[4pt]
\text{surface\_to\_volume\_ratio} &= \frac{A}{V} \\[4pt]
\text{equivalent\_diameter} &= \left(\frac{6V}{\pi}\right)^{1/3}
\end{aligned}
$$

$\Psi = 1$ for a perfect sphere and decreases as the object becomes more irregular or elongated.

The three **principal axes** ($a \ge b \ge c$ = `axis1` ≥ `axis2` ≥ `axis3`) come from a principal-component analysis of the mask's voxel coordinates in physical units. The shape descriptors are then:

$$
\begin{aligned}
\text{elongation} &= \frac{a}{b} & \qquad \text{oblateness} &= 1 - \frac{c}{a} \\[4pt]
\text{prolateness} &= 1 - \frac{b}{a} & \qquad \text{extent} &= \frac{V}{a\,b\,c}
\end{aligned}
$$

**Solidity** ($V / V_{\text{convex}}$) measures how much of the object's convex hull it actually fills — low solidity flags protrusions or concavities.

### Contact (always on — mandatory)

For each other cell type `{type}` present in the metadata, two contact columns and a list:

| Column pattern | Meaning |
|---|---|
| `<type>_contact` | `True` if this cell is **pixel-touching** at least one cell of `<type>` (within ~1.7 px) |
| `<type>_contact_on_distance` | `True` if at least one cell of `<type>` is within the **Contact Threshold (µm)** you set |
| `touching_<type>s` | Comma-separated list of `TrackID`s of `<type>` cells within the distance threshold |
| `any_organoid_contact` / `any_organoid_contact_on_distance` | Aggregated over every organoid type |
| `any_immune_cell_contact` / `any_immune_cell_contact_on_distance` | Aggregated over every immune type |

```{tip}
The `*_contact` columns answer **"is it touching?"** (pixel-adjacency, ~1.73 px diagonal — always computed) and the `*_contact_on_distance` columns answer **"is it close enough to interact?"** (within the **Contact Threshold** in µm). The default Contact Threshold is 0 µm, which makes the two nearly equivalent — masks essentially have to touch. Increase it to 1–5 µm to capture near-contact / proximity relationships between cells that don't physically touch.

**Contact Threshold is not a cheap parameter to sweep interactively:** it feeds the downstream interaction analysis, and changing it requires **re-running Feature Extraction** (it changes the contact columns themselves). Plan any tuning as a batch/offline step rather than trial-and-error in a live session.
```

### Organoid Invasiveness (immune only, optional)

When checked (and Contact is also on), for each organoid type `{org}`:

| Column | Meaning |
|---|---|
| `<org>_invasiveness_perc` | Percentage of the immune cell's *surface voxels* (within ~2 µm of its boundary) that are within Contact Threshold of an `<org>` cell |
| `<org>_invasiveness` | `True` when `<org>_invasiveness_perc ≥ 50%` |
| `any_org_invasiveness_perc` | Maximum across all organoid types |
| `any_org_invasiveness` | `True` when any organoid type's invasiveness is ≥ 50% |

This is a **surface-contact fraction**, not a strict "immune voxels inside organoid volume" check — it measures how much of the immune cell's boundary is wrapped by organoid tissue:

$$
\text{invasiveness\_perc} = 100 \times \frac{\#\{\text{immune surface voxels contacting the organoid}\}}{\#\{\text{immune surface voxels}\}}
$$

$$
\text{invasive} = \big(\text{invasiveness\_perc} \ge 50\%\big)
$$

where "surface voxels" are the immune cell's boundary shell (within ~2 µm of its edge).

### Death (when a dead channel exists and the family is checked)

| Column | Meaning |
|---|---|
| `percentage_dead_mask` | Fraction (0–1) of the cell's voxels that overlap the dead mask |
| `nr_dead_mask_pixels` | Absolute voxel count of overlap = `nr_pixels × percentage_dead_mask` |
| `increase_dead_mask` | Change in `nr_dead_mask_pixels` compared to the same cell 10 timepoints earlier — a simple death-rate proxy |
| `smoothed_percentage_dead_mask`, `smoothed_nr_dead_mask_pixels`, `smoothed_increase_dead_mask` | Rolling-mean-smoothed versions of the death signals (≈10-minute window), used by Death Dynamics to damp single-frame segmentation noise |
| `dead` | Boolean flag. Once `percentage_dead_mask` crosses the threshold for the first time, the cell is marked dead from that timepoint onward — the flag is **sticky** and never flips back to alive |

`percentage_dead_mask` is the fraction of the cell's own voxels that overlap the dead mask, and the pixel count scales it by the cell size:

$$
\text{percentage\_dead\_mask} = \frac{\#\{\text{cell voxels overlapping the dead mask}\}}{\#\{\text{cell voxels}\}}
$$

$$
\text{nr\_dead\_mask\_pixels} = \text{nr\_pixels} \times \text{percentage\_dead\_mask}
$$

The three death signals capture death differently: `percentage_dead_mask` is **size-independent** (good when objects differ in size) but saturates; `nr_dead_mask_pixels` is an **absolute** count (cannot go negative, pairs well with an absolute threshold — see [Active Killing](population_analysis/active_killing.md)); and `mean_dead_dye` (from the Intensity family) is the **mean dead-channel intensity** across the whole mask, which is the right choice when you have no dead-cell segmentation, or when the dye is diffuse and fills the whole cell rather than forming a discrete region.

If you set the threshold to 0, the `dead` column is not added.

```{note}
**The dead channel does not have to carry a death dye.** These columns measure how much of an object is occupied by, or how brightly it carries, whatever signal you declared as `dead_channel`. Any reporter that **switches on and stays on** — a differentiation marker, an activation or stress reporter — works here, and the downstream [Death Dynamics](population_analysis/index) and Active Killing analyses work unchanged. Read "dead" as "past the threshold you set".

A **fluctuating** reporter that goes on and off, such as calcium, is a different case: a sticky flag is meaningless for it. Keep its channel out of `dead_channel`, and use its intensity as a behavioural feature in [Single Cell](single_cell/index) analysis instead.
```

## Tips & best practices

- **Always use the Dead Threshold preview before running.** The green/red overlay is the fastest way to catch a wrong threshold — much faster than rerunning extraction and inspecting a CSV.
- **Use Apply to all category for multi-sample experiments.** Once one organoid sub-tab is correctly tuned, propagate the settings to the others; otherwise it is easy to forget a sub-tab and end up with inconsistent features.
- **Workers parallelise image-based work, not movement.** Movement and death-flag calculation are CPU-cheap and serial. The bulk of wall time is in morphology + intensity + contact + dead-mask measurement, which is what the worker count actually speeds up.
- **Reruns are incremental.** Per-sample intermediate CSVs (intensity, contact, morphology, dead mask) are reused if they already exist on disk and you choose **Skip** in the overwrite dialog — useful when you only need to recompute one feature family. Choose **Overwrite** when you change a threshold and need a clean recompute.
- **Active Killing is opt-in for a reason.** Baseline Feature Extraction is fast; Active Killing has to scan every contact event for every immune track per sample and is much slower. Run it only when you actually need killing-event analysis.
- **Active Killing reads the *filtered* CSV when available.** If you have already run [Filtering](filtering.md), Active Killing will use the filtered feature table rather than the raw one — so it sees the same set of tracks you'll be analysing downstream.

## See also

- [Filtering](filtering.md) — the next step. Drops low-quality tracks and produces QC plots.
- [Population Analysis](population_analysis/index) — uses the signal and contact columns for population dynamics, interaction, invasiveness and active killing.
- [Single Cell](single_cell/index) — classifies per-cell behavioural states from the movement / contact / morphology table.
- [Output Directory & File Layout](../plugin_essentials/output_layout) — where the CSVs live.
