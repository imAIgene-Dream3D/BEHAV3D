# 🔬 State Classification

The first inner sub-tab of **Analysis → 🧬 Single Cell**. It assigns every cell, at every timepoint, to a recurring **behavioural state** — a characteristic "mode" of movement, intensity, morphology and contact (for example *slow & in contact* vs *fast & free*). States are discovered with a **Hidden Markov Model (HMM)**, which models behaviour as a sequence that tends to stay in one state for a while before switching, rather than treating each timepoint independently.

The analysis runs on the **cell type chosen in the dropdown at the top of the Single Cell sub-tab** (immune / other only) and reads that cell type's per-track features table from [Feature Extraction](../feature_extraction).

![State Classification sub-tab](../../_static/screenshots/state_classification_tab.png)

```{note}
*Screenshot placeholder.*
```

## Two ways to run

At the very top is a checkbox: **Apply existing behavioral state classification**.

- **Unchecked (default)** — *train* a fresh model on the current cell type. This is the **Step 1 — State Clustering** workflow described below.
- **Checked** — *apply* a previously saved model to this dataset instead of fitting a new one. A small panel appears where you **Browse…** to a saved HMM deployment artifact (a `.pkl` file) and click **▶ Apply saved HMM artifact**. Use this to classify a new experiment with exactly the same state definitions you established earlier, so results are comparable across datasets.

The rest of this page describes training a new model.

## Step 1 — State Clustering

This is where you choose what the model "sees" and how many states to fit. The controls are grouped into collapsible sections.

### Feature Selection

This section is populated **from the actual feature columns** found in the selected cell type's features CSV, so you only ever see features you really computed.

- **Timepoint features** — checkboxes for the per-timepoint measurements from Feature Extraction, grouped by family (movement, intensity, morphology, contact, death, …). Tick the ones that define the behaviour you care about. If no features appear, run [Feature Extraction](../feature_extraction) for this cell type first.
- **Window features** — features computed over a short rolling window rather than a single frame:
  - **Window size** (default `5`, range 1–500) — how many timepoints each window spans. Keep it around 5 for motility behaviour; **set it to 1 when the events you care about happen at a single timepoint** (e.g. calcium-intensity peaks), so a window doesn't smear them out. The HMM still assigns a state *per timepoint* regardless of this setting — that is what a hidden Markov model does; the window only controls how the window-based features are aggregated.
  - **net_displacement**, **straightness**, **mean_square_displacement** — tick the window-based motility summaries you want added. These are computed here, from the per-timepoint feature table — they are **not** columns in the Feature Extraction CSV.

```{tip}
Start small and biological. A handful of well-chosen features (e.g. speed, a contact column, and one window feature such as straightness) usually gives cleaner, more interpretable states than throwing in every available column.
```

### Feature Processing

Optional clean-up applied to the chosen features before fitting:

| Control | Default | Meaning | When to use |
|---|---|---|---|
| **Log scaling** (per-feature checkboxes) | none | Optionally apply `log1p` to the ticked features before the model standardises all continuous features. | Use **selectively** for strongly right-skewed non-negative measurements, not blanket across all features. Use **📊 Preview feature distributions** to see which features are skewed first. |
| **Smooth window** | 5 | Rolling-average smoothing applied to features (1 = no smoothing). Usually matched to the Window size. | Damps single-frame segmentation errors (e.g. a cell that looks elongated for one frame through under-segmentation). Since no sticky HMM is used by default, this smoothing does that job. |
| **Low percentile cap** | 0.00 | Clip values below this quantile (0 = off). | Tame extreme low outliers. |
| **High percentile cap** | 1.00 | Clip values above this quantile (1 = off). | Off by default — in practice percentile clipping tended to hurt the HMM more than it helped. Reach for it only when tracking errors produce large outliers; inspect the distribution first. |

```{tip}
**Log scaling — `speed` is the classic case.** Speed is heavily zero-inflated (most tracked objects barely move most of the time) with a few high-speed excursions. Left untransformed, that skew lets the rare large values dominate any distance-based clustering and drown out the common near-zero values. Log-transforming compresses the long tail and gives the low-speed majority proportionally more influence. Check the histogram (📊 button) before deciding — don't apply it to features that aren't skewed.
```

### Binary Group Selection

Checkboxes for the **categorical / true-false** columns in the features table (for example a `dead` flag, an `*_contact` column, or an active-killing flag). These are deliberately **kept out of the HMM fit** — clustering on binary columns tends to create artificial states — and instead **subdivide** each behavioural state afterwards. A state such as *slow-moving* can be split into *slow-moving & in contact* vs *slow-moving & not in contact*. This is what later produces the **full behavioural clusters** (see Step 2).

### Number of states

**n_states** (default `4`, range 2–50) — the number of behavioural states the HMM fits. The value is used directly.

### ⚙ Advanced Configuration

Collapsed by default. Most users never need to touch these — sensible defaults are pre-filled.

| Control | Default | Meaning |
|---|---|---|
| **Start offset** | 1 | Number of leading timepoints to skip per track before classifying, backfilling the skipped frames from the next timepoint's state. **Leave at 1; changing it is not recommended.** The first timepoint of a track has no speed by definition (there is no previous position), so without this offset every track would start in a "static" state. Window features also use an *expanding* rolling window (timepoints 1–2, then 1–3, …) until the full window is available. |
| **Skipped frames** | backfill | What to do with those skipped frames: `backfill` assigns them the first classified state; `leave_unassigned` leaves them out. |
| **Covariance type** | full | Shape of each state's feature distribution (`full`, `diag`, `spherical`, `tied`). `full` is the most flexible; `diag` is a lighter, more stable choice when you have many features or few cells. |
| **n_iter** | 200 | Maximum training iterations. Increase if the model reports it has not converged. |
| **tol** | 0.001 | Convergence tolerance; training stops when the fit improves by less than this. |
| **min_covar** | 0.001 | Floor on the variance of each state — guards against numerically degenerate states. |
| **Sticky HMM** | off | When on, makes states more **persistent** (less flickering between states from frame to frame). |
| **kappa / alpha** | 8.0 / 1.0 | Strength of the stickiness (shown only when Sticky HMM is on). Higher `kappa` = states stay put longer. |
| **Random seed** | 123 | Fixes the random initialisation so a re-run with identical settings gives identical states. |

```{tip}
If your states flicker rapidly along a track (a cell flips state every frame), enable **Sticky HMM** and/or raise the **Smooth window** rather than reducing the number of states.
```

### How the states are computed

After the chosen per-timepoint and window features are assembled, smoothed and standardised, they are modelled with a **Gaussian Hidden Markov Model**. An HMM treats each track as a sequence that occupies one of `n_states` hidden behavioural states at every timepoint, where each state has its own multivariate-Gaussian feature distribution (shape set by *Covariance type*) and a transition matrix gives the probability of switching between states. Because transitions are penalised, the model naturally produces **temporally coherent** state assignments — a cell tends to stay in a behaviour for a while rather than flickering — which is exactly what makes it suited to behaviour over time rather than clustering each frame independently. The fitted model then assigns every timepoint its most likely state (Viterbi decoding).

The **window features** summarise motion over each rolling window of length $w$. With $\mathbf{p}(t)$ the position and $\mathbf{s}(i)=\mathbf{p}(i)-\mathbf{p}(i-1)$ the steps inside the window:

$$
\begin{aligned}
\text{net\_displacement} &= \lVert \mathbf{p}(t) - \mathbf{p}(t-w+1) \rVert \\[4pt]
\text{straightness} &= \frac{\text{net\_displacement}}{\sum_i \lVert \mathbf{s}(i) \rVert}
\end{aligned}
$$

**Straightness** ranges from 0 (a cell that wanders but ends where it started) to 1 (perfectly straight travel), and the window **mean_square_displacement** is the mean squared distance of each in-window position from the window's start — together these separate directed migration from local scanning.

### Running it

**▶ Run State Classification** fits the model in the background (a progress bar and **Log** show what is happening; the rest of the GUI stays responsive). The **👁** button next to it becomes available once a model exists and lets you re-open the classified result. The fitted, classified data is written under the cell type's `behavioral_states` folder in the output directory.

```{note}
You can run State Classification interactively with the button, or add it to the [Processing Queue](../../plugin_essentials/processing_queue) with the **+🛒** button next to it (the queue exposes it as **🔬 State Clustering**, with separate **Train** / **Apply State Classifier** queue steps for the apply-existing workflow).
```

## Step 2 — Rename Clusters

Freshly fitted states are numbered, not named. This step lets you give them meaningful biological labels. The buttons enable themselves once the matching clusters exist in the model.

| Button | Renames | Available when |
|---|---|---|
| **✏ Rename Primary Dynamic State Clusters** | The raw HMM states (the *intrinsic* behavioural states). | A model has been fitted. |
| **✏ Rename Full Behavioral Clusters (Binary Groups)** | The states **after** they were subdivided by the binary groups you selected in Step 1. | Available after the intrinsic states are renamed (the full clusters are created when intrinsic states are combined with the binary groups). |

Each button opens a dialog listing the current clusters with editable names; the new names are saved back into the classified data so all downstream reports and backprojection use them. **Giving two clusters the same name merges them** — a handy way to collapse states that mean the same thing biologically (for example merging several motion states into one *scanner*). A status line above the buttons tells you how many intrinsic and full clusters were found.

```{important}
Rename the **primary dynamic states first**. The full behavioural clusters are built from the intrinsic states, so the "Full Behavioral Clusters" button only becomes meaningful once the intrinsic states exist (and ideally have been named).
```

## Step 3 — Reports

Two report buttons summarise the classification across your samples. Each has a **👁** button to reopen its PDF.

- **▶ State Composition Report** — what fraction of time each state occupies, broken down by sample and (optionally) by metadata groupings. The **Group composition plots by** list below it lets you pick one or more metadata columns to split the composition by (Ctrl/Cmd-click for several).
- **▶ State Transition Report** — how often cells switch between states (the transition structure of the behaviour).

Reports are saved as PDFs and can be reopened at any time with their **👁** button or from the shared **Results** panel.

```{note}
Run a report only after you are happy with the classification (and ideally after renaming), since the report uses whatever state labels are currently stored.
```

## Step 4 — Backprojection

The final step paints the behavioural states **back onto the raw images**, so you can verify frame by frame that the computed labels match what the cells are actually doing. It is built directly into this sub-tab and works on the cell type selected at the top of Single Cell.

### Live overlay in napari

The **Live Napari Layer Backprojection** panel overlays coloured state labels onto the selected sample's image.

| Control | Default | Meaning |
|---|---|---|
| **Sample** | — All samples — | Which sample to overlay. **— All samples —** uses the first available sample. |
| **Color by** | full_behavioral_cluster | Which state label to colour by: `full_behavioral_cluster` (states subdivided by binary groups) or `intrinsic_behavioral_cluster` (the raw HMM states). |
| **Opacity** | 80 % | Opacity of the coloured overlay (10–100 %). |

Click **▶ Show State Backprojection in Napari** to load the overlay — this produces a napari layer only and writes nothing to disk.

### Export to PDF / MP4

Open the collapsible **⚙ Export Options** to render the same overlay to file.

| Control | Default | Meaning |
|---|---|---|
| **DPI (PDF)** | 150 | Rendering resolution for the PDF (50–600). |
| **PDF** | on | Produce a PDF. |
| **MP4** | off | Produce a movie (useful for presentations and time-lapse playback). |

Click **▶ Export State Backprojection** to run the export in the background; the **Log** reports progress and where the files were written (under `analysis/<cell_type>/behavioral_states/backprojection/`).

```{tip}
Use Backprojection as a sanity check: scrub through time and the colours should change when the cells visibly change behaviour. If they don't, revisit the features or the number of states above. Rename your clusters (Step 2) before exporting so legends are publication-ready.
```

## Outputs

**The classified data** is written under the cell type's folder in the output directory:

```text
<output_dir>/analysis/<cell_type>/behavioral_states/
```

It is stored as an `.h5ad` data file that holds, for every cell at every timepoint, its assigned **intrinsic state** and its **full behavioural cluster** — using your chosen names after renaming. This is the file the Reports and the **Step 4 — Backprojection** step read.

The HMM fit diagnostics, including the state-feature heatmap used to interpret each state, are saved at:

```text
<output_dir>/analysis/<cell_type>/behavioral_states/processing/hmm_behavioral_classification/quality_control/raw/behavioral_clustering_diagnostics.pdf
```

State Classification uses an **HMM on per-timepoint and rolling-window features**. The separate **Track Classification** tab contains whole-trajectory clustering, including the legacy DTW workflow; the two analyses answer different questions and should not be described interchangeably.

The two reports are each saved as a PDF (the composition report also writes a CSV of the underlying curves).

**The State Composition Report** answers *"how much time do cells spend in each state?"*. It contains:

- **Relative state composition per sample** (stacked bars showing the proportion of each state).
- **Absolute cell-count composition per sample** (the same, but as raw counts).
- **Pooled per-sample summary bars** of overall state proportions.
- Any **metadata groupings** you chose in the *Group composition plots by* list.

**The State Transition Report** answers *"how do cells move between states?"*. It contains:

- **Transition-matrix heatmaps** as both raw counts and row-normalised probabilities, plus versions that exclude self-transitions (state → same state) so the switching pattern stands out.
- **All-pairs Sankey diagrams** visualising the flow between states.
- **Rankings of the most common state sequences** (n-grams) cells follow over time.

```{note}
Reopen either report at any time with its **👁** button or from the shared **Results** panel rather than by typing a path.
```

## Tips & best practices

- **Choose features deliberately.** The states are only as meaningful as the features you feed in. Begin with a small, interpretable set and expand only if states are blurry.
- **Use a fixed seed for reproducibility.** Leave **Random seed** at its default (or note the value you used) so re-runs and collaborators get the same states.
- **Save and reuse a model across experiments.** Once you have a state definition you trust, apply it to new datasets via **Apply existing behavioral state classification** so states stay comparable.
- **Pair with Backprojection.** After classifying, use **Step 4 — Backprojection** above to paint the states back onto the raw images and sanity-check that they correspond to real behaviour.

## See also

- [Feature Extraction](../feature_extraction) — computes the per-timepoint features the classifier consumes.
- [Filtering](../filtering) — remove short / low-quality tracks before classifying.
- [Track Classification](track_classification.md) — groups whole trajectories by their behaviour over time (DTW), with a trainable classifier.
