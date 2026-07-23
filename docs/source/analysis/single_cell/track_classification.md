# 🛤️ Track Classification

The second inner sub-tab of **Analysis → 🧬 Single Cell**. Where [State Classification](state_classification.md) labels behaviour **timepoint by timepoint**, Track Classification groups **whole trajectories** into clusters based on how their behaviour unfolds over time — so two cells that follow a similar behavioural "story" end up in the same trajectory cluster, even if they are never in the exact same state at the exact same frame.

It runs on the **cell type chosen in the dropdown** at the top of the Single Cell sub-tab (immune / other only). The default method compares the **sequence of behavioural states** along each track, so you normally run [State Classification](state_classification.md) for that cell type first.

![Track Classification sub-tab](../../_static/screenshots/track_classification_tab.png)

```{note}
*Screenshot placeholder.*
```

## How it works

Trajectories are compared with **Dynamic Time Warping (DTW)**, which aligns two sequences even when they differ slightly in length or timing, then groups the most similar ones with hierarchical clustering. Two clustering engines are available:

- **Categorical DTW (default)** — each track becomes its **sequence of behavioural-state labels** (the per-timepoint states from State Classification). DTW compares state *ordering*, then **hierarchical clustering** groups similar sequences. Because it compares *which states occur and in what order* rather than raw speed, it will not split, say, fast vs. slow scanners that follow the same behavioural story.
- **Original feature-based BEHAV3D DTW** (legacy) — each track becomes a **multi-dimensional numeric trajectory** (displacement, speed, dead-dye, contacts). DTW runs on those scaled features, the distances are embedded with **UMAP**, and **K-means** cuts the clusters. This reproduces the classic BEHAV3D trajectory analysis.

## Choosing a clustering method

There are three ways to assign trajectory clusters. Which one you use depends on what you want to compare and how far you are through the pipeline.

```{mermaid}
flowchart TD
    A[Track Classification] --> B{Apply a saved classifier?}
    B -->|Yes| C[Apply pretrained trajectory classifier<br/>load .pkl + states .h5ad]
    B -->|No| D{State Classification done<br/>for this cell type?}
    D -->|No| E[Only Original BEHAV3D DTW available<br/>Step 1 locked into original mode]
    D -->|Yes| F{Compare tracks on…}
    F -->|behavioural state story<br/>recommended| G[Categorical DTW<br/>Run Track Clustering]
    F -->|numeric movement/contact<br/>legacy| H[Original BEHAV3D DTW<br/>tick checkbox in Advanced]
    G --> I[Steps 2–5: rename, classify,<br/>plots, backprojection]
    H --> J[UMAP cluster CSVs + PDFs<br/>Steps 2–5 need categorical .h5ad]
```

| Method | What it does | Requires first | Full Steps 2–5 in the plugin? |
|---|---|---|---|
| **Categorical DTW** (default) | Compares each track's sequence of behavioural-state labels | [State Classification](state_classification.md) for this cell type | ✅ Yes |
| **Original feature-based BEHAV3D DTW** (legacy) | Clusters on raw movement/contact feature trajectories; reproduces the original BEHAV3D paper's track clustering | [Filtering](../filtering.md) (filtered track-features CSV) | ⚠ No — clustering only; rename/classify/plots/backprojection need the categorical `.h5ad` |
| **Apply pretrained classifier** | Assigns a saved set of clusters to a new dataset | A saved classifier `.pkl` + matching states `.h5ad` | Skips Steps 1–3 |

```{note}
**Applying a saved classifier instead of clustering.** The checkbox **Apply pretrained trajectory classifier** at the very top hides Steps 1–3 and shows a panel where you **Browse…** to a saved classifier (`.pkl`) and its matching state data (`.h5ad`; auto-filled if State Classification output exists for the cell type), then click **▶ Apply Pretrained Classifier**. Use this to reproduce the same trajectory clusters on a new dataset.
```

```{tip}
**If State Classification hasn't been run** for the selected cell type, the sub-tab shows a warning and **locks Step 1 into Original BEHAV3D DTW** (the checkbox is ticked and greyed out), because the categorical method needs the behavioural-states file. Steps 2–5 stay disabled until you run State Classification, after which the sub-tab automatically reverts to the categorical method.
```

## Step 1 — Track Clustering

| Control | Default | Meaning |
|---|---|---|
| **Trajectory size** | 100 | Number of timepoints each trajectory is resampled to before comparison. Set it to its minimum to switch to **Variable-length** mode (compare tracks at their native lengths). |
| **N clusters** | 6 | How many trajectory clusters to cut the hierarchy into. |
| **Random seed** | 123 | Fixes the random parts of clustering so re-runs with identical settings are reproducible. |

Click **▶ Run Track Clustering** to run it (in the background, with a progress bar and **Log**). Next to it: **+🛒** queues the step for a batch run, and **👁** opens the result once it exists.

### ⚙ Advanced Configuration (DTW)

Collapsed by default; the defaults are sensible for most data.

| Control | Default | Meaning |
|---|---|---|
| **Window** | blank (unconstrained) | Sakoe–Chiba band width — the maximum time shift DTW may use when aligning two tracks. Smaller values forbid large temporal shifts (faster, stricter). |
| **Max dist** | blank (off) | Early-abandon threshold: stop computing a pair's distance once it exceeds this value (speed-up for large datasets). |
| **Penalty** | blank (off) | Extra cost added for each warping step, discouraging excessive stretching. |
| **Psi** | blank (off) | Psi-relaxation: lets DTW ignore a few timepoints at the start/end of each track (useful when tracks begin or end mid-behaviour). |
| **Linkage** | average | Hierarchical-clustering linkage (`average`, `complete`, `single`). |
| **Trim mode** | last | When a track is longer than *Trajectory size*, whether to trim from the `last` or `first` timepoints. |
| **Missing policy** | Keep as category | How to treat missing timepoints: `Keep as category` treats "missing" as its own state, or `Drop missing timepoints`. |
| **Parallel DTW computation** | on | Use multiple CPU cores for the (expensive) pairwise distance computation. |
| **Save distance matrix CSV** | off | Also write the full pairwise DTW distance matrix to CSV (large for many tracks). |

### Original feature-based BEHAV3D DTW

Ticking **Use original feature-based BEHAV3D DTW clustering** (at the bottom of Advanced) switches Step 1 to the legacy engine. When it is on, the **Advanced DTW** panel is hidden and the **UMAP** controls appear instead, the run button relabels to **▶ Run Original BEHAV3D DTW**, and the checkbox is mirrored at the top of Step 1 so you can switch back.

| Control | Default | Meaning |
|---|---|---|
| **Trajectory size** | 100 | Track length for feature-based clustering. In this mode tracks are cut to **exactly** this length (minimum *and* maximum), so all trajectories are equal length before DTW. |
| **N clusters** | 6 | Number of trajectory clusters (K-means on the UMAP embedding). |
| **UMAP n_neighbors** | 15 | UMAP neighbourhood size for the embedding (larger → more global structure). |
| **UMAP min_dist** | 0.1 | UMAP minimum distance (smaller → tighter, more clumped embedding). |

Run it with **▶ Run Original BEHAV3D DTW**. It reads the **filtered track-features CSV** (from [Filtering](../filtering.md)), applies the original BEHAV3D feature scaling, and writes UMAP cluster tables and diagnostic PDFs (UMAP plot, feature heatmap, cluster-percentage bars) under `analysis/<cell_type>/results/`.

```{important}
Original mode performs **clustering only**. The integrated Steps 2–5 (rename, train/apply classifier, exemplar plots, backprojection) read the categorical trajectory `.h5ad` from `behavorial_trajectories/`, which only the **Categorical DTW** run produces. After an original-only run, inspect the results through the **Results** panel or the CSV/PDF files under `results/`; run the categorical method if you need the full downstream workflow.
```

## Step 2 — Rename Track Clusters

Freshly computed clusters are numbered. **✏ Rename Track Clusters** opens a dialog where you give each trajectory cluster a meaningful name (e.g. *super-engager*, *engager*, *killer*, *seeder*, *static*); the names are written back into the cluster data so downstream plots and **Step 5 — Backprojection** use them. Giving two clusters the same name merges them. The button enables itself once clustering has produced data, and a status line reports how many trajectories were loaded. **👁** reopens the renamed result.

## Step 3 — Classify Tracks

Once you are happy with the clusters, you can train a classifier so the **same** trajectory clusters can be assigned to new data without re-running the full DTW clustering. This step has two groups: **Train RF Classifier on Named Clusters** and **Apply Classifier to New Data**.

| Button | What it does | Enabled when |
|---|---|---|
| **▶ Train RF Classifier** | Fits a random-forest classifier on the current named clusters. | Clustering has been run. |
| **▶ Apply Classifier** | Assigns clusters to the data using the trained classifier (browse to a classifier and its state data, or use the auto-filled paths). | A trained classifier exists. |

Both have **+🛒** (queue) and **👁** (view) buttons.

### ⚙ Advanced Configuration (Track Classifier)

| Control | Default | Meaning |
|---|---|---|
| **n_estimators** | 100 | Number of trees in the random forest. More trees = steadier but slower. |
| **Test holdout** | 0.20 | Fraction of trajectories held out to estimate classifier accuracy (0.05–0.5). |

## Step 4 — Create Plots

| Control | Default | Meaning |
|---|---|---|
| **Exemplars / cluster** | 10 | How many representative trajectories to show per cluster. |
| **Overview statebars** | on | Include the per-cluster state-bar overview pages. |
| **Backprojection PDFs** | off | Also render per-exemplar backprojection figures to PDF. |
| **Backprojection MP4** | off | Also render exemplar backprojection movies. |

- **▶ Create Exemplar PDFs** — produces a PDF of representative tracks per cluster (optionally with state bars and backprojection figures/movies).
- **▶ Create Diagnostics** — produces clustering-quality diagnostic plots.

Each has a **👁** button to reopen its PDF.

## Step 5 — Backprojection

The final step paints the **trajectory clusters back onto the raw images**, so you can confirm that cells assigned to the same cluster really do behave alike. It is built directly into this sub-tab and works on the cell type selected at the top of Single Cell. It behaves exactly like State Backprojection, only it colours by trajectory cluster rather than per-timepoint state.

### Live overlay in napari

The **Live Napari Layer Backprojection (Tracks)** panel overlays coloured cluster labels onto the selected sample's image.

| Control | Default | Meaning |
|---|---|---|
| **Sample** | — All samples — | Which sample to overlay. **— All samples —** uses the first available sample. |
| **Color by** | behavioral_trajectory_cluster | Which trajectory-cluster label to colour by (`behavioral_trajectory_cluster`, `dtaidistance_cluster`, or `track_cluster`, depending on which clustering you ran). |
| **Opacity** | 80 % | Opacity of the coloured overlay (10–100 %). |

Click **▶ Show Track Backprojection in Napari** to load the overlay — this produces a napari layer only and writes nothing to disk.

### Export to PDF / MP4

Open the collapsible **⚙ Export Options** to render the same overlay to file.

| Control | Default | Meaning |
|---|---|---|
| **DPI (PDF)** | 150 | Rendering resolution for the PDF (50–600). |
| **PDF** | on | Produce a PDF. |
| **MP4** | off | Produce a movie. |

Click **▶ Export Track Backprojection** to run the export in the background; the **Log** reports where the files were written.

## Outputs

Track Classification writes its results under:

```text
<output_dir>/analysis/<cell_type>/behavorial_trajectories/
```

You will find there, depending on which steps you ran:

- The **track-cluster data** as an `.h5ad` file (one trajectory per row, with its cluster label).
- A **track-cluster table** (`BEHAV3D_<cell_type>_track_clusters.csv`) once the classifier is applied.
- The trained **classifier** (`classifier_<cell_type>.pkl`).
- **Exemplar** PDFs (`exemplar_tracks*.pdf`) and **diagnostics** PDFs (`diagnostics*.pdf`).
- Optionally the **DTW distance matrix** CSV (if you ticked *Save distance matrix CSV*).

The **Original feature-based BEHAV3D DTW** engine instead writes its UMAP cluster tables (`BEHAV3D_<cell_type>_UMAP_clusters.csv`, `..._combined_track_features_clustered.csv`, `..._UMAP_cluster_percentages.csv`) and diagnostic PDFs under `analysis/<cell_type>/results/`.

```{tip}
The folder name on disk is `behavorial_trajectories`. The easiest way to reopen any of these results is the **👁** buttons or the shared **Results** panel rather than browsing the path by hand.
```

## Tips & best practices

- **Run State Classification first** (for the categorical method). The default clustering compares behavioural-state sequences, so it needs the per-timepoint states to exist for the cell type.
- **Start with fewer clusters.** Begin around the default *N clusters* and increase only if distinct trajectory types are being merged. Clusters can always be **merged afterwards** during renaming (give two clusters the same name), so err on the side of a few too many rather than too few.
- **Leave Linkage on `average`.** `complete` gives comparable results and is worth trying; **`single` rarely works well** for these distances. Agglomerative clustering is preferred over k-means here, with the caveat that the resulting UMAP embedding can look poor even when the clusters themselves are sensible.
- **Use Variable-length mode for uneven tracks.** If your tracks differ a lot in length and resampling distorts them, switch *Trajectory size* to its minimum (Variable-length) so DTW compares native lengths.
- **Save the distance matrix only when you need it.** It grows with the square of the number of tracks and is rarely needed for routine analysis.
- **Queue the heavy steps.** DTW clustering and classifier training are CPU-intensive — use the **+🛒** buttons to run them unattended behind your other pipeline steps.

## See also

- [State Classification](state_classification.md) — the per-timepoint behavioural-state workflow this builds on.
- [Feature Extraction](../feature_extraction.md) — computes the per-track features the clustering consumes.
