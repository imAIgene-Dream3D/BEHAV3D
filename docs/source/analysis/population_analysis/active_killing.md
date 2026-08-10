# 🎯 Active Killing

Where Death Dynamics asks how a population's signal rises overall, **Active Killing attributes that rise to individual effector cells** — it flags the timepoints at which the target an effector touched shows a signal increase large enough to count as a killing event.

```{important}
**Active Killing is configured and run in the [Feature Extraction](../feature_extraction) tab**, not in the Analysis tab, because it needs the per-timepoint contact and signal columns while they are being computed. Look for the collapsible **▶ Extended Analysis — Active Killing (Immune Cells)** section at the bottom of that tab. It appears only when the metadata contains at least one **immune** cell type, and it requires that cell type's combined feature CSV to exist already.
```

```{note}
**The effector must be declared as an `im_` (immune) cell type.** The effector dropdown is built from immune types only, so an `ot_` (other) population cannot be selected as a killer — although it *can* be selected as a target, and it *can* act as an effector in Interaction Analysis.

Nothing in the calculation depends on the population being immune. The prefix is a **processing category, not a biological claim**: if your killer is a neutrophil, a non-immune cytotoxic line, a parasite, or anything else that engages and damages another object, declare it with the `im_` prefix and Active Killing works normally. The same prefix is what makes **invasiveness** features available and **Movement** mandatory.
```

## How it works in plain terms

For each contact event, the detector looks at every target the effector actually touched, independently, anchored at the moment contact started:

> *From the start of contact to N timepoints later, does the signal on this specific target rise enough to count as killing?*

Whichever touched target clears its own threshold by the largest margin is reported as the event's target. There is no sample-wide or cross-target averaging involved — each target is judged only against its own signal.

If yes (and contact lasted long enough), the whole observed contact duration is flagged `is_active_killing = True`.

```{note}
**Caveat:** when several effectors contact the same target at once, all of them may be flagged as active killers for that event, even if only one actually did the killing — the detector attributes a signal rise to every qualifying contact, not to a single culprit.
```

### The killing rule (how it is computed)

Each touched target is judged **only against its own signal**, anchored at the timepoint the contact **started** ($t_c$) — there is no background rate and no averaging across targets or across the sample. Let $W$ be the observation window and $D$ the chosen signal column. For one touched target:

$$
D_{\text{start}} = D(t_c)
\qquad
D_{\text{end}} = D(t_c + W)
\qquad
\Delta D = D_{\text{end}} - D_{\text{start}}
$$

If the track ends before $t_c + W$, $D_{\text{end}}$ is read at the target's last available timepoint. The increase $\Delta D$ is compared to a threshold $\theta$ set by the mode:

$$
\theta =
\begin{cases}
\theta_{\text{abs}} & \text{(absolute mode)} \\[4pt]
D_{\text{start}} \times (m - 1) & \text{(multiplier mode)}
\end{cases}
$$

Multiplier mode is therefore equivalent to requiring $D_{\text{end}} > D_{\text{start}} \times m$ — the signal must grow past $m\times$ its value at contact start. A baseline of exactly $D_{\text{start}} = 0$ is replaced by $0.1$ so the threshold isn't trivially zero. The target is flagged as killed when

$$
\Delta D > \theta
\qquad\text{with}\qquad
\text{killing\_efficiency} = \frac{\Delta D}{\theta}
$$

Every touched target is scored this way, and the one with the **highest** `killing_efficiency` is reported as the event's target (`targeted_track_id`, set only when it is actually flagged as killed). That single per-event verdict is then written onto **every** timepoint of the observed contact, and the event is only counted once the contact has lasted at least the **minimum contact duration**.

## Parameters

| Control | Default | Range | Meaning |
|---|---|---|---|
| **Immune cell type** | first immune type | dropdown | Which effector tracks to analyse. Immune (`im_`) types only — see the note above. |
| **Observation window** | 5 | 1 – 100 timepoints | How many frames after contact starts to measure the signal rise on each touched target. |
| **Death signal column** | `percentage_dead_mask` | dropdown of `percentage_dead_mask`, `mean_dead_dye`, `nr_dead_mask_pixels` | Which target column is read as the signal. If your `dead_channel` carries a reporter other than a death dye, this is that reporter's column. |
| **Killing threshold multiplier** | 1.5 | 0.1 – 20.0 | If absolute mode is off: a target's signal must reach at least `signal_at_contact_start × multiplier` by the end of the window to count as killing. Scales with each target's own starting signal (a signal of exactly 0 is treated as 0.1 to avoid a trivial threshold). |
| **Use absolute threshold instead of multiplier** | OFF | checkbox | When on, the multiplier is replaced by a fixed value (next field). Recommended together with `nr_dead_mask_pixels`, since a flat pixel-count cutoff is easier to reason about than one on a fraction/intensity scale. |
| **Absolute threshold** | 0.0 | 0.0 – 10000.0 | Fixed minimum signal increase (only used when "Use absolute threshold" is on). |
| **Min contact duration** | 1 | 1 – 50 timepoints | Minimum consecutive timepoints an effector must be in contact with the same target before a killing event can be counted. |
| **Top-N killers to display** | 5 | 1 – 50 | Used by the preview button below. |

## Choosing the window and threshold

None of these values are universal — they depend on your **biology** and your **imaging cadence**.

- **Observation window** — how long, after contact, you expect the signal to register. **Set it as a duration, then convert.** Decide how long the event plausibly takes in minutes, divide by your `time_interval`, and enter that many timepoints: a 5-timepoint window is 10 minutes at a 2-minute interval but 2.5 minutes at 30 seconds. Never carry a timepoint count over from another experiment without redoing that conversion. Try a few lengths around your estimate and check which gives sensible results.
- **Minimum contact duration** — is a single-timepoint touch enough for killing to plausibly begin, or must contact persist to count as a real attempt rather than a passing brush? Again, judge against your interval.
- **Absolute threshold vs. multiplier** — the **absolute** threshold (a fixed increase in the death signal, e.g. +20–30 dead pixels) is the general recommendation and pairs naturally with the pixel-count death signals. The **multiplier** (fold-change over each cell's own pre-contact baseline) is for two situations: a single target line where baseline differences don't matter, or wells with a mix of dying and non-dying cells where cell-to-cell baseline variation makes a single fixed threshold unreliable. With **two target lines that have different baseline death rates the multiplier becomes unreliable** — the same fold-change is a very different absolute increase per line — so prefer the absolute threshold there.

```{tip}
**Calibrating an absolute pixel threshold to "one dead cell".** The right absolute value depends directly on cell size and pixel resolution, so recompute it for your own setup rather than reusing someone else's number. As a worked estimate, a cell of diameter $d$ imaged at resolution $r$ (µm/pixel) fills roughly

$$
\text{area} \approx \frac{\pi}{r^{2}} \left(\frac{d}{2}\right)^{2} \ \text{pixels}
$$

e.g. a ~10 µm cell at 1 µm/pixel is about $\pi \times 5^2 \approx 80$ pixels if the whole cell fills with dye. "One dead cell" could therefore mean anywhere from ~20 pixels (partial/rim staining, smaller cells, coarser resolution) to ~80+ pixels (full-cell fill, larger cells, finer resolution). Use this as a starting point, then fine-tune by eye against a timepoint you trust.
```

## Buttons

- **👁 Load Top Killers in Viewer** — adds a Points layer per top-killing track at every timepoint flagged as killing, so you can scrub through the movie and verify visually.
- **▶ Run Active Killing Analysis** — runs the detector for every sample and writes the output CSVs.
- **+🛒** — queues the analysis with the current parameters.

## 🎯 Active Killing outputs

Written under `<output_dir>/analysis/<immune_type>/active_killing/`:

| File | Contents |
|---|---|
| `BEHAV3D_<immune_type>_advanced_track_features.csv` | The full effector feature table with extra columns: `is_active_killing`, `killing_efficiency`, `targeted_track_id`, `contact_event_id`, and a `death_signal_increase_<N>tp` column where N is your observation window. |
| `active_killing_per_timepoint_<immune_type>.csv` | One row per contact-event timepoint, with the event's classification and the threshold used. |
| `active_killing_summary_<immune_type>.csv` | Per-sample aggregates: number of active-killing timepoints, mean killing efficiency, active-killing rate. |
| `contact_events_<immune_type>.csv` | One row per contact event (start / end timepoint, duration, target track IDs). |
| `plots/combined_killing_efficiency_distribution.png` | Histogram of killing efficiency across all active-killing timepoints. |
