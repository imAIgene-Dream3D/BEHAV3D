# Historical experiment profiles for assistant guidance

These profiles are evidence from completed or previously configured BEHAV3D
experiments. They are **not defaults** and must not be copied solely because a
new experiment has a similarly named cell type.

The assistant should use them only when a researcher asks for a previous
example, a similar experiment, or the values used in an earlier analysis. A
useful answer must:

1. Name the historical profile and the source of the value.
2. Compare acquisition resolution, frame interval, object scale, motion, signal
   quality, and the selected method with the current experiment.
3. Explain what must be measured or previewed before adapting the value.
4. Never issue form-edit actions from a historical value alone.

Source priority is:

- live metadata or the experiment CSV for acquisition facts;
- the experiment README for biological intent and operational definitions;
- the YAML for saved configuration;
- discovered output files for evidence that a module actually ran.

When sources disagree, state the disagreement. Do not silently choose the more
convenient value. Older YAML files may use legacy labels or units, so map every
setting to the current live control before proposing it.

## IVM HIV T-cell motility

**Evidence:** `IVM_HIV/metadata.csv`,
`IVM_HIV/behav3d_parameters.yml`, and
`IVM_HIV/README_BEHAV3D_IVM_HIV.md`.

- Three intravital movies contain infected and uninfected T-cell populations.
- The CSV records 1.15 um XY pixels, 4.0 um Z spacing, and 15 s frames for all
  three samples. The README describes one sample as 10 s, which conflicts with
  the CSV; use the loaded metadata unless the researcher confirms otherwise.
- The saved segmentation configuration uses ConvPaint Probability + Watershed,
  VGG features, multi-channel input, three examples per sample, stack
  normalization, and a shared classifier.
- Historical instance settings differ by population: HIV cells used mask/seed
  probabilities 0.85/0.89, EDT 2.74, and minimum size 31; CMTMR cells used
  0.60/0.90, EDT 1.87, and minimum size 10.
- Both populations used btrack with global optimization. Saved maximum search
  radii were 12 and 10 physical-distance units. Optimizer distance threshold was
  26 for both; time thresholds were 6 and 4 frames.
- Tracks were filtered to a fixed 15-frame window. A shared three-state HMM used
  cumulative displacement, displacement, and speed; trajectory clustering also
  used 15 frames, three clusters, and average linkage.

Use this profile to illustrate why motile T-cell settings can differ even
within one movie. Calibrate the search radius from speed and cadence, and
calibrate optimizer thresholds from observed spatial and missing-frame gaps.

## CD4/CD8 active killing against 13T

**Evidence:** `cd4cd8/metadata_cd8cd4.csv`,
`cd4cd8/behav3d_parameters_clean_cd4cd8.yml`, and
`cd4cd8/README_BEHAV3D_CD4CD8_13T_ActiveKilling.md`.

- Four movies use 1.01 um XY, 1.05 um Z, and 2 min frames.
- Blue and green acquisition slots swap CD4/CD8 identity between paired wells;
  the biological identity is stored in metadata rather than inferred from the
  channel name.
- The saved method is APOC Mask + EDT/Watershed. The YAML contains more than one
  post-processing snapshot, so its exact EDT/minimum-size numbers are historical
  tuning evidence, not a reproducible final preset.
- The organoid classifier used all four channels as contextual input; blue,
  green, and dead models used their own channels. This is a study-specific
  example of shared channels being useful only when the researcher confirms
  that they carry target or boundary context.
- The organoid used propagation. Blue and green T cells used matched btrack
  settings: maximum search radius 100, optimizer distance threshold 60, time
  threshold 3 frames, and the false-positive/initiation/termination/linking
  hypotheses.
- The saved Active Killing definition used a 7-frame observation window, a
  relative 2x rise in dead-mask percentage, and at least 3 contact frames.

Use this profile for dye-swap study design, matched settings across biological
groups, and an example of a cadence-dependent Active Killing definition. Do not
present 100, 60, 3, 7, or 2x as general T-cell defaults.

## Multi-organoid safety profiling

**Evidence:** `safety_profiling/metadata.csv`,
`safety_profiling/behav3d_parameters_multiorganoid_clean.yml`, and
`safety_profiling/README_BEHAV3D_SafetyProfiling_Exp010.md`.

- Four 3D movies compare tumor 27T, healthy MDO, and TEG effectors at 1.77 um
  isotropic spacing and 2 min frames.
- The main comparison is paired 27T versus MDO within the combined wells. The
  controls have only one or two wells, so broader comparisons are descriptive.
- APOC was used, but the saved configuration contains multiple classifier
  tuning rounds. Exact classifier features and post-processing thresholds must
  not be called the final settings without result/model provenance.
- 27T and MDO used propagation; TEG used btrack.
- Feature extraction used symmetric organoid death thresholds and computed
  invasiveness for TEG, not for the organoid targets.
- The experiment README defines Active Killing as a relative 1.5x rise over a
  5-frame window after at least one frame of contact. This describes the study
  definition; it is not proof that the analysis was run.

Use this profile for tumor-versus-healthy safety comparisons, symmetric target
configuration, and the distinction between configured analysis and discovered
results.

## Near-static pancreatic islet calcium reporter

**Evidence:** `Explants_Calcium/metadata_4samples.csv`,
`Explants_Calcium/behav3d_parameters.yml`, and
`Explants_Calcium/README_Calcium_panet islets.md`.

- Four samples compare control and T-cell co-culture conditions across two
  experiments. Acquisition is 0.33 um XY, 2.0 um Z, 5 s frames, and 32 frames.
- Segmentation was generated externally with Cellpose-SAM and imported.
- Near-static cells with intermittent reporter visibility used Reporter
  Propagation: detections were pooled across time, grouped by spatial overlap,
  and one canonical mask was propagated through the movie.
- The README records a historical 100-voxel noise cutoff and 10% overlap
  grouping rule. Those are example values for this resolution and object scale,
  not universal Reporter Propagation defaults.
- Filtering retained the full 32-frame duration.
- The saved five-state HMM used one reporter-intensity fold-change feature,
  `q75_mean_intensity_ch0_fold_change`, with a smoothing window of 1. The
  five-cluster trajectory analysis used the full 32-frame sequence and average
  linkage.

Use this profile to explain why intensity can be the primary behavioral signal
when propagation makes movement and morphology nearly constant. Do not use its
HMM feature set for motility experiments.

## Method-level lessons from the integrated Methods description

- Method choice follows signal quality, object morphology, compute, and whether
  sparse labels or a zero-shot method are appropriate.
- Pixel classifiers produce class probabilities and require instance recovery;
  Cellpose and Cellpose-SAM directly produce instances.
- Pixel-classifier annotations should span samples and intensity ranges.
- Spatial thresholds must be interpreted with image spacing, and temporal
  thresholds with frame cadence.
- Propagation is appropriate for large, slow, deformable objects; Reporter
  Propagation is for near-static objects with intermittent detection; btrack is
  for motile objects and may use a global hypothesis optimizer after Step 1 is
  reviewed.
- Feature and analysis choices must follow the biological readout. Motility
  features are appropriate for the IVM profile, while reporter intensity is the
  informative signal in the calcium profile.
