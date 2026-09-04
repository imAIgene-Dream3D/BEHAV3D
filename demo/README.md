# BEHAV3D Explorer — online demo (maintainer guide)

This folder lets anyone try the **real** BEHAV3D Explorer GUI in a browser, with no installation,
for free, on Google Colab.

```
      Colab VM (free, ~12.7 GB RAM, optional T4)
 ┌───────────────────────────────────────────────────────────┐
 │  napari (PyQt5) + BEHAV3DWidget                           │
 │        │ renders into                                     │
 │  Xvfb :99  (software OpenGL / llvmpipe)                   │
 │        │ captured by                                      │
 │  x11vnc :5900  ──►  websockify/noVNC :6080                │
 └────────────────────────────────┬──────────────────────────┘
                                  │ Colab port proxy (or a Cloudflare quick tunnel)
                                  ▼
                     the visitor's browser tab
```

The visitor sees the identical GUI to a local install — the same window, menus, dock widget and
layer list. Nothing in `behav3d/` is reimplemented; only `$DISPLAY` differs. napari is started
through the repository's own launcher, `napari/.config/launch_napari.py --internal`.

### Files here

| File | What it is |
|---|---|
| `colab/BEHAV3D_Explorer_Colab.ipynb` | What visitors open. Three steps, "Run all". |
| `colab/colab_setup.py` | All the machinery: apt, environment, data, display stack, launch, tunnel fallback, health checks. |
| `colab/prepare_demo.py` | Rewrites the absolute paths inside a demo bundle for the machine it landed on. |
| `colab/apt_packages.txt` | The headless-Qt / OpenGL / VNC package list, shared by Colab and the Docker test. |
| `build_env.sh` | Builds the prebuilt environment tarball (and can run a local dry-run of the GUI). |

**Cost: $0.** Colab's free tier runs the demo; Zenodo and Hugging Face dataset repos host the
files for free. Visitors need a Google account (unavoidable at zero cost).

---

## Step 1 — Build the environment tarball

Solving `installation/environment.yml` inside Colab takes 10+ minutes and breaks whenever
conda-forge moves. Solve it once here instead, inside Ubuntu 22.04 (Colab's own base), and ship the
result as a [conda-pack](https://conda.github.io/conda-pack/) archive that restores in ~2 minutes.

```bash
./demo/build_env.sh
```

* Needs Docker. Takes 20–40 min the first time.
* Produces `dist/behav3d-env.tar.gz` (~2.5 GB for the CPU build).
* `--cuda` builds against CUDA 12.8 instead (matches `CUDA_VERSION` in
  `installation/install_behav3d.py`). It is ~5 GB, so visitors wait longer for a GPU they may not
  get. **Recommendation: ship the CPU build**, and keep a CUDA build around for workshop days.
* The build fails loudly if `napari`, `torch`, `cellpose`, `pyopencl` or `qtpy` cannot be imported —
  better here than in front of a visitor.

## Step 2 — Host the tarball

A free Hugging Face account can host **dataset** repos (only *Spaces* need PRO), files up to 50 GB,
on a fast CDN. GitHub Releases cap at 2 GB per file, which this tarball exceeds.

```bash
pip install huggingface_hub
huggingface-cli login
huggingface-cli repo create behav3d-demo-env --type dataset
huggingface-cli upload <your-org>/behav3d-demo-env dist/behav3d-env.tar.gz \
    behav3d-env.tar.gz --repo-type dataset
```

(`huggingface_hub` 0.34+ renames this command to `hf`, e.g. `hf upload ...`; the CI workflow uses that newer name.)

The direct URL is then:

```
https://huggingface.co/datasets/<your-org>/behav3d-demo-env/resolve/main/behav3d-env.tar.gz
```

## Step 3 — Assemble the demo bundle

Run BEHAV3D locally once on the demo sample so the heavy steps are already done, keeping
**everything under one common folder** (the path rewriter finds the old root as the common parent
of all recorded paths):

```
behav3d_demo/
├── raw/                     # the ~400 MB images
├── metadata.csv             # produced by the Data Preparation tab
└── output/                  # segmentation, tracking, features, analysis results
    └── behav3d_parameters.yml
```

Trim it so it stays pleasant on a free runtime: one sample, cropped in Z/T if needed, target
≤ 500 MB raw. Then:

```bash
tar -czf behav3d_demo.tar.gz behav3d_demo/
```

Sanity-check the rewriter before uploading anything:

```bash
python demo/colab/prepare_demo.py --root /path/to/behav3d_demo --dry-run
```

It should report the old root, the number of paths it would rewrite, and **no missing files**.

## Step 4 — Publish the bundle on Zenodo

Zenodo is free, gives 50 GB per record and a citable DOI, and you already use it for the Cellpose
models (record 18872978). Create a new record, upload `behav3d_demo.tar.gz`, publish, then take the
direct file URL:

```
https://zenodo.org/records/<RECORD_ID>/files/behav3d_demo.tar.gz?download=1
```

## Step 5 — Fill in the two URLs

In `demo/colab/colab_setup.py`, replace the two placeholders:

```python
ENV_URL  = "https://huggingface.co/datasets/<your-org>/behav3d-demo-env/resolve/main/behav3d-env.tar.gz"
DEMO_URL = "https://zenodo.org/records/<RECORD_ID>/files/behav3d_demo.tar.gz?download=1"
```

Both are also overridable at runtime via `BEHAV3D_ENV_URL` / `BEHAV3D_DEMO_URL`, which is how the
notebook lets visitors point at their own data. The bootstrap refuses to run while a placeholder is
still in place, with a message saying which step to go back to.

## Step 6 — Dry run locally (do this before touching Colab)

This catches every Qt / xcb / OpenGL problem on your own machine, in the same Ubuntu base Colab
uses, in about a minute:

```bash
./demo/build_env.sh --test
```

Then open <http://localhost:6080/vnc.html>. The napari window with the BEHAV3D Explorer dock widget
must appear. If it does not, read `/tmp/behav3d_logs/napari.log` inside the container.

## Step 7 — Test on Colab

Push the branch, then open:

```
https://colab.research.google.com/github/imAIgene-Dream3D/BEHAV3D/blob/main/demo/colab/BEHAV3D_Explorer_Colab.ipynb
```

Check, in order:

1. Run all completes in under ~5 minutes on a cold runtime.
2. The new tab shows the GUI and it is **interactive** — drag a slider, open a menu, switch tabs.
3. The two paths printed by Step 2 load in the Data Preparation tab.
4. Repeat once with a GPU runtime and once CPU-only (`Runtime ▸ Change runtime type`).
5. Run the `start_cloudflared()` cell and confirm that route works too.
6. Watch RAM in the Colab resource panel with the demo loaded — stay under ~10 GB.

## Step 8 — Add the badge

```markdown
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/imAIgene-Dream3D/BEHAV3D/blob/main/demo/colab/BEHAV3D_Explorer_Colab.ipynb)
```

## Step 9 — Maintenance

Rebuild and re-upload the tarball whenever `installation/environment.yml` changes; nothing else
moves. `.github/workflows/build-demo-env.yml` does this on demand and pushes the result to the
Hugging Face dataset repo (needs an `HF_TOKEN` repository secret with write access).

---

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `Could not load the Qt platform plugin "xcb"` | a missing `libxcb-*` | add it to `colab/apt_packages.txt`; run napari with `QT_DEBUG_PLUGINS=1` to see which one |
| Tab opens but stays black / never connects | Colab's proxy dropped WebSockets | run the `start_cloudflared()` cell |
| `no python at /opt/behav3d/bin/python` | the tarball is not a conda-pack archive, or it was built for a different base | rebuild with `demo/build_env.sh` |
| napari exits immediately | see `colab_setup.tail("napari", 60)` | usually a missing library or a bad `PYTHONPATH` |
| Session dies while loading the images | out of RAM | crop the demo bundle further |
| Everything is very slow in 3D | software OpenGL (llvmpipe), expected | stay in 2D; a GPU runtime does not help, since the X display has no GPU |

## Other hosting routes considered

* **GitHub Codespaces** — works well with the `desktop-lite` devcontainer feature (8 GB, no GPU,
  spends the *visitor's* free 120 core-hours). A good second front door; the same
  `colab_setup.py --local` entry point drives it.
* **Hugging Face Docker Space** — the only genuinely no-login option, but creating a Docker Space
  now requires PRO ($9/month) and a Space is a *single shared container*: two simultaneous visitors
  would fight over one mouse.
* **mybinder.org** — free and login-free, but capped at 2 GB RAM; napari + torch + the dataset do
  not fit.
* **usegalaxy.eu** — napari is already an interactive tool there, sessions are per-user and last up
  to 30 days. Free and durable, but it needs the Galaxy team to accept a BEHAV3D wrapper.
