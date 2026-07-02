# Installation

BEHAV3D EXPLORER runs on Windows, macOS and Linux. The easiest way to install BEHAV3D EXPLORER is using our automated installer that handles everything, but manual installation steps are also provided.


## 📁 Files Overview

 File | Description |
|------|-------------|
| `install_behav3d.py` | Main Python installation script (cross-platform) |
| `install_behav3d_windows.bat` | Windows batch wrapper (double-click to run) |
| `install_behav3d_linux.sh` | macOS/Linux shell wrapper |
| `install_behav3d_macOS.command` | macOS double-click launcher |
| `environment.yml` | Conda environment specification |
| `requirements.txt` | Pip requirements (alternative to conda) |

## 🚀 Quick Start

### Windows
**Option 1 (Recommended):** Double-click `install_behav3d_windows.bat`

**Option 2:** Open Command Prompt and run:
```bat
cd path\to\BEHAV3D\installation
python install_behav3d.py
```

### macOS
**Option 1:** Double-click `install_behav3d_macOS.command`

**Option 2:** Open Terminal and run:
```bash
cd path/to/BEHAV3D/installation
chmod +x install_behav3d_linux.sh
./install_behav3d.sh
```

### Linux
```bash
cd path/to/BEHAV3D/installation
chmod +x install_behav3d_linux.sh
./install_behav3d.sh
```

What the installer does automatically:

- ✅ Install Miniforge (conda) if not found
- ✅ Detect your GPU (NVIDIA CUDA / Apple Silicon MPS / CPU-only)
- ✅ Create the conda environment with all dependencies
- ✅ Install PyTorch with the appropriate backend
- ✅ Install Cellpose for cell segmentation


## ⚙️ Installation Options

The installer supports several command-line options:

| Option | Description |
|--------|-------------|
| `-n NAME` | Set custom environment name (default: `behav3d`) |
| `--cpu-only` | Force CPU-only installation (no GPU support) |
| `--reinstall` | Remove and reinstall if environment already exists |
| `--keep-existing` | Keep existing environment, only update PyTorch/Cellpose |
| `--pytorch-only` | Only install PyTorch/Cellpose (assumes env exists) |
| `--check` | Check system information without installing |

### Examples

```bash
# Full installation with default settings
python install_behav3d.py

# Install with custom environment name
python install_behav3d.py -n my_behav3d_env

# Force CPU-only installation
python install_behav3d.py --cpu-only

# Check system info (CUDA, GPU, etc.)
python install_behav3d.py --check

# Reinstall from scratch
python install_behav3d.py --reinstall
```


## 📋 Requirements

- **Python**: 3.12 (recommended)
- **Miniforge** (recommended, includes `mamba`) or Conda/Miniconda
- **CUDA**: 12.1+ (optional, for GPU acceleration)

### Package Manager Priority

The installer automatically detects and prefers the fastest available package manager:

1. **mamba** (fastest, included with Miniforge/Mambaforge)
2. **micromamba** (lightweight standalone alternative)
3. **conda** (fallback)

> All `mamba` commands can be replaced with `conda` if using Miniconda/Anaconda.

### Supported Platforms
- ✅ Windows 10/11 (x64)
- ✅ macOS (Intel and Apple Silicon)
- ✅ Linux (x64)


## 🔧 Manual Installation

If the installer fails or you want full control:

### Step 1: Install Miniforge

Download [Miniforge](https://github.com/conda-forge/miniforge/releases/latest/) (recommended, includes `mamba`) or [Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

```{tip}
We recommend using `mamba` instead of `conda` for faster dependency resolution. Miniforge includes `mamba` by default. All `mamba` commands below can be replaced with `conda` if using Miniconda/Anaconda.
```

### Step 2: Create the env

```bash
# Create environment from yml file
cd /path/to/BEHAV3D
mamba env create -f environment.yml
# Activate environment before next stps
mamba activate behav3d
```

### Step 3: Install Cellpose and PyTorch

```bash
# Cellpose
pip install cellpose>=3.0

# First remove any preexisting torch installs
mamba remove pytorch torchvision torchaudio

# Then pick ONE of:
# CPU only (all platforms):
mamba install pytorch=2.4.1 torchvision=0.19.1 torchaudio=2.4.1 -c pytorch

# CUDA 12.1 (Windows/Linux with NVIDIA GPU):
mamba install pytorch=2.4.1 torchvision=0.19.1 torchaudio=2.4.1 pytorch-cuda=12.1 -c pytorch -c nvidia
```

### Step 4: Install BEHAV3D EXPLORER itself

```bash
pip install -e .
# Register Jupyter kernel
python -m ipykernel install --user --name=behav3d --display-name "behav3d"
```


## ❓ Troubleshooting

### "CUDA not detected" although you have an NVIDIA card
1. Make sure NVIDIA drivers are installed (run `nvidia-smi` from a regular terminal — it should print your card).
2. Run `python install_behav3d.py --check` to see what the installer detects.
3. If still missing, reinstall with `python install_behav3d.py --reinstall`.

### Conda Not Found
The installer will automatically download and install Miniforge if conda is not found.

### Environment already exists
```bash
python install_behav3d.py --reinstall       # wipe and recreate
python install_behav3d.py --keep-existing   # update PyTorch/Cellpose only
python install_behav3d.py -n behav3d_v2     # create alongside under a new name
```

### Import errors after installation
```bash
mamba activate behav3d
pip install --upgrade behav3d
```

### OpenMP/MKL Conflicts (nomkl fails to install)
If you see `LibMambaUnsatisfiableError` when installing `nomkl` (because `igraph` requires `blas==mkl`):
```bash
# Switch BLAS backend to openblas first
mamba install -c conda-forge igraph "blas=*=openblas"
# Then install nomkl
mamba install nomkl
```

### Cellpose import errors
If `pip install cellpose` produces import errors:

```bash
pip uninstall cellpose -y
mamba install -c conda-forge cellpose
```

The automated installer handles this automatically.


## 📞 Support

If you encounter issues:
1. Check the troubleshooting section above
2. Run `python install_behav3d.py --check` and include the output
3. Open an issue on the BEHAV3D GitHub repository


## Verifying your install/ How to Run

### Napari GUI (Recommended)

BEHAV3D EXPLORER provides a full graphical interface through the napari viewer:

**Quick launch:**
1. Double-click the launcher for your platform:
   - Windows: `napari/run_behav3d_windows.bat`
   - macOS: `napari/run_behav3d_macOS.command`
   - Linux: `napari/run_behav3d_linux.sh` (run from a terminal: `chmod +x napari/run_behav3d_linux.sh && ./napari/run_behav3d_linux.sh`)
2. The BEHAV3D EXPLORER plugin will open automatically inside napari

> **Note:** The launcher reads its environment configuration from `napari/.config/behav3d_env.json`, which is generated by the installer.

**Manual launch:**
```bash
mamba activate behav3d
napari
```

Then open the BEHAV3D EXPLORER plugin from the **Plugins** menu. You should see a dock widget with seven tabs:

📋 Data Preparation · 👁 Visualization · 🦠 Segmentation · 📍 Tracking · 🧪 Feature Extraction · 🧹 Filtering · 📊 Analysis

…plus the **Processing Queue** panel at the bottom. If all seven tabs render, you're good to go.

### Jupyter Notebook (Alternative)

> **Note:** The Jupyter notebook provides an alternative, script-based interface for advanced users who prefer more control over individual pipeline steps.

1. Install [VS Code](https://code.visualstudio.com/Download) with the Python and Jupyter extensions
2. Open the BEHAV3D folder
3. Open `notebooks/run_behav3d.ipynb`
4. Select kernel: Python Environments > behav3d

Or in a web browser:
```bash
mamba activate behav3d
jupyter notebook notebooks/run_behav3d.ipynb
```

## Required input files

To start a pipeline you need:

- A populated `metadata.csv` (see [Data Preparation](preprocessing/data_preparation) for the full schema).
- Per sample, a raw microscopy image: `.czi`, `.tiff`, `.zarr`, `.lif` or `.ims`.

```{tip}
On Windows we recommend [Modern CSV](https://www.moderncsv.com/) if you need to open `metadata.csv` outisde BEHAV3D EXPLORER, it doesn't silently reformat values the way Excel does.
```

Optional pre-existing data you can plug in via the import options:

- Already-segmented label images (see [Import segmentation](processing/segmentation/import)).
- Already-tracked label images + track CSVs (see [Import tracking](processing/tracking/import)).

## Next step

→ [Getting Started](getting_started)
