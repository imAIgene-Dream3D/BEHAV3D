# BEHAV3D Installation Guide

This folder contains all the necessary files to install BEHAV3D and its dependencies.

## 📁 Files Overview

| File | Description |
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
```cmd
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

1. **mamba** (fastest — included with Miniforge/Mambaforge)
2. **micromamba** (lightweight standalone alternative)
3. **conda** (fallback)

> All `mamba` commands can be replaced with `conda` if using Miniconda/Anaconda.

### Supported Platforms
- ✅ Windows 10/11 (x64)
- ✅ macOS (Intel and Apple Silicon)
- ✅ Linux (x64)

## 🔧 Manual Installation

If the automatic installer doesn't work, you can install manually:

### Using Mamba/Conda (Recommended)
```bash
# Create environment from file (use mamba for speed, or conda)
mamba env create -f environment.yml

# Activate environment
mamba activate behav3d

# Install PyTorch with CUDA support (if GPU available)
python install_behav3d.py --pytorch-only

# Note: If installing PyTorch/Cellpose manually without the installer script:
# pip install napari-convpaint
```

### Using Pip
```bash
# Create virtual environment
python -m venv behav3d_env

# Activate environment
# Windows:
behav3d_env\Scripts\activate
# macOS/Linux:
source behav3d_env/bin/activate

# Install PyTorch first (with CUDA if available)
# See https://pytorch.org/get-started/locally/ for the correct command

# Install requirements
pip install -r requirements.txt

# Install ConvPaint
pip install napari-convpaint
```

## 🎯 What Gets Installed

The installer sets up:

1. **Core Scientific Computing**: NumPy, Pandas, SciPy, Dask
2. **Visualization**: Matplotlib, Seaborn, Plotly
3. **Image Processing**: scikit-image, tifffile, imageio, Pillow
4. **Deep Learning**: PyTorch (with CUDA if available), Cellpose
5. **Tracking**: TrackPy, LAP solver
6. **Napari**: For visualization and interaction
7. **Jupyter**: Notebook support

## ❓ Troubleshooting

### CUDA/GPU Not Detected
1. Ensure NVIDIA drivers are installed
2. Run `python install_behav3d.py --check` to see detected hardware
3. Try reinstalling with `--reinstall` flag

### Conda Not Found
The installer will automatically download and install Miniforge if conda is not found.

### Environment Already Exists
Use one of these options:
- `--reinstall` to remove and recreate
- `--keep-existing` to update only PyTorch/Cellpose
- `-n new_name` to create with a different name

### Import Errors After Installation
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
The automated installer handles this automatically.

## 📞 Support

If you encounter issues:
1. Check the troubleshooting section above
2. Run `python install_behav3d.py --check` and include the output
3. Open an issue on the BEHAV3D GitHub repository
