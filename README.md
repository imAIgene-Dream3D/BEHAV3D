# BEHAV3D

A Python package for analyzing cell behavior in fluorescent 3D imaging.

## Installation

### Quick Install (Recommended)

The easiest way to install BEHAV3D is using our automated installer that handles everything:

**Windows:**
1. Download this repository (Code > Download ZIP) and extract it
2. Double-click `install_behav3d.bat`
3. Follow the prompts

**macOS / Linux:**
```bash
cd /path/to/BEHAV3D
chmod +x install_behav3d.sh
./install_behav3d.sh
```

The installer will automatically:
- ✅ Install Miniforge (conda) if not found
- ✅ Detect your GPU (NVIDIA CUDA / Apple Silicon MPS / CPU-only)
- ✅ Create the conda environment with all dependencies
- ✅ Install PyTorch with the appropriate backend
- ✅ Install Cellpose for cell segmentation

### Manual Installation

If you prefer manual installation:

**Step 1: Install Miniforge**

If you don't have conda/mamba installed, download [Miniforge](https://github.com/conda-forge/miniforge/releases/latest/) (recommended, includes `mamba`) or [Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

> **Note:** We recommend using `mamba` instead of `conda` for faster dependency resolution. Miniforge includes `mamba` by default. All `mamba` commands below can be replaced with `conda` if using Miniconda/Anaconda.

**Step 2: Create the environment**

```bash
cd /path/to/BEHAV3D

# Create environment from yml file
mamba env create -f environment.yml
mamba activate behav3d

# Install Cellpose
pip install cellpose>=3.0

# Install PyTorch (choose ONE based on your system):

# First remove any preexisting pytorch installations:
mamba remove pytorch torchvision torchaudio

# CPU only (all platforms):
mamba install pytorch=2.4.1 torchvision=0.19.1 torchaudio=2.4.1 -c pytorch

# CUDA 12.1 (Windows/Linux with NVIDIA GPU):
mamba install pytorch=2.4.1 torchvision=0.19.1 torchaudio=2.4.1 pytorch-cuda=12.1 -c pytorch -c nvidia

# Install BEHAV3D package
pip install -e .

# Register Jupyter kernel
python -m ipykernel install --user --name=behav3d --display-name "behav3d"
```

***If pytorch installation is not working or you run into errors (OpenMP related), try this:***

**Option 1:** Install nomkl **before** installing pytorch (avoids conflicts with OpenMP)
```
# First switch BLAS backend to openblas (fixes igraph/nomkl conflict)
mamba install -c conda-forge igraph "blas=*=openblas"
# Then install nomkl
mamba install nomkl
```
**Option 2:** Remove cellpose and re-install it with conda:
```
pip uninstall cellpose -y
mamba install -c conda-forge cellpose
```

### Installation Options

```bash
# Check system info before installing
python install_behav3d.py --check

# Force CPU-only installation (even if GPU detected)
python install_behav3d.py --cpu-only

# Only install PyTorch/Cellpose (if environment already exists)
python install_behav3d.py --pytorch-only
```

## How to Run

BEHAV3D is run through Jupyter notebooks. We recommend using **Visual Studio Code**:

1. Install [VS Code](https://code.visualstudio.com/Download)
2. Install the Python and Jupyter extensions
3. Open the BEHAV3D folder
4. Open `notebooks/run_behav3d.ipynb`
5. Select kernel: Python Environments > behav3d

Alternatively, run in a web browser:
```bash
mamba activate behav3d
jupyter notebook notebooks/run_behav3d.ipynb
```

### Required files

[REQUIRED]
- A [metadata .csv](./configs/metadata.csv) that contains all samples that should be analyzed. To manage the metadata.csv in windows we recommend using https://www.moderncsv.com/ as it doesn't modify the format upon saving.
- Per sample, a raw microscopy image (.czi, .tiff, .zarr)

[Optional]
- Already segmented data for each sample
- Already tracked data for each sample
  
## Metadata.csv structure

### Required
| Column Name                    | Explanation                                                                                         |
|-------------------------------|-----------------------------------------------------------------------------------------------------|
| sample_name                   | Name to assign to the sample; used for naming output files/folders.                                |
| organoid_line                 | Name of organoid cell line used in the experiment (e.g. 10T, 162M, etc.). Analysis visual results will be split on organoid lines                                     |
| tcell_line                    | Name or type of T cell used (e.g. TEG/CAR-T/etc.). Analysis visual results will be split on tcell lines                                 |
| exp_nr                        | Experiment number for tracking and reproducibility.                                                |
| well                          | Well identifier on the experimental plate (e.g. `well01`).                                          |
| tcell_channel                 | Index of the fluorescence channel representing T cells in the raw image. (Channel indexing starts at 0)                                      |
| live_channel                  | Index of the fluorescence channel representing alive cells in the raw image. (Channel indexing starts at 0)             |
| dead_channel                  | Index of the fluorescence channel representing dead signal in the raw image. (Channel indexing starts at 0)                                  |
| contact_threshold             | Distance threshold (in `distance_unit`) to define T cell-organoid or T cell - T cell contact.                        |
| pixel_distance_xy             | Physical pixel spacing in the XY plane (usually in micrometers per pixel).                         |
| pixel_distance_z              | Physical spacing between Z slices (z-step), in same unit as `pixel_distance_xy`.                   |
| distance_unit                 | Unit for all provided spatial distances (e.g. `μm` for micrometers).                                        |
| time_interval                 | Time between consecutive frames (usually in seconds).                                              |
| time_unit                     | Unit of the time interval (e.g. `s` for seconds, 'm' for minutes, 'h' for hours).                                                  |
| raw_image_path                | Full path to the original CZI or raw imaging file.                                                 |

### Optional
| Column Name                    | Explanation                                                                                         |
|-------------------------------|-----------------------------------------------------------------------------------------------------|
| tcell_segments_image_path     | Path to the segmented T cell image (e.g. .zarr file).                                              |
| tcell_tracks_image_path       | Path to the image showing T cell tracks over time.                                                 |
| tcell_tracks_csv_path         | Path to the CSV file containing T cell tracking data.                                              |
| organoid_segments_image_path  | Path to the segmented organoid image (e.g. .zarr file).                                            |
| organoid_tracks_image_path    | Path to the image showing organoid tracks over time.                                               |
| organoid_tracks_csv_path      | Path to the CSV file containing organoid tracking data.                                            |

\
\
\\

## Cellpose Segmentation

To use **Cellpose** for 3D segmentation with pretrained models using the `run_behav3d.ipynb` notebook:

### Pretrained Models
You can find and download pretrained T cell and Organoid Cellpose models in [Zenodo](https://zenodo.org/records/18872978).


### Prerequisites

1.  **Metadata**: Ensure your `metadata.csv` is correctly populated with `sample_name`, `raw_image_path`, and cell type categories (Organoid, Immune, Other).
2.  **Image Format**: Input images should be converted to **Zarr** format (using the "Convert to Zarr" button in the notebook) for efficient multiprocessing.


### Workflow in `run_behav3d.ipynb`

#### 1. Configure Channel Labels
Before running segmentation, you must map your image channels to cell type labels.
-   Locate the **Configure Channel Labels** section.
-   Select **"Same labels for all samples"** or **"Per-sample configuration"**.
-   Assign the correct label (e.g., `tcell`, `organoid`, `dead`) to each channel index.
-   Click **"Save Channel Configuration"**.

#### 2. Run Cellpose Segmentation
There are separate panels for different cell categories:

##### **Organoid Segmentation**
-   **Pretrained Model**: Specify the path of your pretrained Cellpose model file and load it.
-   **Input Channels**: Select which channels (e.g., `organoid1`) to segment.
-   Click **"Apply"**.

##### **Immune and Other Cell Segmentation**
-   **Pretrained Model**: Specify the path of your pretrained Cellpose model file and load it.
-   **Input Channels**: Select which channel (e.g., `tcell1`) to segment.
-   Click **"Apply"**.


### Visualization and Verification

After segmentation, use the **Visualize the sample** button:
1.  Select the sample from the dropdown.
2.  Click **"Visualize the sample"**.
3.  This launches **Napari**: overlay segmented and original raw channels for easy inspection.

---

*Note: To train a new Cellpose model, refer to the `train_behav3d_cellpose.ipynb` notebook.*
