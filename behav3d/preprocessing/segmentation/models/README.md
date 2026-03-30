# Cellpose Segmentation

This guide explains how to use **Cellpose** for 3D segmentation with pretrained models using the `run_behav3d.ipynb` notebook.

## Prerequisites

1.  **Metadata**: Ensure your `metadata.csv` is correctly populated with `sample_name`, `raw_image_path`, and cell type categories (Organoid, Immune, Other).
2.  **Image Format**: Input images should be converted to **Zarr** format (using the "Convert to Zarr" button in the notebook) for efficient multiprocessing.


## Workflow in `run_behav3d.ipynb`

### 1. Configure Channel Labels
Before running segmentation, you must map your image channels to cell type labels.
-   Locate the **Configure Channel Labels** section.
-   Select **"Same labels for all samples"** or **"Per-sample configuration"**.
-   Assign the correct label (e.g., `tcell`, `organoid`, `dead`) to each channel index.
-   Click **"Save Channel Configuration"**.

### 2. Run Cellpose Segmentation
There are separate panels for different cell categories:

#### **Organoid Segmentation**
-   **Pretrained Model**: Specify the path of your pretrained Cellpose model file and load it.
-   **Input Channels**: Select which channels (e.g., `organoid1`) to segment.
-   Click **"Apply"**.

#### **Immune and Other Cell Segmentation**
-   **Pretrained Model**: Specify the path of your pretrained Cellpose model file and load it.
-   **Input Channels**: Select which channel (e.g., `tcell1`) to segment.
-   Click **"Apply"**.


## Visualization and Verification

After segmentation, use the **Visualize the sample** button:
1.  Select the sample from the dropdown.
2.  Click **"Visualize the sample"**.
3.  This launches **Napari**: overlay segmented and original raw channels for easy inspection.

---

*Note: To train a new Cellpose model, refer to the `train_behav3d_cellpose.ipynb` notebook.*
