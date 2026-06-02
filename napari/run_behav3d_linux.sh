#!/bin/bash
# BEHAV3D Napari Launcher (Linux)
# Run: chmod +x run_behav3d_linux.sh && ./run_behav3d_linux.sh

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
echo ""
echo "  ========================================"
echo "    BEHAV3D - Launching Napari"
echo "  ========================================"
echo ""

# Suppress PyOpenCL compiler cache warnings
export PYOPENCL_NO_CACHE=1
# export PYOPENCL_COMPILER_OUTPUT=0

python3 "$SCRIPT_DIR/.config/launch_napari.py" || python "$SCRIPT_DIR/.config/launch_napari.py"
