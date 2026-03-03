#!/bin/bash
# BEHAV3D Napari Launcher (Linux)
# Run: chmod +x run_behav3d.sh && ./run_behav3d.sh

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
echo ""
echo "  ========================================"
echo "    BEHAV3D - Launching Napari"
echo "  ========================================"
echo ""

python3 "$SCRIPT_DIR/launch_napari.py" || python "$SCRIPT_DIR/launch_napari.py"
