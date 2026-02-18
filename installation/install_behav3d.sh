#!/bin/bash
# =============================================================================
# BEHAV3D Installer - Unix Shell Script (macOS/Linux)
# =============================================================================
# This script wraps the Python installer for easy execution.
#
# Usage:
#   ./install_behav3d.sh                    # Full installation (prompts for env name)
#   ./install_behav3d.sh -n myenv           # Install with custom environment name
#   ./install_behav3d.sh --cpu-only         # Force CPU-only
#   ./install_behav3d.sh --reinstall        # Remove and reinstall if exists
#   ./install_behav3d.sh --keep-existing    # Keep existing, update PyTorch/Cellpose
#   ./install_behav3d.sh --check            # Check system info only
#
# Examples:
#   ./install_behav3d.sh -n my_behav3d_env
#   ./install_behav3d.sh -n test_env --reinstall --cpu-only
# =============================================================================

set -e

echo ""
echo "============================================================"
echo "              BEHAV3D Installation Script"
echo "============================================================"
echo ""

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Check if Python is available
if command -v python3 &> /dev/null; then
    PYTHON_CMD="python3"
elif command -v python &> /dev/null; then
    PYTHON_CMD="python"
else
    echo "[ERROR] Python not found."
    echo ""
    echo "Please install Python 3.10+ or Miniforge:"
    echo "  macOS:  brew install python3"
    echo "          or: https://github.com/conda-forge/miniforge"
    echo "  Linux:  sudo apt install python3 (Debian/Ubuntu)"
    echo "          or: https://github.com/conda-forge/miniforge"
    exit 1
fi

# Initialize variables
ENV_NAME="behav3d"
ARGS="$@"

# If no arguments provided, prompt for environment name
if [ -z "$ARGS" ]; then
    echo "============================================================"
    echo "              Environment Configuration"
    echo "============================================================"
    echo ""
    echo "Enter the name for the conda environment."
    echo "Press ENTER to use the default name: behav3d"
    echo ""
    read -p "Environment name [behav3d]: " USER_ENV_NAME
    
    if [ -n "$USER_ENV_NAME" ]; then
        ENV_NAME="$USER_ENV_NAME"
    fi
    
    ARGS="-n $ENV_NAME"
    echo ""
    echo "[INFO] Using environment name: $ENV_NAME"
else
    # Try to extract environment name from args
    if echo "$ARGS" | grep -qE "\-n|\-\-name"; then
        ENV_NAME=$(echo "$ARGS" | sed -n 's/.*\(-n\|--name\)[[:space:]]*\([^[:space:]]*\).*/\2/p')
        if [ -z "$ENV_NAME" ]; then
            ENV_NAME="behav3d"
        fi
    fi
fi

echo ""

# Run the Python installer
$PYTHON_CMD "$SCRIPT_DIR/install_behav3d.py" $ARGS

exit_code=$?

if [ $exit_code -ne 0 ]; then
    echo ""
    echo "[ERROR] Installation failed with error code $exit_code"
    exit $exit_code
fi

echo ""
