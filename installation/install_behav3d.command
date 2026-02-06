#!/bin/bash
# =============================================================================
# BEHAV3D Installer - macOS Double-Click Script (.command)
# =============================================================================
# This script can be double-clicked in Finder to run the installer.
# It will search for Python/Conda and prompt for configuration.
#
# Options (if run from terminal):
#   -n NAME           Set custom environment name (default: behav3d)
#   --cpu-only        Force CPU-only installation
#   --reinstall       Remove and reinstall if environment exists
#   --keep-existing   Keep existing environment, only update PyTorch/Cellpose
#   --check           Only check system information
#
# Examples:
#   ./install_behav3d.command
#   ./install_behav3d.command -n my_behav3d_env
#   ./install_behav3d.command --cpu-only
# =============================================================================

# Change to the directory where this script is located
cd "$(dirname "$0")"
SCRIPT_DIR="$(pwd)"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

print_header() {
    echo ""
    echo -e "${CYAN}============================================================${NC}"
    echo -e "${CYAN}              $1${NC}"
    echo -e "${CYAN}============================================================${NC}"
    echo ""
}

print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_header "BEHAV3D Installation Script"

# Initialize variables
PYTHON_CMD=""
ENV_NAME="behav3d"
ARGS="$@"

# =============================================================================
# Find Python/Conda
# =============================================================================

print_info "Searching for Python/Conda installation..."
echo ""

# Check common Miniforge/Conda locations
CONDA_LOCATIONS=(
    "$HOME/miniforge3"
    "$HOME/mambaforge"
    "$HOME/miniconda3"
    "$HOME/anaconda3"
    "/opt/homebrew/Caskroom/miniforge/base"
    "/opt/homebrew/Caskroom/miniconda/base"
    "/usr/local/Caskroom/miniforge/base"
    "/usr/local/Caskroom/miniconda/base"
)

for location in "${CONDA_LOCATIONS[@]}"; do
    if [ -f "$location/bin/python" ]; then
        PYTHON_CMD="$location/bin/python"
        print_success "Conda Python found at: $location"
        break
    fi
done

# If not found, check PATH
if [ -z "$PYTHON_CMD" ]; then
    if command -v python3 &> /dev/null; then
        PYTHON_CMD="python3"
        print_success "Python3 found in PATH"
    elif command -v python &> /dev/null; then
        PYTHON_CMD="python"
        print_success "Python found in PATH"
    fi
fi

# If still not found, offer to install Miniforge
if [ -z "$PYTHON_CMD" ]; then
    print_header "Python/Conda Not Found"
    
    echo "BEHAV3D requires Python. The easiest solution is to install"
    echo "Miniforge, which includes Python and conda package manager."
    echo ""
    read -p "Would you like to download and install Miniforge? [Y/n]: " INSTALL_CHOICE
    
    if [[ "$INSTALL_CHOICE" =~ ^[Nn]$ ]]; then
        echo ""
        print_error "Installation cancelled. Please install Python or Miniforge manually:"
        echo "  https://github.com/conda-forge/miniforge/releases/latest/"
        echo ""
        echo "Press any key to close..."
        read -n 1
        exit 1
    fi
    
    # Detect architecture
    ARCH=$(uname -m)
    if [ "$ARCH" = "arm64" ]; then
        MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh"
    else
        MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh"
    fi
    
    print_info "Downloading Miniforge installer..."
    INSTALLER_PATH="/tmp/miniforge_installer.sh"
    
    if command -v curl &> /dev/null; then
        curl -L -o "$INSTALLER_PATH" "$MINIFORGE_URL"
    elif command -v wget &> /dev/null; then
        wget -O "$INSTALLER_PATH" "$MINIFORGE_URL"
    else
        print_error "Neither curl nor wget found. Cannot download installer."
        echo "Press any key to close..."
        read -n 1
        exit 1
    fi
    
    if [ ! -f "$INSTALLER_PATH" ]; then
        print_error "Failed to download Miniforge installer."
        echo "Press any key to close..."
        read -n 1
        exit 1
    fi
    
    print_info "Installing Miniforge (this may take a few minutes)..."
    bash "$INSTALLER_PATH" -b -p "$HOME/miniforge3"
    
    if [ -f "$HOME/miniforge3/bin/python" ]; then
        PYTHON_CMD="$HOME/miniforge3/bin/python"
        print_success "Miniforge installed successfully!"
        rm -f "$INSTALLER_PATH"
        
        # Initialize conda for the shell
        "$HOME/miniforge3/bin/conda" init bash zsh 2>/dev/null
        
        print_warning "Note: You may need to restart your terminal after installation"
        print_warning "for conda commands to work properly."
        echo ""
    else
        print_error "Miniforge installation may have failed."
        echo "Please try installing manually from:"
        echo "  https://github.com/conda-forge/miniforge/releases/latest/"
        echo ""
        echo "Press any key to close..."
        read -n 1
        exit 1
    fi
fi

# =============================================================================
# Verify Python works
# =============================================================================

echo ""
print_info "Using Python: $PYTHON_CMD"
$PYTHON_CMD --version

if [ $? -ne 0 ]; then
    print_error "Python found but failed to execute."
    echo "Press any key to close..."
    read -n 1
    exit 1
fi

# =============================================================================
# Prompt for environment name if no arguments provided
# =============================================================================

if [ -z "$ARGS" ]; then
    print_header "Environment Configuration"
    
    echo "Enter the name for the conda environment."
    echo "Press ENTER to use the default name: behav3d"
    echo ""
    read -p "Environment name [behav3d]: " USER_ENV_NAME
    
    if [ -n "$USER_ENV_NAME" ]; then
        ENV_NAME="$USER_ENV_NAME"
    fi
    
    ARGS="-n $ENV_NAME"
    echo ""
    print_info "Using environment name: $ENV_NAME"
else
    # Try to extract environment name from args for the final message
    if echo "$ARGS" | grep -qE "\-n|\-\-name"; then
        # Parse the -n or --name argument
        ENV_NAME=$(echo "$ARGS" | sed -n 's/.*\(-n\|--name\)[[:space:]]*\([^[:space:]]*\).*/\2/p')
        if [ -z "$ENV_NAME" ]; then
            ENV_NAME="behav3d"
        fi
    fi
fi

# =============================================================================
# Run the Python installer
# =============================================================================

print_header "Starting BEHAV3D Installation"

$PYTHON_CMD "$SCRIPT_DIR/install_behav3d.py" $ARGS
EXIT_CODE=$?

if [ $EXIT_CODE -ne 0 ]; then
    echo ""
    print_error "Installation failed with error code $EXIT_CODE"
    echo ""
    echo "If you see import errors, you may need to restart this script"
    echo "after Miniforge installation completes."
    echo ""
    echo "Press any key to close..."
    read -n 1
    exit $EXIT_CODE
fi

# =============================================================================
# Success message
# =============================================================================

echo ""
print_header "Installation Complete!"

echo "To use BEHAV3D:"
echo "  1. Open a new Terminal window"
echo "  2. Run: conda activate $ENV_NAME"
echo "  3. Navigate to the notebooks folder and run Jupyter"
echo ""
echo "Press any key to close..."
read -n 1
