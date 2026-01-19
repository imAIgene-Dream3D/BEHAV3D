#!/usr/bin/env python3
"""
BEHAV3D Installer Script
========================

This script handles the complete installation of BEHAV3D, including:
1. Installing Miniconda/Miniforge if conda is not found
2. Detecting CUDA/GPU availability
3. Creating the conda environment with appropriate PyTorch version
4. Installing Cellpose with or without GPU support

Usage:
    python install_behav3d.py              # Full installation (default env: behav3d)
    python install_behav3d.py -n myenv     # Install with custom environment name
    python install_behav3d.py --cpu-only   # Force CPU-only installation
    python install_behav3d.py --pytorch-only  # Only install PyTorch/Cellpose (env exists)
    python install_behav3d.py --check      # Check system info only

Cross-platform: Windows, macOS, Linux
"""

import subprocess
import sys
import os
import platform
import shutil
import tempfile
import urllib.request
import argparse
from pathlib import Path

# =============================================================================
# CONFIGURATION
# =============================================================================

DEFAULT_ENV_NAME = "behav3d"
PYTHON_VERSION = "3.12"
CUDA_VERSION = "12.1"  # Default CUDA version for PyTorch

# Global variable that will be set from command line args
ENV_NAME = DEFAULT_ENV_NAME

# Miniforge download URLs (more reliable than Miniconda for conda-forge)
MINIFORGE_URLS = {
    "Windows": "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe",
    "Darwin": "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-{arch}.sh",
    "Linux": "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh",
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

class Colors:
    """ANSI color codes for terminal output"""
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

def print_header(text):
    print(f"\n{Colors.HEADER}{Colors.BOLD}{'='*60}{Colors.ENDC}")
    print(f"{Colors.HEADER}{Colors.BOLD}{text:^60}{Colors.ENDC}")
    print(f"{Colors.HEADER}{Colors.BOLD}{'='*60}{Colors.ENDC}\n")

def print_step(text):
    print(f"{Colors.CYAN}[STEP]{Colors.ENDC} {text}")

def print_success(text):
    print(f"{Colors.GREEN}[SUCCESS]{Colors.ENDC} {text}")

def print_warning(text):
    print(f"{Colors.WARNING}[WARNING]{Colors.ENDC} {text}")

def print_error(text):
    print(f"{Colors.FAIL}[ERROR]{Colors.ENDC} {text}")

def print_info(text):
    print(f"{Colors.BLUE}[INFO]{Colors.ENDC} {text}")

def get_platform():
    """Get the current platform."""
    return platform.system()

def get_architecture():
    """Get CPU architecture."""
    machine = platform.machine().lower()
    if machine in ('arm64', 'aarch64'):
        return 'arm64'
    return 'x86_64'

def run_command(cmd, shell=True, capture=False, check=True):
    """Run a shell command."""
    try:
        if capture:
            result = subprocess.run(cmd, shell=shell, capture_output=True, text=True, check=check)
            return result.stdout.strip()
        else:
            subprocess.run(cmd, shell=shell, check=check)
            return True
    except subprocess.CalledProcessError as e:
        if capture:
            return None
        raise e

# =============================================================================
# CONDA DETECTION & INSTALLATION
# =============================================================================

def find_conda():
    """Find conda executable."""
    # Check common locations
    conda_names = ['conda', 'mamba']
    
    for name in conda_names:
        path = shutil.which(name)
        if path:
            return path
    
    # Check common installation paths
    home = Path.home()
    common_paths = [
        home / "miniconda3" / "condabin" / "conda",
        home / "miniforge3" / "condabin" / "conda",
        home / "anaconda3" / "condabin" / "conda",
        home / "mambaforge" / "condabin" / "conda",
        Path("C:/Users") / os.environ.get("USERNAME", "") / "miniconda3" / "condabin" / "conda.bat",
        Path("C:/Users") / os.environ.get("USERNAME", "") / "miniforge3" / "condabin" / "conda.bat",
        Path("C:/ProgramData/miniconda3/condabin/conda.bat"),
        Path("C:/ProgramData/miniforge3/condabin/conda.bat"),
    ]
    
    for p in common_paths:
        if p.exists():
            return str(p)
    
    return None

def install_miniforge():
    """Download and install Miniforge."""
    print_step("Installing Miniforge (conda-forge distribution)...")
    
    system = get_platform()
    arch = get_architecture()
    
    if system == "Darwin":
        url = MINIFORGE_URLS[system].format(arch=arch)
    else:
        url = MINIFORGE_URLS[system]
    
    # Download installer
    print_info(f"Downloading from: {url}")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        if system == "Windows":
            installer_path = Path(tmpdir) / "miniforge_installer.exe"
        else:
            installer_path = Path(tmpdir) / "miniforge_installer.sh"
        
        urllib.request.urlretrieve(url, installer_path)
        print_success("Download complete")
        
        # Install
        install_path = Path.home() / "miniforge3"
        
        if system == "Windows":
            print_step("Running installer (this may take a few minutes)...")
            cmd = f'start /wait "" "{installer_path}" /InstallationType=JustMe /RegisterPython=0 /S /D={install_path}'
            run_command(cmd)
        else:
            print_step("Running installer...")
            os.chmod(installer_path, 0o755)
            cmd = f'bash "{installer_path}" -b -p "{install_path}"'
            run_command(cmd)
        
        # Initialize conda
        conda_path = install_path / "condabin" / ("conda.bat" if system == "Windows" else "conda")
        
        if conda_path.exists():
            print_success(f"Miniforge installed at: {install_path}")
            
            # Initialize shell
            print_step("Initializing conda for your shell...")
            if system == "Windows":
                run_command(f'"{conda_path}" init cmd.exe powershell', check=False)
            else:
                run_command(f'"{conda_path}" init bash zsh', check=False)
            
            return str(conda_path)
        else:
            print_error("Installation failed - conda not found after install")
            return None

# =============================================================================
# GPU/CUDA DETECTION
# =============================================================================

def check_nvidia_gpu():
    """Check if NVIDIA GPU is available."""
    system = get_platform()
    
    # Try nvidia-smi
    try:
        if system == "Windows":
            result = run_command("nvidia-smi", capture=True, check=False)
        else:
            result = run_command("nvidia-smi", capture=True, check=False)
        
        if result and "NVIDIA" in result:
            # Extract GPU name
            lines = result.split('\n')
            for line in lines:
                if "NVIDIA" in line and "Driver" not in line:
                    gpu_name = line.strip()
                    return True, gpu_name
            return True, "NVIDIA GPU detected"
    except:
        pass
    
    return False, None

def check_cuda_available():
    """Check if CUDA is available via PyTorch (if installed)."""
    try:
        result = run_command(
            f'python -c "import torch; print(torch.cuda.is_available())"',
            capture=True, check=False
        )
        return result == "True"
    except:
        return False

def check_apple_silicon():
    """Check if running on Apple Silicon (M1/M2/M3)."""
    if get_platform() == "Darwin" and get_architecture() == "arm64":
        return True
    return False

def get_gpu_info():
    """Get comprehensive GPU information."""
    info = {
        "platform": get_platform(),
        "architecture": get_architecture(),
        "nvidia_gpu": False,
        "nvidia_gpu_name": None,
        "apple_silicon": False,
        "cuda_available": False,
        "recommended": "cpu"
    }
    
    # Check NVIDIA
    has_nvidia, gpu_name = check_nvidia_gpu()
    info["nvidia_gpu"] = has_nvidia
    info["nvidia_gpu_name"] = gpu_name
    
    # Check Apple Silicon
    info["apple_silicon"] = check_apple_silicon()
    
    # Determine recommendation
    if has_nvidia:
        info["recommended"] = "cuda"
    elif info["apple_silicon"]:
        info["recommended"] = "mps"  # Metal Performance Shaders
    else:
        info["recommended"] = "cpu"
    
    return info

# =============================================================================
# ENVIRONMENT CREATION
# =============================================================================

def check_env_exists(conda_path, env_name):
    """Check if conda environment exists."""
    try:
        result = run_command(f'"{conda_path}" env list', capture=True)
        return env_name in result
    except:
        return False

def create_conda_environment(conda_path, env_file="environment.yml"):
    """Create the conda environment from yml file."""
    print_step(f"Creating conda environment '{ENV_NAME}'...")
    
    script_dir = Path(__file__).parent
    env_path = script_dir / env_file
    
    if not env_path.exists():
        print_error(f"Environment file not found: {env_path}")
        return False
    
    try:
        # Check if mamba is available (faster)
        mamba_path = shutil.which("mamba")
        installer = mamba_path if mamba_path else conda_path
        
        # Use -n flag to override the name in the YAML file
        cmd = f'"{installer}" env create -f "{env_path}" -n {ENV_NAME} -y'
        print_info(f"Running: {cmd}")
        run_command(cmd)
        print_success(f"Environment '{ENV_NAME}' created successfully")
        return True
    except subprocess.CalledProcessError as e:
        print_error(f"Failed to create environment: {e}")
        return False

def get_conda_run_prefix(conda_path, env_name):
    """Get the command prefix to run commands in conda environment."""
    system = get_platform()
    if system == "Windows":
        return f'"{conda_path}" run -n {env_name}'
    else:
        return f'"{conda_path}" run -n {env_name}'

# =============================================================================
# PYTORCH INSTALLATION
# =============================================================================

def install_pytorch(conda_path, gpu_info, force_cpu=False):
    """Install PyTorch with appropriate backend."""
    print_step("Installing PyTorch...")
    
    run_prefix = get_conda_run_prefix(conda_path, ENV_NAME)
    system = get_platform()
    
    if force_cpu:
        backend = "cpu"
        print_info("Installing CPU-only version (forced)")
    elif gpu_info["recommended"] == "cuda":
        backend = "cuda"
        print_info(f"Installing CUDA version (GPU: {gpu_info['nvidia_gpu_name']})")
    elif gpu_info["recommended"] == "mps":
        backend = "mps"
        print_info("Installing MPS version (Apple Silicon)")
    else:
        backend = "cpu"
        print_info("Installing CPU-only version (no GPU detected)")
    
    try:
        if backend == "cuda":
            # Install PyTorch with CUDA from pytorch channel
            cmd = f'{run_prefix} pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121'
        elif backend == "mps" or system == "Darwin":
            # macOS - use default PyTorch (includes MPS support)
            cmd = f'{run_prefix} pip install torch torchvision torchaudio'
        else:
            # CPU only
            cmd = f'{run_prefix} pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu'
        
        print_info(f"Running: {cmd}")
        run_command(cmd)
        print_success("PyTorch installed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print_error(f"Failed to install PyTorch: {e}")
        return False

# =============================================================================
# CELLPOSE INSTALLATION
# =============================================================================

def install_cellpose(conda_path, gpu_info, force_cpu=False):
    """Install Cellpose with appropriate backend."""
    print_step("Installing Cellpose...")
    
    run_prefix = get_conda_run_prefix(conda_path, ENV_NAME)
    
    try:
        # Cellpose automatically uses GPU if PyTorch has CUDA support
        # Pin to specific version for compatibility
        cmd = f'{run_prefix} pip install cellpose==3.1.1.2'
        print_info(f"Running: {cmd}")
        run_command(cmd)
        print_success("Cellpose installed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print_error(f"Failed to install Cellpose: {e}")
        return False

# =============================================================================
# BEHAV3D PACKAGE INSTALLATION
# =============================================================================

def install_behav3d_package(conda_path):
    """Install BEHAV3D package in development mode."""
    print_step("Installing BEHAV3D package...")
    
    run_prefix = get_conda_run_prefix(conda_path, ENV_NAME)
    script_dir = Path(__file__).parent
    
    try:
        cmd = f'{run_prefix} pip install -e "{script_dir}"'
        print_info(f"Running: {cmd}")
        run_command(cmd)
        print_success("BEHAV3D package installed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print_error(f"Failed to install BEHAV3D: {e}")
        return False

# =============================================================================
# VERIFICATION
# =============================================================================

def verify_installation(conda_path):
    """Verify the installation."""
    print_step("Verifying installation...")
    
    run_prefix = get_conda_run_prefix(conda_path, ENV_NAME)
    
    checks = [
        ("Python", "python --version"),
        ("NumPy", 'python -c "import numpy; print(numpy.__version__)"'),
        ("PyTorch", 'python -c "import torch; print(f\'PyTorch {torch.__version__}\')"'),
        ("CUDA Available", 'python -c "import torch; print(f\'CUDA: {torch.cuda.is_available()}\')"'),
        ("Cellpose", 'python -c "from importlib.metadata import version; print(f\'Cellpose {version(\\\"cellpose\\\")}\')"'),
        ("Napari", 'python -c "import napari; print(f\'Napari {napari.__version__}\')"'),
    ]
    
    all_passed = True
    for name, cmd in checks:
        try:
            result = run_command(f'{run_prefix} {cmd}', capture=True, check=False)
            if result:
                print_success(f"{name}: {result}")
            else:
                print_warning(f"{name}: Not available or failed")
                all_passed = False
        except Exception as e:
            print_warning(f"{name}: Check failed - {e}")
            all_passed = False
    
    return all_passed

# =============================================================================
# MAIN INSTALLATION FLOW
# =============================================================================

def print_system_info(gpu_info):
    """Print system information."""
    print_header("SYSTEM INFORMATION")
    print_info(f"Operating System: {get_platform()}")
    print_info(f"Architecture: {get_architecture()}")
    print_info(f"Python: {sys.version.split()[0]}")
    
    if gpu_info["nvidia_gpu"]:
        print_success(f"NVIDIA GPU: {gpu_info['nvidia_gpu_name']}")
    else:
        print_info("NVIDIA GPU: Not detected")
    
    if gpu_info["apple_silicon"]:
        print_success("Apple Silicon: Detected (MPS acceleration available)")
    
    print_info(f"Recommended backend: {gpu_info['recommended'].upper()}")

def main():
    global ENV_NAME
    
    parser = argparse.ArgumentParser(
        description="BEHAV3D Installation Script",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python install_behav3d.py              # Full installation (env: behav3d)
  python install_behav3d.py -n myenv     # Install with custom environment name
  python install_behav3d.py --cpu-only   # Force CPU-only installation
  python install_behav3d.py --pytorch-only  # Only install PyTorch/Cellpose
  python install_behav3d.py --check      # Check system info only
        """
    )
    parser.add_argument("-n", "--name", type=str, default=DEFAULT_ENV_NAME,
                       help=f"Name for the conda environment (default: {DEFAULT_ENV_NAME})")
    parser.add_argument("--cpu-only", action="store_true", 
                       help="Force CPU-only installation (no CUDA)")
    parser.add_argument("--pytorch-only", action="store_true",
                       help="Only install PyTorch and Cellpose (environment must exist)")
    parser.add_argument("--check", action="store_true",
                       help="Only check system information")
    parser.add_argument("--skip-behav3d", action="store_true",
                       help="Skip installing BEHAV3D package")
    parser.add_argument("--reinstall", action="store_true",
                       help="Remove and reinstall environment if it exists (no prompt)")
    parser.add_argument("--keep-existing", action="store_true",
                       help="Keep existing environment and only update PyTorch/Cellpose (no prompt)")
    
    args = parser.parse_args()
    
    # Set global environment name from args
    ENV_NAME = args.name
    
    print_header("BEHAV3D INSTALLER")
    print("Cross-platform installation for Windows, macOS, and Linux")
    print("="*60)
    
    # Get GPU info
    gpu_info = get_gpu_info()
    print_system_info(gpu_info)
    
    if args.check:
        return 0
    
    # Find or install conda
    print_header("CONDA SETUP")
    conda_path = find_conda()
    
    if conda_path:
        print_success(f"Conda found: {conda_path}")
    else:
        print_warning("Conda not found on system")
        response = input("\nWould you like to install Miniforge (recommended conda distribution)? [Y/n]: ")
        if response.lower() != 'n':
            conda_path = install_miniforge()
            if not conda_path:
                print_error("Failed to install Miniforge. Please install manually.")
                print_info("Visit: https://github.com/conda-forge/miniforge")
                return 1
            print_warning("Please restart your terminal and run this script again.")
            print_info("The conda installation needs a shell restart to take effect.")
            return 0
        else:
            print_error("Conda is required for BEHAV3D installation.")
            return 1
    
    # Check if environment exists
    env_exists = check_env_exists(conda_path, ENV_NAME)
    
    print_header("ENVIRONMENT SETUP")
    print_info(f"Environment name: {ENV_NAME}")
    
    if args.pytorch_only:
        if not env_exists:
            print_error(f"Environment '{ENV_NAME}' does not exist. Run full installation first.")
            return 1
        print_info(f"Environment '{ENV_NAME}' found. Installing PyTorch and Cellpose only.")
    else:
        if env_exists:
            print_warning(f"Environment '{ENV_NAME}' already exists.")
            
            # Handle based on command line flags or prompt user
            if args.reinstall:
                # --reinstall flag: remove and recreate without prompting
                print_step(f"Removing existing environment (--reinstall flag)...")
                run_command(f'"{conda_path}" env remove -n {ENV_NAME} -y')
                env_exists = False
            elif args.keep_existing:
                # --keep-existing flag: skip recreation without prompting
                print_info("Keeping existing environment (--keep-existing flag).")
                print_info("Will update PyTorch and Cellpose only.")
            else:
                # Interactive prompt
                print()
                print("Options:")
                print("  [1] Remove and reinstall environment (clean install)")
                print("  [2] Keep environment, only update PyTorch/Cellpose")
                print("  [3] Cancel installation")
                print()
                response = input("Choose an option [1/2/3]: ").strip()
                
                if response == '1':
                    print_step(f"Removing existing environment...")
                    run_command(f'"{conda_path}" env remove -n {ENV_NAME} -y')
                    env_exists = False
                elif response == '2':
                    print_info("Keeping existing environment. Will update PyTorch/Cellpose.")
                elif response == '3':
                    print_info("Installation cancelled by user.")
                    return 0
                else:
                    print_warning(f"Invalid option '{response}'. Defaulting to keep existing environment.")
                    print_info("Keeping existing environment. Will update PyTorch/Cellpose.")
        
        if not env_exists:
            if not create_conda_environment(conda_path):
                print_error("Failed to create conda environment.")
                return 1
    
    # Install PyTorch
    print_header("PYTORCH INSTALLATION")
    if not install_pytorch(conda_path, gpu_info, force_cpu=args.cpu_only):
        print_error("Failed to install PyTorch.")
        return 1
    
    # Install Cellpose
    print_header("CELLPOSE INSTALLATION")
    if not install_cellpose(conda_path, gpu_info, force_cpu=args.cpu_only):
        print_error("Failed to install Cellpose.")
        return 1
    
    # Install BEHAV3D package
    if not args.skip_behav3d and not args.pytorch_only:
        print_header("BEHAV3D PACKAGE INSTALLATION")
        if not install_behav3d_package(conda_path):
            print_warning("BEHAV3D package installation failed (non-critical)")
    
    # Verify installation
    print_header("VERIFICATION")
    verify_installation(conda_path)
    
    # Final message
    print_header("INSTALLATION COMPLETE")
    print_success("BEHAV3D has been installed successfully!")
    print()
    print("To activate the environment:")
    if get_platform() == "Windows":
        print(f"  conda activate {ENV_NAME}")
    else:
        print(f"  conda activate {ENV_NAME}")
    print()
    print("To start using BEHAV3D:")
    print("  1. Activate the environment")
    print("  2. Navigate to the notebooks folder")
    print("  3. Run: jupyter notebook")
    print()
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
