#!/usr/bin/env python3
"""
BEHAV3D Napari Launcher
=======================
Single script that handles both:
1. Orchestrating the environment launch (reading config, running conda run)
2. Actually starting napari with the BEHAV3D plugin (when run inside the env)

Usage:
    python launch_napari.py             -> Launcher mode (runs outside env)
    python launch_napari.py --internal  -> Payload mode (runs inside env)
"""
import sys
import os
import argparse
import platform
from pathlib import Path

# =============================================================================
# PAYLOAD MODE (Runs inside the environment)
# =============================================================================
def run_napari_payload():
    """Import napari, create viewer, add BEHAV3D widget, and start event loop."""
    print("Starting BEHAV3D Napari Plugin...")
    # Suppress PyOpenCL compiler cache warnings
    os.environ.setdefault('PYOPENCL_NO_CACHE', '1')
    # os.environ.setdefault('PYOPENCL_COMPILER_OUTPUT', '0')
    try:
        import napari
        from behav3d.napari._widget import BEHAV3DWidget
    except ImportError as e:
        print(f"ERROR: Could not import required modules: {e}")
        print("Ensure you are running in the 'behav3d' environment.")
        sys.exit(1)

    # Create the viewer
    viewer = napari.Viewer(title="BEHAV3D")
    
    # Add our dock widget
    widget = BEHAV3DWidget(viewer)
    viewer.window.add_dock_widget(widget, name="BEHAV3D Pipeline", area="right")
    
    # Start the event loop
    napari.run()


# =============================================================================
# LAUNCHER MODE (Runs outside, orchestrates subprocess)
# =============================================================================

def _find_env_python(pkg_manager, env_name):
    """Find the Python binary inside a conda/mamba/micromamba environment.
    
    For micromamba we avoid 'micromamba run' because it generates a wrapper
    script that uses 'exec --', which fails on Ubuntu where /bin/sh is dash.
    Instead, we locate the environment's Python binary directly.
    """
    pkg_path = Path(pkg_manager)
    pkg_stem = pkg_path.stem.lower()  # e.g. 'micromamba', 'conda', 'mamba'
    system = platform.system()
    
    # --- Strategy 1: resolve env prefix via '<pkg> env list --json' ---
    try:
        import subprocess, json
        result = subprocess.run(
            [str(pkg_path), "env", "list", "--json"],
            capture_output=True, text=True, check=False
        )
        if result.returncode == 0:
            data = json.loads(result.stdout)
            envs = data.get("envs", [])
            for env_path_str in envs:
                env_path = Path(env_path_str)
                if env_path.name == env_name:
                    if system == "Windows":
                        py = env_path / "python.exe"
                    else:
                        py = env_path / "bin" / "python"
                    if py.exists():
                        return str(py)
    except Exception:
        pass
    
    # --- Strategy 2: common env locations ---
    home = Path.home()
    candidate_roots = []
    
    if "micromamba" in pkg_stem:
        # micromamba default root prefix
        candidate_roots.append(home / "micromamba" / "envs")
        # MAMBA_ROOT_PREFIX env var
        mamba_root = os.environ.get("MAMBA_ROOT_PREFIX")
        if mamba_root:
            candidate_roots.append(Path(mamba_root))
    
    # miniforge / mambaforge / miniconda / anaconda
    for base in ["miniforge3", "mambaforge", "miniconda3", "anaconda3"]:
        candidate_roots.append(home / base / "envs")
    
    for root in candidate_roots:
        env_dir = root / env_name
        if system == "Windows":
            py = env_dir / "python.exe"
        else:
            py = env_dir / "bin" / "python"
        if py.exists():
            return str(py)
    
    return None


def run_launcher():
    """Read config and spawn subprocess in the correct environment."""
    import json
    import subprocess
    
    CONFIG_NAME = "behav3d_env.json"
    script_dir = Path(__file__).parent.resolve()
    config_path = script_dir / CONFIG_NAME

    if not config_path.exists():
        print(f"ERROR: {CONFIG_NAME} not found in {script_dir}")
        print("Please run the BEHAV3D installer first (install_behav3d.py).")
        input("Press Enter to close...")
        sys.exit(1)

    with open(config_path) as f:
        cfg = json.load(f)

    pkg_manager = cfg.get("pkg_manager", "")
    env_name = cfg.get("env_name", "behav3d")

    if not pkg_manager or not Path(pkg_manager).exists():
        print(f"ERROR: Package manager not found at: {pkg_manager}")
        print("The BEHAV3D environment may have been moved or deleted.")
        print("Please re-run the installer (install_behav3d.py).")
        input("Press Enter to close...")
        sys.exit(1)

    script_path = script_dir / "launch_napari.py"
    
    # ── Decide launch strategy ──
    pkg_stem = Path(pkg_manager).stem.lower()
    
    # For micromamba: bypass 'micromamba run' (broken on Ubuntu/dash)
    # and run the environment's Python directly.
    # For conda/mamba: try direct Python first, fall back to 'conda run'.
    env_python = _find_env_python(pkg_manager, env_name)
    
    if env_python:
        # Direct invocation — most reliable across all platforms
        cmd = [env_python, str(script_path), "--internal"]
        cmd_display = f'"{env_python}" "{script_path}" --internal'
        use_shell = False
    else:
        # Fall back to conda/mamba run (won't work for micromamba on Ubuntu)
        cmd = f'"{pkg_manager}" run --no-capture-output -n {env_name} python "{script_path}" --internal'
        cmd_display = cmd
        use_shell = True
    
    print(f"Launching napari in '{env_name}' environment...")
    print(f"  Command: {cmd_display}")
    print()

    try:
        subprocess.run(cmd, shell=use_shell, check=True)
    except subprocess.CalledProcessError as e:
        print(f"\nERROR: napari exited with error code {e.returncode}")
        input("Press Enter to close...")
        sys.exit(e.returncode)
    except KeyboardInterrupt:
        pass


def main():
    parser = argparse.ArgumentParser(description="BEHAV3D Napari Launcher")
    parser.add_argument("--internal", action="store_true", help="Internal flag for payload mode")
    args = parser.parse_args()

    if args.internal:
        run_napari_payload()
    else:
        run_launcher()

if __name__ == "__main__":
    main()
