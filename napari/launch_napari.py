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
import argparse
from pathlib import Path

# =============================================================================
# PAYLOAD MODE (Runs inside the environment)
# =============================================================================
def run_napari_payload():
    """Import napari, create viewer, add BEHAV3D widget, and start event loop."""
    print("Starting BEHAV3D Napari Plugin...")
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

    # Build the command: "<pkg_manager>" run -n <env_name> python <this_script> --internal
    # We use __file__ to refer to this same script
    script_path = script_dir / "launch_napari.py"
    
    cmd = f'"{pkg_manager}" run -n {env_name} python "{script_path}" --internal'
    print(f"Launching napari in '{env_name}' environment...")
    print(f"  Command: {cmd}")
    print()

    # On Windows, we need to handle signal differently often, but basic run is fine
    try:
        subprocess.run(cmd, shell=True, check=True)
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
