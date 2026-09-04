"""BEHAV3D Explorer - Google Colab bootstrap.

Runs the *real* napari GUI (``napari.Viewer`` + ``BEHAV3DWidget``) on a virtual
X display inside the Colab VM and streams it to the visitor's browser::

    napari (PyQt5) -> Xvfb -> x11vnc -> websockify -> noVNC -> browser tab

Nothing in ``behav3d/`` is modified or reimplemented; only ``$DISPLAY`` differs
from a normal desktop install. napari is started through the repository's own
launcher in payload mode (``napari/.config/launch_napari.py --internal``).

Typical use from the notebook::

    from colab_setup import bootstrap, open_viewer
    bootstrap()        # ~4 min on a cold runtime
    open_viewer()      # opens the GUI in a new browser tab

Every step is idempotent: re-running a cell after a disconnect skips whatever
is already in place.
"""

from __future__ import annotations

import os
import re
import shutil
import socket
import subprocess
import sys
import tarfile
import time
import zipfile
from pathlib import Path

# ==========================================================================
# CONFIGURATION - the two URLs below are the only things a maintainer must
# fill in. See demo/README.md for how to produce and host these two files.
# Both can also be overridden with environment variables, which is what the
# notebook exposes to visitors who want to point at their own data.
# ==========================================================================

#: conda-pack tarball of the Linux ``behav3d`` environment (see demo/build_env.sh).
ENV_URL = os.environ.get(
    "BEHAV3D_ENV_URL",
    "https://huggingface.co/datasets/HF_ORG_PLACEHOLDER/behav3d-demo-env/resolve/main/behav3d-env.tar.gz",
)

#: Demo dataset bundle: raw images + metadata.csv + PRECOMPUTED output/ folder.
DEMO_URL = os.environ.get(
    "BEHAV3D_DEMO_URL",
    "https://zenodo.org/records/ZENODO_RECORD_PLACEHOLDER/files/behav3d_demo.tar.gz?download=1",
)

REPO_URL = os.environ.get("BEHAV3D_REPO_URL", "https://github.com/imAIgene-Dream3D/BEHAV3D.git")
REPO_REF = os.environ.get("BEHAV3D_REPO_REF", "main")

ENV_PREFIX = Path(os.environ.get("BEHAV3D_ENV_PREFIX", "/opt/behav3d"))
REPO_DIR = Path(os.environ.get("BEHAV3D_REPO_DIR", "/content/BEHAV3D"))
DEMO_ROOT = Path(os.environ.get("BEHAV3D_DEMO_ROOT", "/content/behav3d_demo"))
LOG_DIR = Path(os.environ.get("BEHAV3D_LOG_DIR", "/content/behav3d_logs"))

DISPLAY = ":99"
SCREEN = os.environ.get("BEHAV3D_SCREEN", "1920x1080x24")
VNC_PORT = 5900
WEB_PORT = int(os.environ.get("BEHAV3D_WEB_PORT", "6080"))
WINDOW_TITLE = "BEHAV3D Explorer"

#: Everything napari/Qt needs to render without a GPU or a real display.
GUI_ENV = {
    "DISPLAY": DISPLAY,
    "LIBGL_ALWAYS_SOFTWARE": "1",   # llvmpipe software OpenGL
    "GALLIUM_DRIVER": "llvmpipe",
    "QT_X11_NO_MITSHM": "1",        # shared memory is restricted in containers
    "QT_QPA_PLATFORM": "xcb",
    "QT_API": "pyqt5",           # cellpose[gui] can pull PyQt6 in; pin qtpy to PyQt5
    "PYOPENCL_NO_CACHE": "1",       # matches napari/run_behav3d_linux.sh
    "PYOPENCL_COMPILER_OUTPUT": "0",
    "MPLBACKEND": "Agg",
    "PYTHONUNBUFFERED": "1",
}

_PROCS: dict = {}


# ==========================================================================
# small helpers
# ==========================================================================
def _say(msg):
    print(f"[behav3d] {msg}", flush=True)


def _run(cmd, check=True, **kw):
    """Run a command, raising with its output on failure."""
    res = subprocess.run(cmd, shell=isinstance(cmd, str), capture_output=True, text=True, **kw)
    if check and res.returncode != 0:
        raise RuntimeError(
            f"command failed ({res.returncode}): {cmd}\n"
            f"--- stdout ---\n{res.stdout[-2000:]}\n--- stderr ---\n{res.stderr[-2000:]}"
        )
    return res


def _spawn(name, cmd, env=None):
    """Start a long-lived background process, logging to LOG_DIR/<name>.log."""
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = open(LOG_DIR / f"{name}.log", "ab", buffering=0)
    proc = subprocess.Popen(
        cmd,
        stdout=log,
        stderr=subprocess.STDOUT,
        env={**os.environ, **GUI_ENV, **(env or {})},
        start_new_session=True,
    )
    _PROCS[name] = proc
    return proc


def _port_open(port, host="127.0.0.1"):
    with socket.socket() as s:
        s.settimeout(0.4)
        return s.connect_ex((host, port)) == 0


def _wait(predicate, timeout, what, interval=0.5):
    deadline = time.time() + timeout
    while time.time() < deadline:
        if predicate():
            return
        time.sleep(interval)
    raise TimeoutError(f"timed out after {timeout:.0f}s waiting for {what}")


def env_python():
    """Path to the python interpreter of the unpacked ``behav3d`` environment."""
    return ENV_PREFIX / "bin" / "python"


def _download(url, dest):
    """Fetch *url* to *dest* (resumable). A local path in *url* is just copied."""
    src = Path(url)
    if src.exists():
        if src.resolve() != Path(dest).resolve():
            shutil.copy2(src, dest)
        return Path(dest)
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    _say("downloading " + url.split("?")[0].rsplit("/", 1)[-1] + " ...")
    subprocess.run(
        ["wget", "--continue", "--no-verbose", "--show-progress",
         "--progress=dot:giga", "-O", str(dest), url],
        check=True,
    )
    return dest


def _extract(archive, dest):
    Path(dest).mkdir(parents=True, exist_ok=True)
    if Path(archive).suffix == ".zip":
        with zipfile.ZipFile(archive) as z:
            z.extractall(dest)
    else:
        # GNU tar is much faster than tarfile for multi-GB archives
        _run(["tar", "-xf", str(archive), "-C", str(dest)])


# ==========================================================================
# step 1 - system packages
# ==========================================================================
def install_system_packages(package_file=None):
    """apt-get the headless-Qt / OpenGL / VNC stack (Colab runs as root)."""
    if shutil.which("Xvfb") and shutil.which("x11vnc") and shutil.which("websockify"):
        _say("system packages already present - skipping apt")
        return

    path = Path(package_file) if package_file else Path(__file__).with_name("apt_packages.txt")
    if path.exists():
        packages = [
            line.split("#")[0].strip()
            for line in path.read_text().splitlines()
            if line.split("#")[0].strip()
        ]
    else:  # notebook may be run before the repo is checked out
        packages = [
            "xvfb", "x11vnc", "websockify", "novnc", "fluxbox", "x11-utils", "xdotool",
            "libgl1", "libglx-mesa0", "libgl1-mesa-dri", "libegl1", "libglib2.0-0",
            "libdbus-1-3", "libxkbcommon-x11-0", "libxcb-icccm4", "libxcb-image0",
            "libxcb-keysyms1", "libxcb-randr0", "libxcb-render-util0", "libxcb-shape0",
            "libxcb-shm0", "libxcb-sync1", "libxcb-xfixes0", "libxcb-xinerama0",
            "libxcb-xkb1", "libxcb-cursor0", "libsm6", "libice6", "libxi6",
            "libxrender1", "libfontconfig1", "libgomp1",
        ]

    _say(f"installing {len(packages)} system packages (~1 min) ...")
    env = {**os.environ, "DEBIAN_FRONTEND": "noninteractive"}
    _run("apt-get -qq update", env=env)
    _run(["apt-get", "-qq", "install", "-y", "--no-install-recommends"] + packages, env=env)
    _say("system packages ready")


# ==========================================================================
# step 2 - the conda environment (prebuilt tarball, never solved here)
# ==========================================================================
def install_environment(url=None, prefix=None):
    """Download + unpack the conda-pack'd ``behav3d`` environment.

    Solving ``installation/environment.yml`` inside Colab takes 10+ minutes and
    breaks whenever conda-forge moves; a packed tarball restores in ~2 minutes
    and is byte-identical to the environment the maintainer tested.
    """
    url = url or ENV_URL
    prefix = Path(prefix or ENV_PREFIX)
    if env_python().exists():
        _say(f"environment already unpacked at {prefix} - skipping")
        return prefix

    if "PLACEHOLDER" in url:
        raise RuntimeError(
            "ENV_URL still contains a placeholder. Build the tarball with "
            "demo/build_env.sh, upload it (see demo/README.md), then set "
            "BEHAV3D_ENV_URL or edit colab_setup.py."
        )

    archive = Path("/content/behav3d-env.tar.gz")
    _download(url, archive)
    _say(f"unpacking environment into {prefix} (~1-2 min) ...")
    prefix.mkdir(parents=True, exist_ok=True)
    _extract(archive, prefix)

    unpack = prefix / "bin" / "conda-unpack"
    if unpack.exists():
        _run([str(unpack)])          # rewrites the hard-coded build prefixes
    if not env_python().exists():
        raise RuntimeError(
            f"no python at {env_python()} - is the tarball really a conda-pack archive?"
        )

    archive.unlink(missing_ok=True)  # reclaim ~3 GB of Colab disk
    ver = _run([str(env_python()), "-c", "import napari; print(napari.__version__)"]).stdout.strip()
    _say(f"environment ready (napari {ver})")
    return prefix


# ==========================================================================
# step 3 - the BEHAV3D source + the demo dataset
# ==========================================================================
def clone_repo(url=None, ref=None, dest=None):
    """Shallow-clone BEHAV3D. It is used straight from source via PYTHONPATH.

    No ``pip install`` is needed: the launcher imports ``behav3d`` directly and
    adds the dock widget itself, so a path entry is enough (and is instant).
    """
    url, ref, dest = url or REPO_URL, ref or REPO_REF, Path(dest or REPO_DIR)
    if (dest / ".git").exists():
        _say(f"repository already at {dest} - skipping clone")
        return dest
    _say(f"cloning BEHAV3D ({ref}) ...")
    _run(["git", "clone", "--depth", "1", "--branch", ref, url, str(dest)])
    return dest


def fetch_demo_data(url=None, dest=None):
    """Download the demo bundle and repoint its absolute paths at this VM."""
    url, dest = url or DEMO_URL, Path(dest or DEMO_ROOT)
    if (dest / "metadata.csv").exists():
        _say(f"demo data already at {dest} - skipping")
        return dest

    if "PLACEHOLDER" in url:
        raise RuntimeError(
            "DEMO_URL still contains a placeholder. Upload the demo bundle "
            "(see demo/README.md), then set BEHAV3D_DEMO_URL or edit colab_setup.py."
        )

    archive = Path("/content/behav3d_demo_bundle.tar.gz")
    _download(url, archive)
    _say("extracting demo dataset ...")
    dest.mkdir(parents=True, exist_ok=True)
    # A bundle may or may not carry its own top-level folder; handle both.
    _extract(archive, dest.parent if _has_top_level_dir(archive, dest.name) else dest)
    archive.unlink(missing_ok=True)

    prepare = Path(__file__).with_name("prepare_demo.py")
    if not prepare.exists():
        prepare = REPO_DIR / "demo" / "colab" / "prepare_demo.py"
    _say("rewriting dataset paths for this machine ...")
    print(_run([str(env_python()), str(prepare), "--root", str(dest)]).stdout.strip())
    return dest


def _has_top_level_dir(archive, name):
    """True if the tarball already contains a ``<name>/`` top-level folder."""
    try:
        with tarfile.open(archive) as t:
            first = next((m.name for m in t), "")
        return first.split("/")[0] == name
    except Exception:
        return False


# ==========================================================================
# step 4 - virtual display + VNC + noVNC
# ==========================================================================
def start_display():
    """Xvfb -> fluxbox -> x11vnc -> websockify/noVNC, each waited on."""
    if _port_open(WEB_PORT):
        _say(f"display stack already running on port {WEB_PORT} - skipping")
        return

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    _say("starting virtual display ...")
    _spawn("xvfb", ["Xvfb", DISPLAY, "-screen", "0", SCREEN, "-ac",
                    "+extension", "GLX", "+extension", "RANDR", "+render", "-noreset"])
    _wait(lambda: _run(["xdpyinfo", "-display", DISPLAY], check=False).returncode == 0,
          30, "Xvfb")

    # A window manager is what gives napari a title bar, resizing and focus.
    _spawn("fluxbox", ["fluxbox"])
    time.sleep(1)

    _spawn("x11vnc", ["x11vnc", "-display", DISPLAY, "-forever", "-shared", "-nopw",
                      "-rfbport", str(VNC_PORT), "-noxdamage", "-quiet"])
    _wait(lambda: _port_open(VNC_PORT), 30, "x11vnc")

    _spawn("novnc", ["websockify", "--web=/usr/share/novnc",
                     str(WEB_PORT), f"localhost:{VNC_PORT}"])
    _wait(lambda: _port_open(WEB_PORT), 30, "noVNC")
    _say(f"display stack ready (noVNC on :{WEB_PORT})")


# ==========================================================================
# step 5 - napari itself
# ==========================================================================
def launch_napari(repo_dir=None, timeout=300):
    """Start BEHAV3D Explorer through the repo's own launcher (payload mode)."""
    repo_dir = Path(repo_dir or REPO_DIR)
    launcher = repo_dir / "napari" / ".config" / "launch_napari.py"
    if not launcher.exists():
        raise FileNotFoundError(f"launcher not found: {launcher}")

    if _window_id():
        _say("napari is already running - skipping")
        return _PROCS.get("napari")

    _say("launching napari (the first torch/napari import takes ~30-60 s) ...")
    proc = _spawn(
        "napari",
        [str(env_python()), str(launcher), "--internal"],
        env={"PYTHONPATH": str(repo_dir)},   # use BEHAV3D straight from source
    )

    def _ready():
        if proc.poll() is not None:
            raise RuntimeError("napari exited during startup:\n" + tail("napari", 40))
        return bool(_window_id())

    _wait(_ready, timeout, "the BEHAV3D Explorer window", interval=1.0)
    _maximise()
    _say("BEHAV3D Explorer is up")
    return proc


def _window_id():
    res = _run(["xdotool", "search", "--name", WINDOW_TITLE], check=False,
               env={**os.environ, "DISPLAY": DISPLAY})
    ids = [i for i in res.stdout.split() if i.strip()]
    return ids[-1] if ids else None


def _maximise():
    """Fill the virtual screen so the visitor sees the whole GUI immediately."""
    wid = _window_id()
    if not wid:
        return
    w, h = SCREEN.split("x")[:2]
    env = {**os.environ, "DISPLAY": DISPLAY}
    _run(["xdotool", "windowmove", wid, "0", "0"], check=False, env=env)
    _run(["xdotool", "windowsize", wid, w, h], check=False, env=env)
    _run(["xdotool", "windowactivate", wid], check=False, env=env)


# ==========================================================================
# step 6 - getting the pixels to the visitor
# ==========================================================================
_VNC_PATH = "/vnc.html?autoconnect=true&resize=remote&reconnect=true&quality=6"


def open_viewer(in_tab=True, height=800):
    """Open the GUI through Colab's own port proxy (new tab, or inline)."""
    try:
        from google.colab import output  # type: ignore
    except ImportError:
        _say(f"not running in Colab - open http://localhost:{WEB_PORT}{_VNC_PATH}")
        return None
    if not _port_open(WEB_PORT):
        raise RuntimeError("noVNC is not running - call start_display() first")
    if in_tab:
        _say("opening BEHAV3D Explorer in a new browser tab "
             "(allow pop-ups for colab.research.google.com if nothing happens)")
        return output.serve_kernel_port_as_window(WEB_PORT, path=_VNC_PATH)
    return output.serve_kernel_port_as_iframe(WEB_PORT, path=_VNC_PATH, height=str(height))


def start_cloudflared(timeout=60):
    """Fallback route: a free Cloudflare quick tunnel (no account needed).

    Colab's port proxy has had intermittent WebSocket failures, and noVNC needs
    WebSockets. This gives a plain ``https://<random>.trycloudflare.com`` URL
    that works from any browser, including a phone.
    """
    binary = Path("/usr/local/bin/cloudflared")
    if not binary.exists():
        _say("installing cloudflared ...")
        _download(
            "https://github.com/cloudflare/cloudflared/releases/latest/download/"
            "cloudflared-linux-amd64",
            binary,
        )
        binary.chmod(0o755)

    log = LOG_DIR / "cloudflared.log"
    log.unlink(missing_ok=True)
    _spawn("cloudflared", [str(binary), "tunnel", "--no-autoupdate",
                           "--url", f"http://localhost:{WEB_PORT}"])

    pattern = re.compile(r"https://[-\w]+\.trycloudflare\.com")
    deadline = time.time() + timeout
    while time.time() < deadline:
        if log.exists():
            match = pattern.search(log.read_text(errors="ignore"))
            if match:
                url = match.group(0) + _VNC_PATH
                _say("tunnel ready - open this URL:")
                print("\n    " + url + "\n")
                return url
        time.sleep(1)
    raise TimeoutError("cloudflared did not report a URL:\n" + tail("cloudflared", 20))


# ==========================================================================
# orchestration + troubleshooting
# ==========================================================================
def bootstrap(with_data=True):
    """Everything from a cold runtime to a running GUI. Safe to re-run."""
    t0 = time.time()
    steps = [("system packages", install_system_packages),
             ("environment", install_environment),
             ("repository", clone_repo)]
    if with_data:
        steps.append(("demo data", fetch_demo_data))
    steps += [("display", start_display), ("napari", launch_napari)]

    for label, fn in steps:
        t = time.time()
        fn()
        _say(f"  '{label}' done in {time.time() - t:.0f}s")
    _say(f"ready in {time.time() - t0:.0f}s total - now run open_viewer()")
    print_paths()


def print_paths():
    """The two paths a visitor types into the Data Preparation tab."""
    print(
        "\n" + "=" * 68
        + f"\n  Metadata CSV : {DEMO_ROOT / 'metadata.csv'}"
        + f"\n  Output folder: {DEMO_ROOT / 'output'}"
        + "\n" + "=" * 68 + "\n"
    )


def gpu_status():
    """Report whether this runtime got a GPU. BEHAV3D runs fine without one."""
    try:
        return _run([
            str(env_python()), "-c",
            "import torch;print('GPU:', torch.cuda.get_device_name(0)) "
            "if torch.cuda.is_available() else print('No GPU - running on CPU')",
        ]).stdout.strip()
    except Exception as exc:                      # environment not installed yet
        return f"unknown ({exc})"


def tail(name, n=30):
    log = LOG_DIR / f"{name}.log"
    if not log.exists():
        return f"(no log for '{name}')"
    return "\n".join(log.read_text(errors="ignore").splitlines()[-n:])


def status():
    """One-glance health check - the first thing to run when something breaks."""
    xvfb_up = _run(["xdpyinfo", "-display", DISPLAY], check=False).returncode == 0
    print(f"  environment   : {'ok' if env_python().exists() else 'MISSING'} ({ENV_PREFIX})")
    print(f"  repository    : {'ok' if (REPO_DIR / '.git').exists() else 'MISSING'} ({REPO_DIR})")
    print(f"  demo data     : {'ok' if (DEMO_ROOT / 'metadata.csv').exists() else 'MISSING'} ({DEMO_ROOT})")
    print(f"  Xvfb {DISPLAY}     : {'ok' if xvfb_up else 'down'}")
    print(f"  x11vnc :{VNC_PORT}  : {'ok' if _port_open(VNC_PORT) else 'down'}")
    print(f"  noVNC  :{WEB_PORT}  : {'ok' if _port_open(WEB_PORT) else 'down'}")
    print(f"  napari window : {'ok' if _window_id() else 'not found'}")
    for name, proc in _PROCS.items():
        state = "running" if proc.poll() is None else f"exited {proc.returncode}"
        print(f"  process {name:<12}: {state}")


def shutdown():
    """Stop napari and the display stack (bootstrap() brings it all back)."""
    for proc in _PROCS.values():
        if proc.poll() is None:
            proc.terminate()
    time.sleep(2)
    for proc in _PROCS.values():
        if proc.poll() is None:
            proc.kill()
    _PROCS.clear()
    _say("stopped")


if __name__ == "__main__":
    # Also usable outside Colab, e.g. inside the Docker dry-run container:
    #   python colab_setup.py --local
    if "--local" in sys.argv:
        start_display()
        launch_napari()
        _say(f"open http://localhost:{WEB_PORT}{_VNC_PATH}")
        _PROCS["napari"].wait()
    else:
        bootstrap()
