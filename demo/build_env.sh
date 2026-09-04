#!/usr/bin/env bash
# =============================================================================
# Build the prebuilt Linux environment used by the Colab demo.
#
# Solving installation/environment.yml inside Colab takes 10+ minutes and breaks
# whenever conda-forge moves. Instead we solve it ONCE here, inside an Ubuntu
# 22.04 container (same base as Colab), and ship the result as a conda-pack
# tarball that Colab restores in ~2 minutes.
#
#   ./demo/build_env.sh                 # CPU build  (~2.5 GB packed) - default
#   ./demo/build_env.sh --cuda          # CUDA build (~5 GB packed), for GPU runtimes
#   ./demo/build_env.sh --test          # after building: run the GUI locally on :6080
#
# Output: dist/behav3d-env.tar.gz  (upload it - see demo/README.md)
# Requires: Docker. Takes 20-40 min the first time, minutes after that.
# =============================================================================
set -euo pipefail

FLAVOUR="cpu"
OUTPUT_DIR="dist"
RUN_TEST=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cuda)    FLAVOUR="cuda" ;;
    --cpu)     FLAVOUR="cpu" ;;
    --test)    RUN_TEST=1 ;;
    --output)  OUTPUT_DIR="$2"; shift ;;
    -h|--help) sed -n '2,16p' "$0"; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
  shift
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMAGE="behav3d-demo-env:${FLAVOUR}"
cd "$REPO_ROOT"

# CUDA 12.8 matches installation/install_behav3d.py (CUDA_VERSION = "cu128").
if [[ "$FLAVOUR" == "cuda" ]]; then
  TORCH_INDEX="--index-url https://download.pytorch.org/whl/cu128"
else
  TORCH_INDEX="--index-url https://download.pytorch.org/whl/cpu"
fi

# The apt list is shared with the Colab bootstrap so the two never drift.
APT_PACKAGES="$(sed 's/#.*//' demo/colab/apt_packages.txt | tr '\n' ' ' | tr -s ' ')"

echo "==> building ${IMAGE} (torch: ${FLAVOUR})"
DOCKERFILE="$(mktemp)"
trap 'rm -f "$DOCKERFILE"' EXIT

cat > "$DOCKERFILE" <<DOCKERFILE_EOF
# Ubuntu 22.04 base = same libc/mesa generation as Google Colab, so the packed
# environment is guaranteed to run there.
FROM mambaorg/micromamba:1.5.8-jammy
USER root
ENV DEBIAN_FRONTEND=noninteractive MAMBA_ROOT_PREFIX=/opt/conda

RUN apt-get update && apt-get install -y --no-install-recommends \\
        git wget ca-certificates ${APT_PACKAGES} \\
    && rm -rf /var/lib/apt/lists/*

# --- the environment, solved exactly once ---------------------------------
COPY environment.yml /tmp/environment.yml
RUN micromamba create -y -n behav3d -f /tmp/environment.yml && micromamba clean -a -y

# --- PyTorch + Cellpose, matching installation/install_behav3d.py ----------
RUN micromamba run -n behav3d pip install --no-cache-dir \\
        torch torchvision torchaudio ${TORCH_INDEX} \\
 && micromamba run -n behav3d pip install --no-cache-dir "cellpose[gui]==3.1.1.2" \\
 && micromamba run -n behav3d pip install --no-cache-dir conda-pack

# --- fail the build here rather than in front of a visitor -----------------
RUN micromamba run -n behav3d python -c "\\
import napari, torch, cellpose, pyopencl, qtpy; \\
print('napari', napari.__version__, '| torch', torch.__version__, '| qt', qtpy.API_NAME)"

# --- pack it ---------------------------------------------------------------
RUN micromamba run -n behav3d conda-pack -p /opt/conda/envs/behav3d \\
        -o /behav3d-env.tar.gz --ignore-missing-files --n-threads -1 \\
 && ls -lh /behav3d-env.tar.gz
DOCKERFILE_EOF

# Build context is installation/ only: the repo root carries ~1 GB of images and
# history that Docker would otherwise upload to the daemon on every build.
docker build -f "$DOCKERFILE" -t "$IMAGE" installation/

echo "==> extracting the tarball"
mkdir -p "$OUTPUT_DIR"
CONTAINER="$(docker create "$IMAGE")"
docker cp "${CONTAINER}:/behav3d-env.tar.gz" "${OUTPUT_DIR}/behav3d-env.tar.gz"
docker rm -f "$CONTAINER" >/dev/null

SIZE="$(du -h "${OUTPUT_DIR}/behav3d-env.tar.gz" | cut -f1)"
echo "==> ${OUTPUT_DIR}/behav3d-env.tar.gz (${SIZE})"
echo "    Next: upload it (demo/README.md, step 2) and put the URL in"
echo "    demo/colab/colab_setup.py -> ENV_URL"

if [[ "$RUN_TEST" == "1" ]]; then
  cat <<'TEST_EOF'

==> local dry run: starting the display stack + napari inside the image.
    When it says the window is up, open  http://localhost:6080/vnc.html
    (Ctrl-C here to stop.)

TEST_EOF
  docker run --rm -it -p 6080:6080 --shm-size=1g \
    -e BEHAV3D_ENV_PREFIX=/opt/conda/envs/behav3d \
    -e BEHAV3D_REPO_DIR=/work \
    -e BEHAV3D_LOG_DIR=/tmp/behav3d_logs \
    -v "${REPO_ROOT}:/work" \
    "$IMAGE" \
    /opt/conda/envs/behav3d/bin/python /work/demo/colab/colab_setup.py --local
fi
