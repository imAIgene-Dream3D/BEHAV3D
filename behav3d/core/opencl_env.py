"""Process-wide PyOpenCL setup shared by every GPU path (APOC / pyclesperanto).

Call :func:`configure_pyopencl` *before* ``pyopencl``, ``apoc`` or
``pyclesperanto_prototype`` are imported, so the environment variables are in
place by the time PyOpenCL reads them.
"""

import os
import warnings

# Start of the message pyopencl.compiler_output() emits when
# PYOPENCL_COMPILER_OUTPUT is off (pyopencl/__init__.py).
_COMPILER_OUTPUT_MESSAGE = "Non-empty compiler output encountered"


def silence_pyopencl_compiler_warnings():
    """Hide PyOpenCL's cosmetic "Non-empty compiler output encountered" warning.

    ``PYOPENCL_COMPILER_OUTPUT`` only selects *which* text PyOpenCL emits, not
    whether it warns: ``compiler_output()`` calls ``warnings.warn()`` either way
    whenever an OpenCL driver prints anything while building a kernel (vendor
    banners, register-spill notes, ...), which most drivers do on every build.
    napari routes warnings through its notification manager, so training a pixel
    classifier filled the activity area with this one message. Real build
    failures raise ``RuntimeError`` and are unaffected.

    Matching on the message instead of ``pyopencl.CompilerWarning`` keeps this
    module importable without pulling in PyOpenCL, and leaves the verbose
    variant (``PYOPENCL_COMPILER_OUTPUT=1``) visible when debugging a build.

    Safe to call repeatedly: ``filterwarnings`` de-duplicates identical entries.
    """
    warnings.filterwarnings(
        "ignore",
        message=_COMPILER_OUTPUT_MESSAGE,
        # pyopencl.CompilerWarning subclasses UserWarning, and filters match
        # subclasses, so this needs no import of pyopencl.
        category=UserWarning,
    )


def configure_pyopencl():
    """Set the PyOpenCL environment BEHAV3D needs and silence build chatter."""
    # Defeat buggy PyOpenCL compiler caching that causes TypeErrors on some
    # systems (e.g. enqueue_knl_predict() argument-count mismatches from a
    # stale cached kernel).
    os.environ['PYOPENCL_NO_CACHE'] = '1'
    # Only surface compiler output when something actually failed.
    os.environ['PYOPENCL_COMPILER_OUTPUT'] = '0'
    silence_pyopencl_compiler_warnings()
