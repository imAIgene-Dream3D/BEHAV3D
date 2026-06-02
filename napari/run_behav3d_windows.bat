@echo off
:: BEHAV3D Napari Launcher (Windows)
:: Double-click this file to open napari with the BEHAV3D plugin.
:: The correct environment is read from .conf\behav3d_env.json automatically.

title BEHAV3D - Launching Napari...
echo.
echo  ========================================
echo    BEHAV3D - Launching Napari
echo  ========================================
echo.

:: Suppress PyOpenCL compiler cache warnings
set PYOPENCL_NO_CACHE=1
:: set PYOPENCL_COMPILER_OUTPUT=0

:: Run the Python launcher script (located in the .conf subfolder)
python "%~dp0.config\launch_napari.py"

:: If python is not on PATH, try common locations
if errorlevel 1 (
    echo.
    echo Trying alternative Python paths...
    "%USERPROFILE%\miniforge3\python.exe" "%~dp0.config\launch_napari.py"
    if errorlevel 1 (
        "%USERPROFILE%\miniconda3\python.exe" "%~dp0.config\launch_napari.py"
        if errorlevel 1 (
            echo.
            echo ERROR: Could not find Python. Please ensure Python is installed and on your PATH.
            pause
        )
    )
)

:: Keep window open if there was an error or if run via double-click
if errorlevel 1 pause
