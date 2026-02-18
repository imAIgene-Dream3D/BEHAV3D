@echo off
:: BEHAV3D Napari Launcher (Windows)
:: Double-click this file to open napari with the BEHAV3D plugin.
:: The correct environment is read from behav3d_env.json automatically.

title BEHAV3D - Launching Napari...
echo.
echo  ========================================
echo    BEHAV3D - Launching Napari
echo  ========================================
echo.

:: Run the Python launcher script (located in the same directory)
python "%~dp0launch_napari.py"

:: If python is not on PATH, try common locations
if errorlevel 1 (
    echo.
    echo Trying alternative Python paths...
    "%USERPROFILE%\miniforge3\python.exe" "%~dp0launch_napari.py" 2>nul
    if errorlevel 1 (
        "%USERPROFILE%\miniconda3\python.exe" "%~dp0launch_napari.py" 2>nul
        if errorlevel 1 (
            echo.
            echo ERROR: Could not find Python. Please ensure Python is installed and on your PATH.
            pause
        )
    )
)
