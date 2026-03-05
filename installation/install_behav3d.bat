@echo off
REM =============================================================================
REM BEHAV3D Installer - Windows Batch Script
REM =============================================================================
REM This script wraps the Python installer for easy double-click execution.
REM It will search for Python/Conda in common locations.
REM
REM Options (pass as arguments):
REM   -n NAME        Set custom environment name (default: behav3d)
REM   --cpu-only     Force CPU-only installation
REM   --reinstall    Remove and reinstall if environment exists
REM   --keep-existing  Keep existing environment, only update PyTorch/Cellpose
REM   --check        Only check system information
REM
REM Examples:
REM   install_behav3d.bat
REM   install_behav3d.bat -n my_behav3d_env
REM   install_behav3d.bat --cpu-only
REM   install_behav3d.bat -n test_env --reinstall
REM =============================================================================

setlocal EnableDelayedExpansion

echo.
echo ============================================================
echo              BEHAV3D Installation Script
echo ============================================================
echo.

REM Get the directory where this script is located
set "SCRIPT_DIR=%~dp0"

REM Initialize Python path variable
set "PYTHON_CMD="

REM First, try to find conda/miniforge and use its Python
echo [INFO] Searching for Python/Conda installation...
echo.

REM Check common Miniforge/Conda locations
set "CONDA_LOCATIONS=%USERPROFILE%\miniforge3 %USERPROFILE%\mambaforge %USERPROFILE%\miniconda3 %USERPROFILE%\anaconda3 C:\ProgramData\miniforge3 C:\ProgramData\miniconda3"

for %%L in (%CONDA_LOCATIONS%) do (
    if exist "%%L\python.exe" (
        set "PYTHON_CMD=%%L\python.exe"
        echo [FOUND] Conda Python at: %%L
        goto :found_python
    )
)

REM Check if Python is in PATH (but actually works - not the MS Store alias)
where python >nul 2>&1
if %ERRORLEVEL% EQU 0 (
    REM Test if it actually runs
    python --version >nul 2>&1
    if %ERRORLEVEL% EQU 0 (
        set "PYTHON_CMD=python"
        echo [FOUND] Python in PATH
        goto :found_python
    )
)

REM Check common Python installation paths
set "PYTHON_LOCATIONS=%LOCALAPPDATA%\Programs\Python\Python312 %LOCALAPPDATA%\Programs\Python\Python311 %LOCALAPPDATA%\Programs\Python\Python310 C:\Python312 C:\Python311 C:\Python310"

for %%L in (%PYTHON_LOCATIONS%) do (
    if exist "%%L\python.exe" (
        set "PYTHON_CMD=%%L\python.exe"
        echo [FOUND] Python at: %%L
        goto :found_python
    )
)

REM Python not found - offer to download Miniforge
echo.
echo ============================================================
echo   Python/Conda not found on this system
echo ============================================================
echo.
echo BEHAV3D requires Python. The easiest solution is to install
echo Miniforge, which includes Python and conda package manager.
echo.
echo Would you like to download and install Miniforge now?
echo.
set /p INSTALL_MINIFORGE="Enter Y to install, N to cancel [Y/n]: "

if /i "%INSTALL_MINIFORGE%"=="n" (
    echo.
    echo Installation cancelled. Please install Python or Miniforge manually:
    echo   https://github.com/conda-forge/miniforge/releases/latest/
    echo.
    pause
    exit /b 1
)

REM Download and install Miniforge
echo.
echo [INFO] Downloading Miniforge installer...
echo.

set "MINIFORGE_URL=https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe"
set "MINIFORGE_INSTALLER=%TEMP%\Miniforge3-installer.exe"

REM Use PowerShell to download
powershell -Command "& {[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; Invoke-WebRequest -Uri '%MINIFORGE_URL%' -OutFile '%MINIFORGE_INSTALLER%'}"

if not exist "%MINIFORGE_INSTALLER%" (
    echo [ERROR] Failed to download Miniforge installer.
    echo Please download manually from:
    echo   https://github.com/conda-forge/miniforge/releases/latest/
    pause
    exit /b 1
)

echo [INFO] Installing Miniforge (this may take a few minutes)...
echo.

REM Run silent install
start /wait "" "%MINIFORGE_INSTALLER%" /InstallationType=JustMe /RegisterPython=0 /AddToPath=1 /S /D=%USERPROFILE%\miniforge3

REM Check if installation succeeded
if exist "%USERPROFILE%\miniforge3\python.exe" (
    set "PYTHON_CMD=%USERPROFILE%\miniforge3\python.exe"
    echo.
    echo [SUCCESS] Miniforge installed successfully!
    echo.
    
    REM Clean up installer
    del "%MINIFORGE_INSTALLER%" >nul 2>&1
    
    goto :found_python
) else (
    echo [ERROR] Miniforge installation may have failed.
    echo Please try installing manually from:
    echo   https://github.com/conda-forge/miniforge/releases/latest/
    pause
    exit /b 1
)

:found_python
echo.
echo [INFO] Using Python: %PYTHON_CMD%
echo.

REM Verify Python works
"%PYTHON_CMD%" --version
if %ERRORLEVEL% NEQ 0 (
    echo [ERROR] Python found but failed to execute.
    pause
    exit /b 1
)

REM Check if any arguments were provided
set "ARGS=%*"
set "ENV_NAME=behav3d"

REM If no arguments, prompt for environment name
if "%ARGS%"=="" (
    echo.
    echo ============================================================
    echo              Environment Configuration
    echo ============================================================
    echo.
    echo Enter the name for the conda environment.
    echo Press ENTER to use the default name: behav3d
    echo.
    set /p "USER_ENV_NAME=Environment name [behav3d]: "
    
    if not "!USER_ENV_NAME!"=="" (
        set "ENV_NAME=!USER_ENV_NAME!"
    )
    
    set "ARGS=-n !ENV_NAME!"
    echo.
    echo [INFO] Using environment name: !ENV_NAME!
) else (
    REM Try to extract environment name from args for the final message
    echo !ARGS! | findstr /i "\-n \-\-name" >nul
    if !ERRORLEVEL! EQU 0 (
        REM Parse -n or --name argument
        for %%a in (%*) do (
            if "!NEXT_IS_NAME!"=="1" (
                set "ENV_NAME=%%a"
                set "NEXT_IS_NAME=0"
            )
            if "%%a"=="-n" set "NEXT_IS_NAME=1"
            if "%%a"=="--name" set "NEXT_IS_NAME=1"
        )
    )
)

echo.
echo ============================================================
echo              Starting BEHAV3D Installation
echo ============================================================
echo.

REM Run the Python installer
"%PYTHON_CMD%" "%SCRIPT_DIR%install_behav3d.py" %ARGS%

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo [ERROR] Installation failed with error code %ERRORLEVEL%
    echo.
    echo If you see import errors, you may need to restart this script
    echo after Miniforge installation completes.
    pause
    exit /b %ERRORLEVEL%
)

echo.
echo Press any key to close this window...
pause >nul
