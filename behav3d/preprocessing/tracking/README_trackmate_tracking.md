## Instruction to set up BEHAV3D data processing


```
conda create -n BEHAV3D python=3.11
conda install numpy
conda install pandas
conda install pyyaml
conda install h5py
conda install tifffile
conda install scikit-image
conda install -c conda-forge pyimagej
conda install java
```

### Segmentation (Ilastik)
1. Install Ilastik (https://www.ilastik.org/download.html)
2. Install Fiji (https://imagej.net/software/fiji/downloads)
2. Set *ilastik_path* in your 'config.yml' file
```
For Mac:
.../ilastik-1.3.2-OSX.app/Contents/ilastik-release/run_ilastik.sh

For Linux:
.../ilastik-1.3.2-Linux/run_ilastik.sh

For Windows (untested):
.../Program Files/ilastik-1.3.2/ilastik.exe
```

<!-- 
https://forum.image.sc/t/running-the-fiji-trackmate-plugin-from-python/68066
https://forum.image.sc/t/running-trackmate-using-pyimagej-headless-on-mac/69129/2
 -->
