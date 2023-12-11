import setuptools
from setuptools import setup

setup(
    name='BEHAV3D',
    version='3.0',
    description='Module for the analysis of cell behaviour in fluorescent imaging',
    author='SdeBlank',
    author_email='S.deBlank-3@prinsesmaximacentrum.nl',
    license='Apache License 2.0',
    packages=setuptools.find_packages(),
    install_requires=[
          'pyyaml==6.0.1',
          'numpy==1.26.0',
          'pandas==2.1.3',
          'scyjava==1.9.1',
          'pyimagej==1.4.1',
          'tifffile==2023.9.26',
          'scikit-image==0.22.0',
          'scipy==1.11.3',
          'seaborn==0.13.0',
          'plotnine==0.12.4',
          'patchworklib==0.6.3',
          'h5py==3.10.0',
          'ipykernel==6.26.0',
          'napari==0.4.18',
          'pyqt5',
          "umap-learn==0.5.5",
          'scikit-learn==1.3.2'
      ]
)
