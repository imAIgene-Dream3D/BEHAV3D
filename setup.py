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
    package_data={
        "behav3d.napari": ["napari.yaml"],
        "behav3d": [
            "resources/*.png",
            "resources/tutorials/*.png",
            "resources/tracking/*.png",
        ],
    },
    install_requires=[
        "napari>=0.5.0",
        "magicgui>=0.8",
        "qtpy>=2.4",
        "dask",
        "zarr",
        "pandas",
        "pyyaml",
    ],
    entry_points={
        "napari.manifest": [
            "BEHAV3D = behav3d.napari:napari.yaml",
        ],
    },
)
