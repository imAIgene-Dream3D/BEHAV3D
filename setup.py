import setuptools
from setuptools import setup

# https://betterscientificsoftware.github.io/python-for-hpc/tutorials/python-pypi-packaging/
# Use 'pip install -e .' in this folder to install all modules
# https://python-packaging.readthedocs.io/en/latest/dependencies.html
setup(
    name='BEHAV3D',
    version='3.0',
    description='Module for the analysis of cell behaviour in fluorescent imaging',
    author='SdeBlank',
    author_email='S.deBlank-3@prinsesmaximacentrum.nl',
    license='Apache License 2.0',
    packages=setuptools.find_packages(),
    install_requires=[
          'pyyaml'
          'numpy',
          'pandas',
          'scyjava',
          'pyimagej',
          'tifffile',
          'scikit-image',
          'scipy',
          'matplotlib',
          'seaborn',
          'plotnine',
          'h5py'
      ]
)
