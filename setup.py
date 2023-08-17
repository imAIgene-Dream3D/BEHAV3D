import setuptools
from setuptools import setup

#https://betterscientificsoftware.github.io/python-for-hpc/tutorials/python-pypi-packaging/
# Use 'pip install -e .' in this folder to install all modules

setup(
    name='BEHAV3D',
    version='3.0',
    description='Module for the analysis of cell behaviour in fluorescent imaging',
    author='SdeBlank',
    packages=setuptools.find_packages()
)
