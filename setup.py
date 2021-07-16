#!/usr/bin/env python
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Nov 2020                           #
############################################################


# Always prefer setuptools over distutils
import os
from setuptools import setup, find_packages


# Grab from README file: long_description
with open("docs/README.md", "r") as f:
    long_description = f.read()

# Grab from version.py file: version and description
with open("mintpy/version.py", "r") as f:
    lines = f.readlines()
    # version
    line = [line for line in lines if line.startswith("def get_release_info")][0].strip()
    version = line.replace("'",'"').split('"')[1].split('v')[1]
    # description
    line = [line for line in lines if line.startswith("description")][0].strip()
    description = line.replace("'",'"').split('"')[1]
    # website
    line = [line for line in lines if line.startswith("website")][0].strip()
    website = line.replace("'",'"').split('"')[1]


def do_setup():
    setup(
        name="mintpy",
        version=version,
        description=description,
        long_description=long_description,
        long_description_content_type="text/markdown",
        url=website,
        author="Zhang Yunjun, Heresh Fattahi",
        author_email="yunjunzgeo@gmail.com",

        classifiers=[
            "Development Status :: 4 - Beta",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            "Operating System :: OS Independent",
            "Programming Language :: Python",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
        ],
        download_url=("https://github.com/insarlab/MintPy/archive/v{}.tar.gz".format(version)),
        keywords="InSAR, deformation, time-series, volcano, earthquake, tectonics, geodesy, geophysics, remote-sensing",

        # package discovery
        packages=find_packages(),
        scripts=['mintpy/smallbaselineApp.py', 'mintpy/view.py', 'mintpy/tsview.py'],

        # dependencies
        python_requires=">=3.6",
        install_requires=[
            "cartopy",
            "cvxopt",
            "dask>=1.0",
            "dask-jobqueue>=0.3",
            "defusedxml",
            "h5py",
            "joblib",
            "lxml",
            "matplotlib",
            "numpy",
            "pyproj",
            "setuptools",
            "scikit-image",
            "scipy",
            # pyaps dependencies
            "cdsapi",
            "ecCodes",
            "netcdf4",
            "pygrib",
            "pyhdf",
            # pyresample dependencies
            "pyresample",
            #"openmp",
            "pykdtree",
            "xarray",
            "zarr",
        ],
        dependency_links=[
            "git+https://github.com/insarlab/PySolid.git",
            "git+https://github.com/tylere/pykml.git",
        ],

        # data files
        include_package_data=True,
        package_data={
            "mintpy": [
                "data/*.js",
                "data/*.png",
                "data/colormaps/*.cpt",
                "data/input_files/*.txt",
                "data/input_files/*.template",
                "data/input_files/*.md",
                "defaults/*.cfg",
                "defaults/*.yaml",
                "sh/*.sh",
            ],
        },

        project_urls={
            "Bug Reports": "https://github.com/insarlab/mintpy/issues",
            "Documentation": "https://mintpy.readthedocs.io/",
            "Source": "https://github.com/insarlab/mintpy",
        },
    )


if __name__ == "__main__":
    do_setup()
