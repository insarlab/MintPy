#!/usr/bin/env python
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Nov 2020                           #
############################################################

import os
import sys

# Always prefer setuptools over distutils
from setuptools import find_packages, setup

# Grab version and description from version.py
# link: https://stackoverflow.com/questions/53648900
sys.path.append(os.path.join(os.path.dirname(__file__), "src"))
from mintpy.version import description, version

# Grab long_description from README.md
with open("docs/README.md") as f:
    long_description = f.read()

setup(
    name="mintpy",
    version=version,
    description=description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/insarlab/MintPy",
    download_url=(f"https://github.com/insarlab/MintPy/archive/v{version}.tar.gz"),
    author="Zhang Yunjun, Heresh Fattahi",
    author_email="yunjunz@outlook.com",
    license="GPL-3.0-or-later",
    license_files=("LICENSE",),

    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
    ],
    keywords="InSAR, deformation, time-series, volcano, earthquake, tectonics, geodesy, geophysics, remote-sensing",

    project_urls={
        "Bug Reports": "https://github.com/insarlab/mintpy/issues",
        "Documentation": "https://mintpy.readthedocs.io/",
        "Source": "https://github.com/insarlab/mintpy",
    },

    # dependencies
    python_requires=">=3.6",
    install_requires=[
        "argcomplete",
        "cartopy",
        "cvxopt",
        "dask>=1.0",
        "dask-jobqueue>=0.3",
        "h5py",
        "joblib",
        "lxml",
        "matplotlib",
        "numpy",
        "pre-commit",  # for developers
        "pyaps3>=0.3",
        "pykml>=0.2",
        "pyproj",
        "pyresample",  # pip installed version does not work
        "pysolid",     # pip installed version does not work because Fortran compiler is needed but not available via pip
        "rich",
        "setuptools",
        "scikit-image",
        "scipy",
        "shapely",
        "utm",
    ],
    extras_require={
        "extra": ["gdal"],
        "fractal": ["pyfftw"],
        "gbis": ["geoid"],                          # not available on pypi
        "isce": ["isce"],                           # not available on pypi
        "kite": ["kite"],
        "all": [
            "extra",
            "fractal",
            "gbis",
            "isce",
            "kite",
        ],
    },

    # package discovery
    packages=find_packages("src"),  # include all packages under src
    package_dir={"": "src"},        # tell distutils packages are under src
    entry_points={
        'console_scripts': [
            'mintpy = mintpy.__main__:main',
            'add.py = mintpy.cli.add:main',
            'asc_desc2horz_vert.py = mintpy.cli.asc_desc2horz_vert:main',
            'closure_phase_bias.py = mintpy.cli.closure_phase_bias:main',
            'dem_error.py = mintpy.cli.dem_error:main',
            'dem_gsi.py = mintpy.cli.dem_gsi:main',
            'diff.py = mintpy.cli.diff:main',
            'generate_mask.py = mintpy.cli.generate_mask:main',
            'geocode.py = mintpy.cli.geocode:main',
            'ifgram_inversion.py = mintpy.cli.ifgram_inversion:main',
            'image_math.py = mintpy.cli.image_math:main',
            'image_stitch.py = mintpy.cli.image_stitch:main',
            'info.py = mintpy.cli.info:main',
            'iono_tec.py = mintpy.cli.iono_tec:main',
            'load_data.py = mintpy.cli.load_data:main',
            'load_gbis.py = mintpy.cli.load_gbis:main',
            'local_oscilator_drift.py = mintpy.cli.local_oscilator_drift:main',
            'lookup_geo2radar.py = mintpy.cli.lookup_geo2radar:main',
            'mask.py = mintpy.cli.mask:main',
            'modify_network.py = mintpy.cli.modify_network:main',
            'multilook.py = mintpy.cli.multilook:main',
            'multi_transect.py = mintpy.multi_transect:main',
            'plate_motion.py = mintpy.cli.plate_motion:main',
            'plot_coherence_matrix.py = mintpy.cli.plot_coherence_matrix:main',
            'plot_network.py = mintpy.cli.plot_network:main',
            'plot_transection.py = mintpy.cli.plot_transection:main',
            'prep_aria.py = mintpy.cli.prep_aria:main',
            'prep_cosicorr.py = mintpy.cli.prep_cosicorr:main',
            'prep_fringe.py = mintpy.cli.prep_fringe:main',
            'prep_gamma.py = mintpy.cli.prep_gamma:main',
            'prep_gmtsar.py = mintpy.cli.prep_gmtsar:main',
            'prep_hyp3.py = mintpy.cli.prep_hyp3:main',
            'prep_isce.py = mintpy.cli.prep_isce:main',
            'prep_nisar.py = mintpy.cli.prep_nisar:main',
            'prep_roipac.py = mintpy.cli.prep_roipac:main',
            'prep_snap.py = mintpy.cli.prep_snap:main',
            'reference_date.py = mintpy.cli.reference_date:main',
            'reference_point.py = mintpy.cli.reference_point:main',
            'remove_hdf5_dset.py = mintpy.cli.remove_hdf5_dset:main',
            'remove_ramp.py = mintpy.cli.remove_ramp:main',
            's1ab_range_bias.py = mintpy.cli.s1ab_range_bias:main',
            'save_gbis.py = mintpy.cli.save_gbis:main',
            'save_gdal.py = mintpy.cli.save_gdal:main',
            'save_gmt.py = mintpy.cli.save_gmt:main',
            'save_hdfeos5.py = mintpy.cli.save_hdfeos5:main',
            'save_kite.py = mintpy.cli.save_kite:main',
            'save_kmz.py = mintpy.cli.save_kmz:main',
            'save_kmz_timeseries.py = mintpy.cli.save_kmz_timeseries:main',
            'save_qgis.py = mintpy.cli.save_qgis:main',
            'save_roipac.py = mintpy.cli.save_roipac:main',
            'smallbaselineApp.py = mintpy.cli.smallbaselineApp:main',
            'solid_earth_tides.py = mintpy.cli.solid_earth_tides:main',
            'spatial_average.py = mintpy.cli.spatial_average:main',
            'spatial_filter.py = mintpy.cli.spatial_filter:main',
            'subset.py = mintpy.cli.subset:main',
            'temporal_average.py = mintpy.cli.temporal_average:main',
            'temporal_derivative.py = mintpy.cli.temporal_derivative:main',
            'temporal_filter.py = mintpy.cli.temporal_filter:main',
            'timeseries2velocity.py = mintpy.cli.timeseries2velocity:main',
            'timeseries_rms.py = mintpy.cli.timeseries_rms:main',
            'tropo_gacos.py = mintpy.cli.tropo_gacos:main',
            'tropo_phase_elevation.py = mintpy.cli.tropo_phase_elevation:main',
            'tropo_pyaps3.py = mintpy.cli.tropo_pyaps3:main',
            'tsview.py = mintpy.cli.tsview:main',
            'unwrap_error_bridging.py = mintpy.cli.unwrap_error_bridging:main',
            'unwrap_error_phase_closure.py = mintpy.cli.unwrap_error_phase_closure:main',
            'view.py = mintpy.cli.view:main',
        ]
    },

    # data files
    include_package_data=True,
    package_data={
        "mintpy": [
            "data/*.js",
            "data/*.png",
            "data/colormaps/*.cpt",
            "defaults/*.cfg",
            "defaults/*.yaml",
        ],
    },
    zip_safe=False,
)
