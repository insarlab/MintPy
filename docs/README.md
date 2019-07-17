## MintPy

[![Language](https://img.shields.io/badge/python-3.6%2B-blue.svg)](https://www.python.org/)
[![Docs Status](https://readthedocs.org/projects/mintpy/badge/?version=latest)](https://mintpy.readthedocs.io/?badge=latest)
[![CircleCI](https://img.shields.io/circleci/build/github/insarlab/MintPy.svg?color=green&logo=circleci)](https://circleci.com/gh/insarlab/MintPy)
[![Latest version](https://img.shields.io/badge/latest%20version-v1.1.2-yellowgreen.svg)](https://github.com/insarlab/MintPy/releases)
[![License](https://img.shields.io/badge/license-GPL-yellow.svg)](https://github.com/insarlab/MintPy/blob/master/LICENSE)
[![Forum](https://img.shields.io/badge/forum-Google%20Group-orange.svg)](https://groups.google.com/forum/#!forum/mintpy)

The Miami INsar Time-series software in PYthon (MintPy) is an open-source package for Interferometric Synthetic Aperture Radar time series analysis. It reads the stack of interferograms (coregistered and unwrapped) in [ISCE](https://github.com/isce-framework/isce2), [SNAP](http://step.esa.int/), Gamma or ROI_PAC format, and produces three dimensional (2D in space and 1D in time) ground surface displacement. It includes a routine time series analysis (`smallbaselineApp.py`) and some independent toolbox.

This package was called PySAR before version 1.1.1. For version 1.1.2 and onward, we use MintPy instead.

### 1. [Installation](./installation.md) ###

### 2. Running MintPy ###

MintPy reads a stack of interferograms (unwrapped interferograms, coherence, wrapped interferograms and connecting components from SNAPHU if available) and the geometry files (DEM, lookup table, etc.). You need to give the path to where the files are and MintPy takes care of the rest!

```bash
smallbaselineApp.py                         #run with default template 'smallbaselineApp.cfg'
smallbaselineApp.py <custom_template>       #run with default and custom templates
smallbaselineApp.py -h / --help             #help
smallbaselineApp.py -H                      #print    default template options
smallbaselineApp.py -g                      #generate default template if it does not exist
smallbaselineApp.py -g <custom_template>    #generate/update default template based on custom template

# Run with --start/stop/dostep options
smallbaselineApp.py GalapagosSenDT128.template --dostep velocity  #run at step 'velocity' only
smallbaselineApp.py GalapagosSenDT128.template --end load_data    #end after step 'load_data'
```

#### [Example](./example_dataset.md) on Fernandina volcano, Galápagos with Sentinel-1 data ####

```
wget https://zenodo.org/record/2748487/files/FernandinaSenDT128.tar.xz
tar -xvJf FernandinaSenDT128.tar.xz
cd FernandinaSenDT128/MintPy
smallbaselineApp.py ${MINTPY_HOME}/docs/examples/input_files/FernandinaSenDT128.txt
```

<p align="left">
  <img width="600" src="https://yunjunzhang.files.wordpress.com/2019/06/fernandinasendt128_poi.jpg">
</p>

Inside smallbaselineApp.py, it reads the unwrapped interferograms, references all of them to the same coherent pixel (reference point), calculates the phase closure and estimates the unwrapping errors (if it has been asked for), inverts the network of interferograms into time-series, calculates a parameter called "temporal coherence" which can be used to evaluate the quality of inversion, corrects local oscillator drift (for Envisat only), corrects stratified tropospheric delay (using pyaps or phase-elevation-ratio approach), removes phase ramps (if it has been asked for), corrects DEM error,... and finally estimates the velocity.

Check **./pic** folder for auto-generated figures. More details about this test data are in [here](./example_dataset.md).

#### 2.1 Data visualization ####

Below are some useful scripts for data information and visulization.

```bash
info.py                    #check HDF5 file structure and metadata
view.py                    #2D map view
tsview.py                  #1D point time-series (interactive)   
plot_coherence_matrix.py   #plot coherence matrix for one pixel (interactive)
plot_network.py            #plot network configuration of the dataset    
plot_transection.py        #1D profile (interactive)
save_kmz.py                #generate Google Earth KMZ file in raster image
save_kmz_timeseries.py     #generate Goodle Earth KMZ file in points for time-series (interactive)
```

#### 2.2 Customized processing recipe: [example](https://github.com/insarlab/MintPy/blob/master/sh/compare_velocity_with_diff_tropo.sh) ####

MintPy is a toolbox with a lot of individual utility scripts, highly modulized in python. Check its documentation or simply run it with -h to see its usage, you could build your own customized processing recipe! Here is an example to compare the velocities estimated from displacement time-series with different tropospheric delay corrections: [link](https://github.com/insarlab/MintPy/blob/master/sh/compare_velocity_with_diff_tropo.sh)

### 3. [Documentation](https://mintpy.readthedocs.io/) ###

+ [Tutorials in Jupyter Notebook](./tutorials/README.md)
+ [Example datasets](./example_dataset.md)
+ [Example template files for InSAR processors](./examples/input_files/README.md)
+ [Google Earth KMZ file](./google_earth.md)
+ [Paper figures in Jupyter Notebook](./paper/README.md)

### 4. User Forum ###

Join our google group [https://groups.google.com/forum/#!forum/mintpy](https://groups.google.com/forum/#!forum/mintpy) to ask questions, get notice of latest features pushed to you!

### Contributors ###

* Zhang Yunjun
* Heresh Fattahi
* Falk Amelung
* Scott Baker
* Joshua Zahner
* Alfredo Terreco
* David Grossman
* Yunmeng Cao
* [_other community members_](https://github.com/insarlab/MintPy/graphs/contributors)
