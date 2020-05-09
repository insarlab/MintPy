[![Language](https://img.shields.io/badge/python-3.5%2B-blue.svg)](https://www.python.org/)
[![Docs Status](https://readthedocs.org/projects/mintpy/badge/?version=latest)](https://mintpy.readthedocs.io/?badge=latest)
[![CircleCI](https://img.shields.io/circleci/build/github/insarlab/MintPy.svg?color=green&logo=circleci)](https://circleci.com/gh/insarlab/MintPy)
[![Latest version](https://img.shields.io/badge/latest%20version-v1.2.2-yellowgreen.svg)](https://github.com/insarlab/MintPy/releases)
[![License](https://img.shields.io/badge/license-GPLv3-yellow.svg)](https://github.com/insarlab/MintPy/blob/master/LICENSE)
[![Forum](https://img.shields.io/badge/forum-Google%20Group-orange.svg)](https://groups.google.com/forum/#!forum/mintpy)
[![Citation](https://img.shields.io/badge/doi-10.1016%2Fj.cageo.2019.104331-blue)](https://doi.org/10.1016/j.cageo.2019.104331)

## MintPy ##

The Miami INsar Time-series software in PYthon (MintPy) is an open-source package for Interferometric Synthetic Aperture Radar time series analysis. It reads the stack of interferograms (coregistered and unwrapped) in [ISCE](https://github.com/isce-framework/isce2), [ARIA](https://github.com/aria-tools/ARIA-tools), [FRInGE](https://github.com/isce-framework/fringe), [SNAP](http://step.esa.int/), [GAMMA](https://www.gamma-rs.ch/no_cache/software.html) or ROI_PAC format, and produces three dimensional (2D in space and 1D in time) ground surface displacement in line-of-sight direction. It includes a routine time series analysis (`smallbaselineApp.py`) and some independent toolbox.

This package was called PySAR before version 1.1.1. For version 1.1.2 and onward, we use MintPy instead. 

This is research code provided to you "as is" with NO WARRANTIES OF CORRECTNESS. Use at your own risk.

### 1. [Installation](./installation.md) ###

### 2. Running MintPy ###

#### 2.1 Running routine workflow `smallbaselineApp.py` ####

MintPy reads a stack of interferograms (unwrapped interferograms, coherence and connecting components from SNAPHU if available) and the geometry files (DEM, lookup table, incidence angle, etc.). You need to give the [path to where the files are](dir_structure.md) and MintPy takes care of the rest! 

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

#### [Example](./demo_dataset.md) on Fernandina volcano, Galápagos with Sentinel-1 data ####

```bash
wget https://zenodo.org/record/3635245/files/FernandinaSenDT128.tar.xz
tar -xvJf FernandinaSenDT128.tar.xz
cd FernandinaSenDT128/mintpy
smallbaselineApp.py ${MINTPY_HOME}/docs/examples/input_files/FernandinaSenDT128.txt
```

<p align="left">
  <img width="600" src="https://yunjunzhang.files.wordpress.com/2019/06/fernandinasendt128_poi.jpg">
</p>

Inside smallbaselineApp.py, it reads the unwrapped interferograms, references all of them to the same coherent pixel (reference point), calculates the phase closure and estimates the unwrapping errors (if it has been asked for), inverts the network of interferograms into time-series, calculates the temporal coherence to evaluate the quality of inversion, corrects local oscillator drift (for Envisat only), corrects stratified tropospheric delay (using global atmospheric models or phase-elevation-ratio approach), removes phase ramps (if it has been asked for), corrects DEM error,... and finally estimates the velocity. 

Configuration parameters for each step are initiated with default values in a customizable text file [**smallbaselineApp.cfg**](../mintpy/defaults/smallbaselineApp.cfg). Results are plotted in **./pic** folder.

#### 2.2 Data visualization ####

Below are some useful scripts for data information and visulization.

```bash
info.py                    #check HDF5 file structure and metadata
view.py                    #2D map view
tsview.py                  #1D point time-series (interactive)   
plot_coherence_matrix.py   #plot coherence matrix for one pixel (interactive)
plot_network.py            #plot network configuration of the dataset    
plot_transection.py        #plot 1D profile along a line of a 2D matrix (interactive)
save_kmz.py                #generate Google Earth KMZ file in raster image
save_kmz_timeseries.py     #generate Goodle Earth KMZ file in points for time-series (interactive)
```

#### 2.3 Customized processing recipe: [example](https://github.com/insarlab/MintPy/blob/master/sh/compare_velocity_with_diff_tropo.sh) ####

MintPy is a toolbox with a lot of individual utility scripts. Check its documentation or simply run the script with `-h / --help` to see its usage, you could build your own customized processing recipe! Here is an example to compare the velocities estimated from displacement time-series with different tropospheric delay corrections: [link](https://github.com/insarlab/MintPy/blob/master/sh/compare_velocity_with_diff_tropo.sh)

#### 2.4 Build on top of mintpy python module ####

MintPy is modulized in Python with a lot of utilities class and functions and well commented in the code level. Users who are familiar with Python could build their own functions and modules on top of [`mintpy.objects`](../mintpy/objects) and [`mintpy.utils`](../mintpy/utils). However, we don't have a complete API document website yet (maybe you can contribute this!). Below is an example of reading the 3D matrix of displacement time-series from an HDF5 file.

```python
from mintpy.utils import readfile
ts_data, meta = readfile.read('timeseries_ERA5_ramp_demErr.h5')
```

### 3. [Documentation](https://mintpy.readthedocs.io/) ###

Algorithms implemented in the software are described in details in [Yunjun et al. (2019)](https://doi.org/10.1016/j.cageo.2019.104331).

+ [Tutorials in Jupyter Notebook](https://github.com/insarlab/MintPy-tutorial)
+ [Example datasets](./demo_dataset.md)
+ [Example data directory](./dir_structure.md)
+ [Example template files for InSAR processors](./examples/input_files/README.md)
+ [Google Earth KMZ file](./google_earth.md)

### 4. User Forum ###

Join our google group [https://groups.google.com/forum/#!forum/mintpy](https://groups.google.com/forum/#!forum/mintpy) to ask questions, get notice of latest features pushed to you!

### 5. Citing this work ###

Yunjun, Z., H. Fattahi, F. Amelung (2019), Small baseline InSAR time series analysis: Unwrapping error correction and noise reduction, _Computers & Geosciences_, _133_, 104331, doi:[10.1016/j.cageo.2019.104331](https://doi.org/10.1016/j.cageo.2019.104331), [arXiv](https://eartharxiv.org/9sz6m/), [data & figures](https://github.com/geodesymiami/Yunjun_et_al-2019-MintPy).

In addition to the above, we recommend that you cite the original publications that describe the algorithms used in your specific analysis. They are noted briefly in the [default template file](../mintpy/defaults/smallbaselineApp.cfg) and listed in [here](./references.md).

### Contributors ###

* Zhang Yunjun
* Heresh Fattahi
* Falk Amelung
* Scott Baker
* Joshua Zahner
* Alfredo Terreco
* David Grossman
* Emre Havazli
* Yunmeng Cao
* Andre Theron
* [Other community members.](https://github.com/insarlab/MintPy/graphs/contributors) [[Contribute](./CONTRIBUTING.md)]
