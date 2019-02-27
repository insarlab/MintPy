[![Language](https://img.shields.io/badge/python-3.5%2B-blue.svg)](https://www.python.org/)
[![Latest version](https://img.shields.io/badge/latest%20version-v1.0.0--dev-green.svg)](https://github.com/yunjunz/PySAR/blob/master/docs/download.md)
[![License](https://img.shields.io/badge/license-GPL-yellow.svg)](https://github.com/yunjunz/PySAR/blob/master/LICENSE)
[![Forum](https://img.shields.io/badge/forum-Google%20Group-orange.svg)](https://groups.google.com/forum/#!forum/py-sar)
       
## InSAR time series analysis in Python
   
PySAR is an open-source package in Python for InSAR (Interferometric Synthetic Aperture Radar) time series analysis. It reads the stack of interferograms (coregistered and unwrapped) in ISCE, Gamma or ROI_PAC format, and produces three dimensional (2D in space and 1D in time) ground displacement. It includes a routine time series analysis (pysarApp.py) and some independent toolbox.      
   

### 1. [Download](https://github.com/yunjunz/PySAR/blob/master/docs/download.md)    


### 2. [Installation](https://github.com/yunjunz/PySAR/blob/master/docs/installation.md)   
    
   
### 3. Running PySAR

PySAR reads a stack of interferograms (unwrapped interferograms, coherence, wrapped interferograms and connecting components from SNAPHU if available) and the geometry files (DEM, lookup table, etc.). You need to give the path to where the files are and PySAR takes care of the rest!   
 
    pysarApp.py -h   #see the processing options.   
    pysarApp.py -H   #see the default template options with explanation.   
    pysarApp.py -g   #generate a default template file and see the detailed settings.   

#### [Example](https://github.com/yunjunz/PySAR/wiki/Example): Fernandina volcano, Galápagos with Sentinel-1 data    

    wget https://zenodo.org/record/2571980/files/FernandinaSenDT128.tar.xz
    tar -xvJf FernandinaSenDT128.tar.xz
    cd FernandinaSenDT128/PYSAR
    pysarApp.py FernandinaSenDT128.txt   

<p align="left">
  <img width="600" src="https://github.com/yunjunz/PySAR/blob/master/docs/resources/images/FernandinaSenDT128_POI.jpg">
</p>    

Inside pysarApp.py, it reads the unwrapped interferograms, references all of them to the same coherent pixel (reference point), calculates the phase closure and estimates the unwrapping errors (if it has been asked for), inverts the network of interferograms into time-series, calculates a parameter called "temporal coherence" which can be used to evaluate the quality of inversion, corrects local oscillator drift (for Envisat only), corrects stratified tropospheric delay (using pyaps or phase-elevation-ratio approach), removes phase ramps (if it has been asked for), corrects DEM error,... and finally estimates the velocity.   

Check **./PIC** folder for auto-generated figures. More details about this test data are in [here](https://github.com/yunjunz/PySAR/wiki/Example).     

Check **/exmaple/input_files** folder for example template file for different InSAR processors.

#### Some useful scripts for information and visualization:   

    info.py                    #check HDF5 file structure and metadata
    view.py                    #2D map view
    tsview.py                  #1D point time-series (interactive)   
    transect.py                #1D profile (interactive)
    plot_coherence_matrix.py   #plot coherence matrix for one pixel (interactive)
    plot_network.py            #plot network configuration of the dataset    
    save_kmz.py                #generate Google Earth KMZ file in raster image
    save_kmz_timeseries.py     #generate Goodle Earth KMZ file in points for time-series (interactive)
     
#### Build your own processing recipe   

PySAR is a toolbox with a lot of individual utility scripts, highly modulized in python. Check its documentation or simply run it with -h to see its usage, you could build your own customized processing recipe! Here is an example to compare the velocities estimated from displacement time-series with different tropospheric delay corrections: [link](https://github.com/yunjunz/PySAR/blob/master/sh/compare_velocity_with_diff_tropcor.sh)

   
### 4. Documentation
   
- Manual: [PDF](https://github.com/yunjunz/PySAR/blob/master/docs/Manual-0.4.0_201803.pdf), [HTML](https://github.com/yunjunz/PySAR/blob/master/docs/Manual-0.4.0_201803.html.zip), [Workshop](https://miami.box.com/v/pysar-workshop-2017-miami)     
- Wiki: Check our [Github Wiki](https://github.com/yunjunz/PySAR/wiki) to see the example data, paper references, file naming convention and more.
   
### 5. Google Group

Join our google group [https://groups.google.com/forum/#!forum/py-sar](https://groups.google.com/forum/#!forum/py-sar) to ask questions, get notice of latest features pushed to you!

### Contributors    

* Zhang Yunjun
* Heresh Fattahi
* Falk Amelung
* Scott Baker
* Joshua Zahner
* Alfredo Terreco
* David Grossman
* Yunmeng Cao
* [_other community members_](https://github.com/yunjunz/PySAR/graphs/contributors)


### License
[![FOSSA Status](https://app.fossa.io/api/projects/git%2Bgithub.com%2Fyunjunz%2FPySAR.svg?type=large)](https://app.fossa.io/projects/git%2Bgithub.com%2Fyunjunz%2FPySAR?ref=badge_large)
