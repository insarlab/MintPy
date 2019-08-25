## Tutorials in Jupyter Notebook for InSAR time series analysis

This tutorial walks through the various processing steps of InSAR time series analysis using MintPy.

### Contents ###

1. Small baseline time series analysis with `smallbaselineApp`

   - ISCE/topsStack + MintPy (Sentinel-1 on Fernandina volcano, Galápagos): [nbviewer](https://nbviewer.jupyter.org/github/insarlab/MintPy/blob/master/docs/tutorials/smallbaselineApp.ipynb)
   - ARIA + MintPy (Sentinel-1 on San Francisco): [nbviewer](https://nbviewer.jupyter.org/github/insarlab/MintPy/blob/master/docs/tutorials/smallbaselineApp_aria.ipynb)

2. Visualizations   

   - Interactive time-series with [tsview](https://nbviewer.jupyter.org/github/insarlab/MintPy/blob/master/docs/tutorials/tsview.ipynb)
   - Interactive coherence matrix with [plot_coherence_matrix](https://nbviewer.jupyter.org/github/insarlab/MintPy/blob/master/docs/tutorials/plot_coherence_matrix.ipynb)
   - Interactive transection with [plot_transection](https://nbviewer.jupyter.org/github/insarlab/MintPy/blob/master/docs/tutorials/plot_transection.ipynb)
   - Google Earth [doc](https://mintpy.readthedocs.io/en/latest/google_earth/)

### Useful links ###

1. Single interferogram processing with [ISCE2](https://github.com/isce-framework/isce2-docs/tree/master/Notebooks)

   - Sentinel-1 TOPS mode SAR data with [topsApp](https://nbviewer.jupyter.org/github/isce-framework/isce2-docs/blob/master/Notebooks/TOPS/Tops.ipynb)
   - StripMap mode SAR data with [stripmapApp](https://nbviewer.jupyter.org/github/isce-framework/isce2-docs/blob/master/Notebooks/Stripmap/stripmapApp.ipynb)

2. Manipulate ARIA standard InSAR products with [ARIA-tools](https://github.com/aria-tools/ARIA-tools-docs)

   - Downloading GUNW products using [ariaDownload](https://nbviewer.jupyter.org/github/aria-tools/ARIA-tools-docs/blob/master/JupyterDocs/ariaDownload/ariaDownload_tutorial.ipynb)
   - Preparing GUNW products for time series analysis using [ariaTSsetup](https://nbviewer.jupyter.org/github/aria-tools/ARIA-tools-docs/blob/master/JupyterDocs/ariaTSsetup/ariaTSsetup_tutorial.ipynb)
