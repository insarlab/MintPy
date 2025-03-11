Displacement time-series can be visualized in QGIS InSAR Explorer plugin. The plugin supports time-series data as GRD files or shapefile.

### Setup QGIS and InSAR Explorer ###
1. Download and Install [QGIS](https://qgis.org/en/site/) if you have not done so.
2. Install InSAR Explorer:
  - Install the latest stable version of the [InSAR Explorer plugin](https://plugins.qgis.org/plugins/insar_explorer-dev/) via “Plugins -> Manage and Install Plugins.”
  - Alternatively, download the plugin as a *.zip file from the [InSAR Explorer GitHub repository](https://github.com/luhipi/insar-explorer) and install it through “Plugins -> Manage and Install Plugins -> Install from ZIP.”
3. Launch InSAR Explorer: Access it from the toolbar or through “Plugins -> InSAR Explorer -> InSAR Explorer.”


### Using GRD files ###
1. Export MintPy results to GRD files compatible with the QGIS InSAR Explorer plugin using `save_explorer.py`.

  ```
    $ save_explorer.py geo_timeseries.h5 -v geo_velocity.h5 -o geo_maskTempCoh.h5 -o timeseries/
   ```

2. Load data in QGIS: Open one of the exported GRD files (for example `geo_velocity_mm.h5`), in QGIS.
3. Launch InSAR Explorer and Click on any point to plot the time series.

### Using shapefile ###
1. Export to shapefile can be done using `save_qgis.py` script.
2. Load the shapefile in QGIS.
3. Launch InSAR Explorer and Click on any point to plot the time series.


### More information ###
For more details on using the plugin, please refer to the [InSAR Explore documentation](https://insar-explorer.readthedocs.io/).
