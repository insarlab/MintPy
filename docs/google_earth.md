MintPy use [pyKML](https://pythonhosted.org/pykml/) to generate KMZ (Keyhole Markup Zip) files for easy offline viewing in [Google Earth](https://www.google.com/earth/) via `save_kmz_timeseries.py` and `save_kmz.py` script. Below are screenshots of the displacement time-series and average velocity of [Fernandina volcano estimated from Sentinel-1 data](example_dataset.md).    

### 1. Displacement time-series ###

`save_kmz_timeseries.py` takes 3D displacement time-series file and outputs a KMZ file with interactive time-seires plot.

```bash
save_kmz_timeseries.py geo_timeseries_ECMWF_ramp_demErr.h5 --vel geo_velocity_mask.h5 --tcoh geo_temporalCoherence.h5
```

<p align="center">
  <img src="https://yunjunzhang.files.wordpress.com/2019/02/fernandinasendt128_ge-1.png">
</p>

[Download KMZ file](https://miami.box.com/v/FernandinaSenDT128TS)

### 2. Raster image ###

`save_kmz.py` takes any 2D matrix and outputs a KMZ file with a overlay image.

```bash
save_kmz.py geo_velocity_masked.h5 --wrap --wrap-range -3 3
```

<p align="center">
  <img src="https://yunjunzhang.files.wordpress.com/2019/02/vel_fernandinasendt128_ge.png">
</p>

[Download KMZ file](https://miami.box.com/v/FernandinaSenDT128VEL)

### Notes for developers ###

save_kmz_timeseries.py takes the 3D HDF5 file and outputs a KMZ file at multiple levels of details (LODs). It subsets of the data into regionalized boxes, writes the required KML files, references the appropriate auxiliary resources, and zips all of the files together into a KMZ file. 

The script also embeds a javascript-based interactive plot of the timeseries deformation at each latitude-longitude point included in the original dataset. This allows for the user to select any of the onscreen placemarks and view a timeseries chart displaying the timeseries data from that unique geographic coordinate. The Placemarks are colored onscreen based on their deformation velocity, as computed from a velocity data file.

#### KML and KMZ Performance ####

KML files are capable of holding a very large amount of data, but earth viewer programs, such as Google Earth, often struggle to display extremely large datasets effectively. As such, a few tradeoffs and specialized file structures have been enabled within `save_kmz_timeseries.py` to increase Google Earth's performance with large data sets.

+ Multiple Levels of Detail

By default, MintPy will read in the provided dataset and then subset the data into three separate levels of details: a low resolution, a high-resolution, and a full-resolution. The low-resolution LOD contains 1 point for every 20x20 points in the original file. The high-resolution LOD contains 1 point for every 3x3 points in the original file. And the full-resolution LOD contains every point in the original file. 

Each LOD is displayed at a different zoom-level within Google Earth. On startup, the low-resolution LOD is displayed, while the high-resolution LOD becomes visible around 20km in altitude, and the full-resolution LOD at around 10km in altitude. This ensures that Google Earth only has to load as many Placemark as are on the screen currently, which drastically increases performance as fewer Placemarks are onscreen at higher zoom levels.

The full-resolution LOD is presently calculated and presented only for those actions showing signs of active deformation so as to further increase performance.

+ Regionalized Network Links

To further increase performance, MintPy splits each LOD into 300x300 point subsets known as regions. Each region is written to a separate KML file, and are then referenced via a "Network Link" in a master level KML file. Google Earth specifically has the ability to load conditionally load Network Link elements based on whether or not the coordinates dictating their bounding box are on screen at a given moment, so this method ensures that, at high zoom-levels, only as many placemarks as are onscreen at the time are loaded.

#### Modify save_kmz_timeseries.py ####

If you wish to add or modify the levels of detail or other aspects of the KMZ generation process, the function `generate_network_link(inps, ts_obj, box_list, step, lod, output_file=None)` handles the generation and reference linking of the required files. Be sure to handle the initial reading of the data and filesystem cleanup afterward.
