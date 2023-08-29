MintPy use [pyKML](https://pythonhosted.org/pykml/) to generate KMZ (Keyhole Markup Zip) files for easy offline viewing in [Google Earth](https://www.google.com/earth/) via `save_kmz_timeseries.py` and `save_kmz.py` script. Below are screenshots of the displacement time-series and average velocity of [Fernandina volcano estimated from Sentinel-1 data](demo_dataset.md).

### 1. Displacement time-series ###

`save_kmz_timeseries.py` takes 3D displacement time-series file and outputs a KMZ file with interactive time-seires plot.

<p align="center">
  <img src="https://yunjunzhang.files.wordpress.com/2019/02/fernandinasendt128_ge-1.png">
</p>

[Download KMZ file](https://miami.box.com/v/FernandinaSenDT128TS)

### 2. Raster image ###

`save_kmz.py` takes any 2D matrix and outputs a KMZ file with a overlay image.

<p align="center">
  <img src="https://yunjunzhang.files.wordpress.com/2019/02/vel_fernandinasendt128_ge.png">
</p>

[Download KMZ file](https://miami.box.com/v/FernandinaSenDT128VEL)

### Notes for developers ###

save_kmz_timeseries.py embeds a [dygraphs](http://dygraphs.com) javascript for interactive plot of the time-series deformation at each point. This allows the user to select any placemark onscreen to display the time-series data in an interactive chart. Placemarks are colored based on the velocity.

The script also use the [regions KML feature](https://developers.google.com/kml/documentation/regions) to support very large datasets without sacrificing resolution. It divides the data matrix into regionalized boxes, nests them using network links so that Google Earth could load them in a "smart" way.

**Alert: for very large datasets, the default settings are not generic due to the various computer memories, data sizes and different preferred details. The user is highly recommended to read the following to understand how the regions feature works and adjust parameters accordingly.**

1. Level of Detail (LOD)

The script samples the input 3D dataset at 3 levels of details by default (`--steps` option): low-, moderate- and high-resolution. Each LOD is displayed at a different zoom-level (`--level-of-details` option) within Google Earth. On startup, the low-resolution LOD is displayed; then at ~20km in altitude, the low-resolution LOD disappears and the moderate-resolution LOD becomes visible; similarly, the high-resolution LOD shows at ~10km. In this way, Google Earth only has to load as many placemark as are on the screen currently. This LOD strategy drastically increases performance.

The low- and moderate-resolution LODs cover the entire region, while the high-resolution LOD covers only the actively deforming regions. These regions (red boxes below) are currently identified as boxes having >20% pixels with velocity magnitude > the global velocity median absolute deviation (`mintpy.save_kmz_timeseries.get_boxes4deforming_area`).

<p align="center">
  <img src="https://yunjunzhang.files.wordpress.com/2020/03/defo_area.png">
</p>

2. Region-based Network Links

To further increase performance, the script splits each LOD into 300x300 point subsets known as regions. Each region is written to a separate KML file and referenced via a "Network Link" in a root level KML file. Based on whether the bounding box of each region is currently on screen or not, Google Earth will load them accordingly.
