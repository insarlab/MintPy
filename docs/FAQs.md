## Frequently Asked Questions

### 1. What's the sign convention of the line-of-sight data?

For line-of-sight (LOS) phase in the unit of radians, i.e. 'unwrapPhase' dataset in `ifgramStack.h5` file, positive value represents motion away from the satellite. We assume the "date1_date2" format for the interferogram with "date1" being the earlier acquisition.

For LOS displacement (velocity) in the unit of meters (m/yr), i.e. 'timeseries' dataset in `timeseries.h5` file, positive value represents motion toward the satellite (uplift for pure vertical motion).

### 2. How to prepare the input for MintPy if I am using currently un-supported InSAR software?

The input of MintPy routine workflow (`smallbaselineApp.py`) is a stack of unwrapped interferograms. For "stack", we mean all the interferograms (unwrapped phase and spatial coherence) and geometries (DEM, incidence angle, etc.) have the same spatial extent and same spatial resolution, either in geo-coordinates or radar (range-doppler) coordinates. The input has 2 components: data and attributes.

All inputs are saved into the following HDF5 files during the data loading process in MintPy:

+ `inputs/ifgramStack.h5` for stack of interferograms in geo-/radar-coordinates and attributes.
+ `inputs/geometryGeo.h5` for geometry data in geo-coordinates and attributes.
+ `inputs/geometryRadar.h5` for geometry data in radar-coordinates and attributes.

The currently supported formats are:

+ `ISCE-2` stack processors (`topsStack`, `stripmapStack` and `alosStack`)
+ `ARIA` products pre-processed using ARIA-tools
+ `SNAP` products produced using a preliminary workflow here: https://github.com/insarlab/MintPy/wiki/SNAP-input-data

If not using the above software workflows, below is a brief guide for preparation before running `smallbaseelineApp.py`:

#### Data

For dataset in geo-coordinates [recommended]:

+ for each interferogram, the unwrapped phase
+ for each interferogram, the spatial coherence
+ for each interferogram, the connected components from phase unwrapping (produced by SNAPHU) [optional]
+ for each stack, the DEM (Digital Elevation Model)
+ for each stack, the LOS incidence angle [optional]
+ for each stack, the LOS azimuth angle [optional]
+ for each stack, shadow mask [optional]
+ for each stack, water mask [optional]

For dataset in radar-coordinates, the extra lookup table file(s) is required (_e.g._ lat/lon.rdr for `ISCE-2`, sim_\*.UTM_TO_RDC for `Gamma`, geo_\*.trans for `ROI_PAC`).

All the files above should be in the same spatial extent and same spatial resolution (except for the lookup table in geo-coordinates from Gamma/ROI_PAC). If they are not (e.g. different row/column number, different spatial extent in terms of SNWE, different spatial resolution, etc.), the easiest way is to geocode them with the same output spatial extent and same output spatial resolution.

MintPy read data files via `mintpy.utils.readfile.read()`. It supports the following two types of file formats:

+ binary files with metadata files in the format of `ROIPAC .rsc`, `ISCE .xml`, `Gamma .par`, `ENVI .hdr` and `GDAL .vrt` via `numpy`.
+ GeoTiff and GRD files via `GDAL`.

The same function is used in the visualization as well, thus, you could run `view.py` on your files to check. If they are displayed correctly, they will be read properly by MintPy.

Note that MintPy assumes **date1_date2** convention for interferograms (date1 < date2) with positive phase value in radians for motion away from the satellite.

#### Attributes

For each data file, MintPy requires some attributes/metadata in `ROI_PAC .rsc` format as described [here](https://mintpy.readthedocs.io/en/latest/api/attributes/). The optional attributes are highly recommend. The interferogram specific attributes (DATE12, P_BASELINE_TOP/BOTTOM_HDR) is not needed for geometry files (DEM, incidence angle, etc.).

For the supported InSAR software workflows, we prepare these metadata via `prep_isce.py`, `prep_gamma.py`, etc., to read their native metadata, convert and write them into `ROI_PAC .rsc` style for each data file.

#### Recommendations

Option 1: If MintPy could read your data files, e.g. via testing with `view.py`, we recommend leveraging the existing `readfile.read()` by:

+ writing your own `prep_*.py` to extract the metadata, convert and write them into `ROI_PAC .rsc` convention for each data file.
+ specifying the path of each type of data files in the template file and loading them via `load_data.py`.

This is very similar to the existing ISCE/topsStack + MintPy example on Fernandina volcano.

Option 2: Write a independent script to handle the data reading, attributes/metadata preparation and writing into HDF5 file, like `prep_aria.py`, as shown in the ARIA + MintPy example on San Francisco Bay.
