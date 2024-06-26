MintPy mainly uses attribute names from [ROI_PAC](http://www.geo.cornell.edu/eas/PeoplePlaces/Faculty/matt/pub/winsar/InSAR_textbook_for_web_2014.pdf), with some additional self-generated attributes.

### Required attributes ###

If using ROI_PAC as the InSAR processor, both **baseline parameter RSC** file (i.e. *100416-100901_baseline.rsc*) and **basic metadata file** (i.e. *filt_100416-100901-sim_HDR_4rlks_c10.unw.rsc*) will be imported into MintPy. The following attributes for each interferogram are required in order to run MintPy:

+  LENGTH = number of rows.
+  WIDTH = number of columns.
+  X/Y_FIRST = (for geocoded product) Longitude/easting/X and latitude/northing/Y coordinate in degrees/meters of the upper left corner of the first pixel.
+  X/Y_STEP = (for geocoded product) Ground resolution in degrees/meters in X/Y direction.
+  X/Y_UNIT = (for geocoded product) Coordinate unit in X/Y direction: degrees or meters.
+  LAT/LON_REF1/2/3/4 = Latitude/northing and longitude/easting at corner 1/2/3/4 (in degrees or meters), used in save_unavco, PyAPS (DEM file in radar coord), not accurate; number named in order of first line near/far range, last line near/far range.
+  WAVELENGTH = Radar wavelength in meters.
+  RANGE_PIXEL_SIZE = Slant range pixel size (search for pixel_ratio to convert to ground size, in m), used in dem_error, incidence_angle, multilook, transect.
+  EARTH_RADIUS = Best fitting spheroid radius in meters, used in dem_error, incidence_angle, convert2mat.
+  CENTER_LINE_UTC = Time at middle of interferogram in seconds, used in tropo correction using PyAPS.
+  HEIGHT = Height of satellite in meters, used in dem_error, incidence_angle, convert2mat.
+  STARTING_RANGE = Distance from satellite to first ground pixel in meters, used in incidence_angle calculation
+  PLATFORM = satellite/sensor name, used in Local Oscillator Drift correction for Envisat.
+  ORBIT_DIRECTION = ascending, or descending.
+  ALOOKS/RLOOKS = multilook number in azimuth/range direction, used in weighted network inversion.

The following attributes vary for each interferogram:

+  DATE12 = (date1)-(date2), reference - secondary date of interferogram in 6 digit number.
+  P_BASELINE_TOP_HDR = Perpendicular baseline at top (first line) of interferogram in meters.
+  P_BASELINE_BOTTOM_HDR = Perpendicular baseline at bottom (last line) of interferogram in meters.

### Optional attributes ###

+  ANTENNA_SIDE = -1 for right looking radar, used in save_unavco
+  AZIMUTH_PIXEL_SIZE = Azimuth pixel size at orbital altitude (multiply by Re/(Re+h) for ground size (m), where Re is the local earth radius), used in baseline_error/trop and multilook.
+  HEADING = Spacecraft heading at peg point (degrees), measured from the north with clock-wise as positive, used in asc_desc
+  PRF = Pulse repetition frequency (Hz), used in save_unavco

### Self-generated attributes ###

+  FILE_TYPE = file type.
    - for HDF5 files, it's the root level dataset name, such as `velocity, timeseries, ifgramStack, temporalCoherence, mask, HDFEOS, dem, coherence, etc.`;`
    - for binary files, it's the file extension name, such as `.unw, .cor, .int, .amp, .mli, .dem, .hgt, .unw.conncomp, .UTM_TO_RDC, .trans, etc.`, except for ISCE geometry files, which is the file base name such as `hgt, lat, lon, los, shadowMask, incLocal`.
+  FILE_PATH = absolute file path
+  PROCESSOR = processing software, i.e. isce, aria, snap, gamma, roipac etc.
+  DATA_TYPE = data type, i.e. float32, int16, etc., for isce product read using GDAL
+  BYTE_ORDER = order of bytes, or endianness, for binary files, i.e. big-endian or little-endian. Most InSAR processors use little-endian, except for Gamma, which sometimes use big-endian.
+  BANDS = number of bands, for binary file I/O.
+  INTERLEAVE = band interleave type, i.e. BSQ, BIL, BIP for binary file I/O.
+  NO_DATA_VALUE = No data value, value that should be ignored.
+  UNIT = data unit, i.e. m, m/yr, radian, and 1 for file without unit, such as coherence [[source]](https://github.com/insarlab/MintPy/blob/main/src/mintpy/objects/stack.py#L75)
+  REF_DATE = reference date
+  REF_X/Y = column/row of the reference point
+  REF_LAT/LON = latitude/northing and longitude/easting of the reference point (in degrees or meters)
+  SUBSET_XMIN/XMAX/YMIN/YMAX = start/end column/row number of subset in the original coverage
+  MODIFICATION_TIME = dataset modification time, exists in ifgramStack.h5 file for 3D dataset, used for "--update" option of unwrap error corrections.
+  NCORRLOOKS = number of independent looks, as explained in [SNAPHU](https://web.stanford.edu/group/radar/softwareandlinks/sw/snaphu/snaphu.conf.full)
+  UTM_ZONE = [UTM zone](https://docs.up42.com/data/reference/utm#utm-wgs84), comprises a zone number and a hemisphere, e.g. 11N, 60S, for geocoded file with UTM projection only.
+  EPSG = EPSG code for coordinate systems, for geocoded files only. Check [here](https://docs.up42.com/data/reference/utm#utm-wgs84) for its relationship with UTM zone.
+  CENTER_INCIDENCE_ANGLE = incidence angle in degrees at the scene center, read from the 2D incidence angle matrix, for isce2 files only.

### Reference ###

+ Pritchard et al., (2014), Open-source software for geodetic imaging: ROI_PAC for InSAR and pixel tracking, pp 44-48. [PDF](http://www.geo.cornell.edu/eas/PeoplePlaces/Faculty/matt/pub/winsar/InSAR_textbook_for_web_2014.pdf)
