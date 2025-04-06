We support output geocoded displacement time-series product into [HDF-EOS5](http://hdfeos.org) format via `save_hdfeos5.py`. This is designed to easily share the InSAR time-series product to the broader community.

```bash
save_hdfeos5.py geo_timeseries_ERA5_ramp_demErr.h5 --tc geo_temporalCoherence.h5 --asc geo_avgSpatialCoh.h5 -m geo_maskTempCoh.h5 -g geo_geometryRadar.h5
save_hdfeos5.py timeseries_ERA5_ramp_demErr.h5     --tc temporalCoherence.h5     --asc avgSpatialCoh.h5     -m maskTempCoh.h5     -g inputs/geometryGeo.h5
```

### 1. File structure ###

```
/                             Root level group
Attributes                    metadata in dict
/HDFEOS/GRIDS/timeseries      timeseries group
    /observation
        /displacement         3D array of float32 in size of (n, l, w) in meter
        /date                 1D array of string  in size of (n,     ) in YYYYMMDD format.
        /bperp                1D array of float32 in size of (n,     ) in meter
    /quality
        /mask                 2D array of bool_   in size of (   l, w).
        /temporalCoherence    2D array of float32 in size of (   l, w).
        /avgSpatialCoherence  2D array of float32 in size of (   l, w).
    /geometry
        /height               2D array of float32 in size of (   l, w) in meter.
        /incidenceAngle       2D array of float32 in size of (   l, w) in degree.
        /slantRangeDistance   2D array of float32 in size of (   l, w) in meter.
        /azimuthAngle         2D array of float32 in size of (   l, w) in degree. (optional)
        /shadowMask           2D array of bool    in size of (   l, w).           (optional)
        /waterMask            2D array of bool    in size of (   l, w).           (optional)
        /bperp                3D array of float32 in size of (n, l, w) in meter.  (optional)
```

### 2. Metadata

Besides [the attributes used in MintPy](./api/attributes.md), we add extra metadata inheritated from [UNAVCO InSAR Product Archive](https://winsar.unavco.org/insar/) ([format specification](https://docs.google.com/document/d/1fm6RY8aL4hhRa88M9cd_Ejh6OL3YfibfjN1UQ7TWsmI/edit?usp=sharing)) as below:

#### 2.1 required & manual

The following metadata requires manual specification in the custom template file, e.g. [WellsEnvD2T399.txt](./templates/WellsEnvD2T399.txt).

+   **mission:** short name of the air-/space-borne SAR (constellation) mission, e.g. ALOS, ALOS2, CSK, ENV, ERS, JERS, NISAR, RS1, RS2, S1, TSX, UAV [auto-grabbed for tops/stripmap/alosStack only]
+   **beam_mode:** short name of the beam mode as used by the space agency, e.g. IW for Sentinel-1, SM for stripmap, SL for spotlight, etc. [auto-grabbed for tops/stripmapStack only]
+   **relative_orbit:** relative orbit (track / path) number [auto-grabbed for tops/stripmapStack only]
+   **first/last_frame:** first and last frame number (same if only one frame) [auto-grabbed for tops/alosStack only]

#### 2.2 recommended & manual

+   **beam_swath:** value used by the space agency, e.g. 1/2/3 for Sentinel-1 IW (default: 0)
+   **processing_dem:** DEM data source used during processing, e.g. SRTM, ASTER, NED (default: Unknown) ...
+   **unwrap_method:** method used for phase unwrapping, e.g. snaphu (default: Unknown)
+   **atmos_correct_method:** method/model used for atmospheric correction (default: None)

#### 2.3 auto-grabbed / hardwired by script

+   **first/last_date:** ISO 8601 format (YYYY-MM-DD) [auto-grabbed]
+   **data_footprint:** WKT formatted polygon outlining the area covered by the data [precise; for geocoded file only; auto-grabbed from Y/X_FIRST/STEP]. This is temporary and should be merged into scene_footprint.
+   **scene_footprint:** WKT formatted polygon outlining the area covered by the data [coarse; auto-grabbed from LON/LAT_REF1/2/3/4]
+   **processing_type:** data product processing level, e.g. INTERFEROGRAM, LOS_VELOCITY, LOS_TIMESERIES [hardwired as LOS_TIMESERIES]
+   **history:** creation date and time in ISO 8601 format (YYYY-MM-DD) [auto-grabbed]
+   **processing_software:** method/software used to generate interferograms/offsets (default: isce) [auto-grabbed from PROCESSOR]
+   **post_processing_software:** method/software used to generate time-series or velocity field [hardwired as MintPy].
+   **flight_direction:** flight direction of the satellite platform, A(scending) or D(escending) (default: Unknown) [auto-grabbed from ORBIT_DIRECTION]
+   **look_direction:** R(ight) or L(eft) [auto-grabbed from ANTENNA_SIDE]
+   **polarization:** transmit and received polarization of the radar wave: HH, VV, HH+VV, etc. (default: Unknown) [auto-grabbed from POLARIZATION]
+   **prf:** pulse repetition frequency (default: 0) [auto-grabbed from PRF]
+   **wavelength:** radar wavelength [auto-grabbed from WAVELENGTH]


### 3. Filename convention ###

Inherited from [UNAVCO InSAR Product Archive](https://winsar.unavco.org/insar/) ([format specification](https://docs.google.com/document/d/1fm6RY8aL4hhRa88M9cd_Ejh6OL3YfibfjN1UQ7TWsmI/edit?usp=sharing)), we use the filename convention below:

&lt;SAT>\_&lt;SW>\_&lt;RELORB>\_&lt;FRAME1>(\_&lt;FRAME2>)\_&lt;DATE1>\_&lt;DATE2>(\_&lt;SUB>).he5

E.g. S1_IW12_128_0593_0597_20141213_20170928.he5

  | Items     | Descriptions | Values |
  | --------- | ------------ | -------|
  | &lt;SAT>    | Mission name | ALOS, ALOS2, CSK, ENV, ERS, JERS, NISAR, RS1, RS2, S1, TSX, UAV |
  | &lt;SW>     | Beam mode with swath number   | SM2 (for ENV), IW3 (for S1) |
  | &lt;RELORB> | Relative orbit (track) number | 3 digits with zero padding  |
  | &lt;FRAME1> | Start frame number | 4 digits with zero padding  |
  | &lt;FRAME2> | End frame number   | 4 digits with zero padding; shown only if it's different from \<FRAME1>  |
  | &lt;DATE1>  | Start date         | YYYYMMDD |
  | &lt;DATE2>  | End date           | YYYYMMDD; "XXXXXXXX" if update mode is ON. |
  | &lt;SUB>    | Subset range       | N{:05d}_S{:05d}_W{:05d}_E{:05d} in degrees; number with precision of 3 digits after decimal * 1000; e.g. S00500_N01300_W001200_E005800; shown only if data is cropped. |

### 4. Web Viewer ###

HDF-EOS5 file format is used as the input of the University of Miami's web viewer for InSAR time-series products. Below is a screenshot of the web viewer for the dataset on Kuju volcano from ALOS-1 ascending track 422.

<p align="center"><b>http://insarmaps.miami.edu</b><br></p>

[![InSAR Web Viewer](https://insarlab.github.io/figs/docs/mintpy/insarmaps-KujuAlosA422.png)](http://insarmaps.miami.edu/)
