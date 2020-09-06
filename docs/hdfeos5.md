We support output geocoded displacement time-series product into [HDF-EOS5](http://hdfeos.org) format via `save_hdfeos5.py`. This is designed to easily share the InSAR time-series products to the broader community.

```bash
save_hdfeos5.py geo_timeseries_ECMWF_ramp_demErr.h5 -c geo_temporalCoherence.h5 -m geo_maskTempCoh.h5 -g geo_geometryRadar.h5
```

### File structure ###

```
/                             Root level group
Attributes                    metadata in dict.
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

### Filename convention ###

We use a filename convention inherited from [UNAVCO InSAR Archive format](https://winsar.unavco.org/insar/).

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

### Web Viewer ###

HDF-EOS5 file format is used as the input of the University of Miami's web viewer for InSAR time-series products. Below is a screenshot of the web viewer for the dataset on Kuju volcano from ALOS-1 acending track 422.

<p align="center"><b>http://insarmaps.miami.edu</b><br></p>

[![InSAR Web Viewer](https://yunjunzhang.files.wordpress.com/2019/06/web_viewer_kujualosat422.png)](http://insarmaps.miami.edu/)


