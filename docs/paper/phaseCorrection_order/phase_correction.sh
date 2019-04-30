#! /bin/sh
## Truth
timeseries2velocity.py timeseriesSim.h5 -o velocity_sim.h5

## Phase correction used in PYSAR
diff.py timeseries.h5 ECMWF.h5 -o timeseries_ECMWF.h5
dem_error.py timeseries_ECMWF.h5 -o timeseries_ECMWF_demErr.h5
mv demErr.h5 demErr_pysar.h5
timeseries_rms.py timeseriesResidual.h5 -m waterMask.h5
reference_date.py timeseries_ECMWF_demErr.h5 --ref-date reference_date.txt
remove_ramp.py timeseries_ECMWF_demErr_refDate.h5 -m waterMask.h5
timeseries2velocity.py timeseries_ECMWF_demErr_refDate_plane.h5 -o velocity_pysar.h5

## Phase correction used in GIANT
diff.py timeseries.h5 ECMWF.h5 -o timeseries_ECMWF.h5
remove_ramp.py timeseries_ECMWF.h5 -m waterMask.h5
dem_error.py timeseries_ECMWF_plane.h5 #-g geometryRadar.h5
mv demErr.h5 demErr_giant.h5
timeseries2velocity.py timeseries_ECMWF_plane_demErr.h5 -o velocity_giant.h5

diff.py timeseries_ECMWF_demErr_refDate_plane.h5 timeseries_ECMWF_plane_demErr.h5 -o timeseries_diff.h5
