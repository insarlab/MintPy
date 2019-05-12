#! /bin/sh
## Compare the estimated velocities from time-series with different troposphric delay corrections
## Created by Zhang Yunjun, Sep 6th, 2018

## Run Control: change to 0 to skip the part
run_proc=1
run_plot=1

## Run Customized Processing
if [ $run_proc -eq 1 ]; then
    generate_mask.py temporalCoherence.h5 -m 0.8 -o maskTempCoh.h5

    # Tropospheric Correction with ERA-Interim
    tropo_pyaps.py -f timeseries.h5 -g inputs/geometryRadar.h5 -m ECMWF -w ~/insarlab/WEATHER
    dem_error.py  timeseries_ECMWF.h5 -g inputs/geometryRadar.h5
    remove_ramp.py timeseries_ECMWF_demErr.h5 -m maskTempCoh.h5 -s linear
    timeseries2velocity.py timeseries_ECMWF_demErr_ramp.h5 -o velocity_tropECMWF.h5

    # Tropospheric Correction with MERRA-2
    tropo_pyaps.py -f timeseries.h5 -g inputs/geometryRadar.h5 -m MERRA -w ~/insarlab/WEATHER
    dem_error.py  timeseries_MERRA.h5 -g inputs/geometryRadar.h5
    remove_ramp.py timeseries_MERRA_demErr.h5 -m maskTempCoh.h5 -s linear
    timeseries2velocity.py timeseries_MERRA_demErr_ramp.h5 -o velocity_tropMERRA.h5

    # Tropospheric Correction with Phase/Elevation Ratio
    dem_error.py timeseries.h5 -g -g inputs/geomtropo_phase_elevationetryRadar.h5
    tropo_phase_elevation.py timeseries_demErr.h5 -g inputs/geometryRadar.h5 -m maskTempCoh.h5
    remove_ramp.py timeseries_demErr_tropHgt.h5 -m maskTempCoh.h5 -s linear
    timeseries2velocity.py timeseries_demErr_tropHgt_ramp.h5 -o velocity_tropHgt.h5

    # No tropospheric Correction
    remove_ramp.py timeseries_demErr.h5 -m maskTempCoh.h5 -s linear
    timeseries2velocity.py timeseries_demErr_ramp.h5 -o velocity_tropNone.h5
fi

## Run plotting
opt=' velocity --wrap --wrap-range -1 1 --nodisplay --cbar-nbins 2 --notitle --figsize 3 3 --fontsize 12 --notick --dpi 600 --ref-size 3 --dem inputs/geometryRadar.h5 --dem-nocontour '
proj_name='AlcedoSenDT128'
if [ $run_plot -eq 1 ]; then
    view.py velocity_tropNone.h5  $opt -o ${proj_name}_velocity_0_tropNone.png
    view.py velocity_tropECMWF.h5 $opt -o ${proj_name}_velocity_1_tropECMWF.png
    view.py velocity_tropMERRA.h5 $opt -o ${proj_name}_velocity_2_tropMERRA.png
    view.py velocity_tropHgt.h5   $opt -o ${proj_name}_velocity_3_tropHgt.png
fi
