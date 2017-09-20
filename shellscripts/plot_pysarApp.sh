#! /bin/sh
###############################################################
# Plot pysarApp routine processing results.
# Author: Zhang Yunjun, 2017-07-23
###############################################################


## Change to 0 if you do not want to re-plot loaded dataset again
plot_key_files=1
plot_loaded_data=1
plot_loaded_data_aux=1
plot_timeseries=1
plot_geocoded_data=1
plot_the_rest=1


# Default file name
mask_file='maskTempCoh.h5'
geo_mask_file='geo_maskTempCoh.h5'

## Log File
log_file='plot_pysarApp.log'
echo "touch log file: "$log_file
touch $log_file

## Create PIC folder
if [ ! -d "PIC" ]; then
    echo 'Create PIC folder'
    mkdir PIC
fi


## Plot Key files
opt=' -d demRadar.h5 --mask '$mask_file' -u cm '
#opt=' -d demRadar.h5 --mask '$mask_file' -u cm -m -2 -M 2 '
if [ $plot_key_files -eq 1 ]; then
    view.py --nodisplay velocity.h5           $opt               | tee -a $log_file
    view.py --nodisplay temporalCoherence.h5  -c gray -m 0 -M 1  | tee -a $log_file
    view.py --nodisplay maskTempCoh.h5        -c gray -m 0 -M 1  | tee -a $log_file
    view.py --nodisplay demRadar_error.h5 --mask $mask_file      | tee -a $log_file
fi


## Loaded Dataset
if [ $plot_loaded_data -eq 1 ]; then
    view.py --nodisplay unwrapIfgram.h5 --mask $mask_file  | tee -a $log_file
    view.py --nodisplay coherence.h5 -c gray -m 0 -M 1     | tee -a $log_file
    view.py --nodisplay demRadar.h5                        | tee -a $log_file
    view.py --nodisplay demGeo.h5 --lalo-label             | tee -a $log_file
fi


## Auxliary Files from loaded dataset
if [ $plot_loaded_data_aux -eq 1 ]; then
    view.py --nodisplay unwrapIfgram_stack.h5                        | tee -a $log_file
    view.py --nodisplay mask.h5                    -c gray -m 0 -M 1 | tee -a $log_file
    view.py --nodisplay averageSpatialCoherence.h5 -c gray -m 0 -M 1 | tee -a $log_file
fi


## Time-series files
view='view.py --nodisplay --mask '$mask_file' --noaxis -u cm '
#view='view.py --nodisplay --mask '$mask_file' --noaxis -u cm -m -10 -M 10 '
if [ $plot_timeseries -eq 1 ]; then
    $view timeseries.h5        | tee -a $log_file

    $view timeseries_LODcor_ECMWF.h5                           | tee -a $log_file
    $view timeseries_LODcor_ECMWF_demErr.h5                    | tee -a $log_file
    $view timeseries_LODcor_ECMWF_demErr_refDate.h5            | tee -a $log_file
    $view timeseries_LODcor_ECMWF_demErr_refDate_plane.h5      | tee -a $log_file
    $view timeseries_LODcor_ECMWF_demErr_refDate_quadratic.h5  | tee -a $log_file

    $view timeseries_ECMWF.h5                           | tee -a $log_file
    $view timeseries_ECMWF_demErr.h5                    | tee -a $log_file
    $view timeseries_ECMWF_demErr_refDate.h5            | tee -a $log_file
    $view timeseries_ECMWF_demErr_refDate_plane.h5      | tee -a $log_file
    $view timeseries_ECMWF_demErr_refDate_quadratic.h5  | tee -a $log_file

    $view timeseries_demErr.h5                    | tee -a $log_file
    $view timeseries_demErr_refDate.h5            | tee -a $log_file
    $view timeseries_demErr_refDate_plane.h5      | tee -a $log_file
    $view timeseries_demErr_refDate_quadratic.h5  | tee -a $log_file

    $view timeseries_demErr.h5                            | tee -a $log_file
    $view timeseries_demErr_tropHgt.h5                    | tee -a $log_file
    $view timeseries_demErr_tropHgt_refDate.h5            | tee -a $log_file
    $view timeseries_demErr_tropHgt_refDate_plane.h5      | tee -a $log_file
    $view timeseries_demErr_tropHgt_refDate_quadratic.h5  | tee -a $log_file
fi


## Geo coordinates for UNAVCO Time-series InSAR Archive Product
if [ $plot_geocoded_data -eq 1 ]; then
    view.py --nodisplay --lalo-label geo_incidenceAngle.h5                                    | tee -a $log_file
    view.py --nodisplay --lalo-label geo_maskTempCoh.h5          -c gray -m 0 -M 1            | tee -a $log_file
    view.py --nodisplay --lalo-label geo_temporalCoherence.h5    -c gray -m 0 -M 1            | tee -a $log_file
    view.py --nodisplay --lalo-label geo_velocity.h5             --mask $geo_mask_file -u cm  | tee -a $log_file
    view.py --nodisplay --lalo-label geo_timeseries_LODcor*.h5 --noaxis --mask $geo_mask_file -u cm  | tee -a $log_file
    view.py --nodisplay --lalo-label geo_timeseries_ECMWF*.h5  --noaxis --mask $geo_mask_file -u cm  | tee -a $log_file
    view.py --nodisplay --lalo-label geo_timeseries_demErr*.h5 --noaxis --mask $geo_mask_file -u cm  | tee -a $log_file
fi


if [ $plot_the_rest -eq 1 ]; then
    view.py --nodisplay velocityStd.h5    --mask $mask_file -u cm  | tee -a $log_file
    view.py --nodisplay velocityEcmwf.h5  --mask $mask_file -u cm  | tee -a $log_file
fi


## Move picture files to PIC folder
echo "Move *.png *.pdf into PIC folder"
mv *.png PIC/
mv *.pdf PIC/
mv *.kmz PIC/
mv $log_file PIC/

