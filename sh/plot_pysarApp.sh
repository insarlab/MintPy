#! /bin/sh
###############################################################
# Plot Results from Routine Workflow with pysarApp.py
# Author: Zhang Yunjun, 2017-07-23
# Latest update: 2019-04-01
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
dem_file='./INPUTS/geometryRadar.h5'
if [ ! -f $dem_file ]; then
    dem_file='./INPUTS/geometryGeo.h5'
fi


## Log File
log_file='plot_pysarApp.log'
touch $log_file
echo "\n\n\n\n\n" >> $log_file
echo "#############################  ./plot_pysarApp.sh  ###########################" >> $log_file
date >> $log_file
echo "##############################################################################" >> $log_file
#use "echo 'yoyoyo' | tee -a log" to output message to both screen and file.
#use "echo 'yoyoyo' >> log" to output message to file only.

## Create PIC folder
if [ ! -d "PIC" ]; then
    echo 'Create ./PIC folder'
    mkdir PIC
fi

## common view.py option for all files
view='view.py --nodisplay --dpi 150 --update '

## Plot Key files
opt=' --dem '$dem_file' --mask '$mask_file' -u cm '
#opt=' --dem '$dem_file' --mask '$mask_file' -u cm --vlim -2 2'
if [ $plot_key_files -eq 1 ]; then
    file=velocity.h5;              test -f $file && $view $file $opt               >> $log_file
    file=temporalCoherence.h5;     test -f $file && $view $file -c gray --vlim 0 1 >> $log_file
    file=maskTempCoh.h5;           test -f $file && $view $file -c gray --vlim 0 1 >> $log_file
    file=INPUTS/geometryRadar.h5;  test -f $file && $view $file                    >> $log_file
    file=INPUTS/geometryGeo.h5;    test -f $file && $view $file                    >> $log_file
fi


## Loaded Dataset
if [ $plot_loaded_data -eq 1 ]; then
    file=INPUTS/ifgramStack.h5
    test -f $file && $view $file unwrapPhase-  --zero-mask --wrap >> $log_file
    test -f $file && $view $file unwrapPhase-  --zero-mask        >> $log_file
    test -f $file && $view $file coherence-    --mask no          >> $log_file
fi


## Auxliary Files from loaded dataset
if [ $plot_loaded_data_aux -eq 1 ]; then
    file=avgPhaseVelocity.h5;   test -f $file && $view $file                      >> $log_file
    file=avgSpatialCoh.h5;      test -f $file && $view $file -c gray --vlim 0 1   >> $log_file
    file=maskConnComp.h5;       test -f $file && $view $file -c gray --vlim 0 1   >> $log_file
fi


## Time-series files
opt='--mask '$mask_file' --noaxis -u cm '
#opt='--mask '$mask_file' --noaxis -u cm --vlim -10 10 '
if [ $plot_timeseries -eq 1 ]; then
    file=timeseries.h5;                             test -f $file && $view $file $opt >> $log_file

    #LOD for Envisat
    file=timeseries_LODcor.h5;                      test -f $file && $view $file $opt >> $log_file
    file=timeseries_LODcor_ECMWF.h5;                test -f $file && $view $file $opt >> $log_file
    file=timeseries_LODcor_ECMWF_demErr.h5;         test -f $file && $view $file $opt >> $log_file
    file=timeseries_LODcor_ECMWF_ramp.h5;           test -f $file && $view $file $opt >> $log_file
    file=timeseries_LODcor_ECMWF_ramp_demErr.h5;    test -f $file && $view $file $opt >> $log_file

    #w trop delay corrections
    for trop in '_ECMWF' '_MERRA' '_NARR' '_tropHgt'
    do
        file=timeseries${trop}.h5;                  test -f $file && $view $file $opt >> $log_file
        file=timeseries${trop}_demErr.h5;           test -f $file && $view $file $opt >> $log_file
        file=timeseries${trop}_ramp.h5;             test -f $file && $view $file $opt >> $log_file
        file=timeseries${trop}_ramp_demErr.h5;      test -f $file && $view $file $opt >> $log_file
    done

    #w/o trop delay correction
    file=timeseries_ramp.h5;                        test -f $file && $view $file $opt >> $log_file
    file=timeseries_demErr_ramp.h5;                 test -f $file && $view $file $opt >> $log_file
fi


## Geo coordinates for UNAVCO Time-series InSAR Archive Product
if [ $plot_geocoded_data -eq 1 ]; then
    file=./GEOCODE/geo_maskTempCoh.h5;                   test -f $file && $view $file -c gray  >> $log_file
    file=./GEOCODE/geo_temporalCoherence.h5;             test -f $file && $view $file -c gray  >> $log_file
    file=./GEOCODE/geo_velocity.h5;                      test -f $file && $view $file velocity >> $log_file
    file=./GEOCODE/geo_timeseries_ECMWF_demErr_ramp.h5;  test -f $file && $view $file --noaxis >> $log_file
    file=./GEOCODE/geo_timeseries_ECMWF_demErr.h5;       test -f $file && $view $file --noaxis >> $log_file
    file=./GEOCODE/geo_timeseries_demErr_ramp.h5;        test -f $file && $view $file --noaxis >> $log_file
    file=./GEOCODE/geo_timeseries_demErr.h5;             test -f $file && $view $file --noaxis >> $log_file
fi


if [ $plot_the_rest -eq 1 ]; then
    for trop in 'Ecmwf' 'Merra' 'Narr'
    do
        file=velocity${trop}.h5;    test -f $file && $view $file --mask no >> $log_file
    done
    file=numInvIfgram.h5;           test -f $file && $view $file --mask no >> $log_file
fi


## Move picture files to PIC folder
echo "Move *.png/pdf/kmz files into ./PIC folder."
mv *.png PIC/
mv *.pdf PIC/
mv *.kmz ./GEOCODE/*.kmz PIC/

