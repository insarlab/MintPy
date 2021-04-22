#! /bin/sh
###############################################################
# Plot Results from Routine Workflow with smallbaselineApp.py
# Author: Zhang Yunjun, 2017-07-23
# Latest update: 2021-03-08
###############################################################
# Update the date above to enable auto copyover/overwrite


## Change to 0 if you do not want to re-plot loaded dataset again
plot_key_files=1
plot_loaded_data=1
plot_loaded_data_aux=1
plot_timeseries=1
plot_geocoded_data=1
plot_the_rest=1


# Default file name
mask_file='maskTempCoh.h5'
dem_file='./inputs/geometryRadar.h5'
if [ ! -f $dem_file ]; then
    dem_file='./inputs/geometryGeo.h5'
fi

## Log File
log_file='plot_smallbaselineApp.log'
touch $log_file
printf "\n\n\n\n\n" >> $log_file
echo "########################  ./plot_smallbaselineApp.sh  ########################" >> $log_file
date >> $log_file
echo "##############################################################################" >> $log_file
#use "echo 'yoyoyo' | tee -a log" to output message to both screen and file.
#use "echo 'yoyoyo' >> log" to output message to file only.

## Create pic folder
if [ ! -d "pic" ]; then
    echo 'Create ./pic folder'
    mkdir pic
fi

## common view.py option for all files
view='view.py --nodisplay --dpi 150 --update '

## Plot Key files
opt=' --dem '$dem_file' --mask '$mask_file' -u cm '
#opt=' --dem '$dem_file' --mask '$mask_file' -u cm --vlim -2 2'
if [ $plot_key_files -eq 1 ]; then
    file='velocity.h5';              test -f $file && $view $file $opt               >> $log_file
    file='temporalCoherence.h5';     test -f $file && $view $file -c gray --vlim 0 1 >> $log_file
    file='maskTempCoh.h5';           test -f $file && $view $file -c gray --vlim 0 1 >> $log_file
    file='inputs/geometryRadar.h5';  test -f $file && $view $file                    >> $log_file
    file='inputs/geometryGeo.h5';    test -f $file && $view $file                    >> $log_file
fi


## Loaded Dataset
if [ $plot_loaded_data -eq 1 ]; then
    file='inputs/ifgramStack.h5'
    test -f $file && h5ls $file/unwrapPhase      && $view $file unwrapPhase- --zero-mask --wrap -c cmy >> $log_file
    test -f $file && h5ls $file/unwrapPhase      && $view $file unwrapPhase- --zero-mask               >> $log_file
    test -f $file && h5ls $file/coherence        && $view $file coherence-         --mask no -v 0 1    >> $log_file
    test -f $file && h5ls $file/connectComponent && $view $file connectComponent-  --mask no           >> $log_file

    # phase-unwrapping error correction
    for dset in 'unwrapPhase_bridging' 'unwrapPhase_phaseClosure' 'unwrapPhase_bridging_phaseClosure'; do
        test -f $file && h5ls $file/$dset        && $view $file $dset-             --zero-mask        >> $log_file
    done
fi


## Auxliary Files from loaded dataset
if [ $plot_loaded_data_aux -eq 1 ]; then
    file='avgPhaseVelocity.h5';   test -f $file && $view $file                      >> $log_file
    file='avgSpatialCoh.h5';      test -f $file && $view $file -c gray --vlim 0 1   >> $log_file
    file='maskConnComp.h5';       test -f $file && $view $file -c gray --vlim 0 1   >> $log_file
fi


## Time-series files
opt='--mask '$mask_file' --noaxis -u cm --wrap --wrap-range -10 10 '
if [ $plot_timeseries -eq 1 ]; then
    file='timeseries.h5'; test -f $file && $view $file $opt >> $log_file
    find . -name 'timeseries_*.h5' -exec   $view {}    $opt >> $log_file \;
fi


## Geo coordinates for UNAVCO Time-series InSAR Archive Product
if [ $plot_geocoded_data -eq 1 ]; then
    file='./geo/geo_maskTempCoh.h5';          test -f $file && $view $file -c gray  >> $log_file
    file='./geo/geo_temporalCoherence.h5';    test -f $file && $view $file -c gray  >> $log_file
    file='./geo/geo_velocity.h5';             test -f $file && $view $file velocity >> $log_file
    find . -name './geo/geo_timeseries_*.h5' -exec             $view {}    $opt     >> $log_file \;
fi


if [ $plot_the_rest -eq 1 ]; then
    for tropo in 'ERA5' 'ERAI' 'ECMWF' 'MERRA' 'NARR'; do
        file='velocity'${tropo}'.h5';  test -f $file && $view $file --mask no >> $log_file
    done
    file='numInvIfgram.h5';            test -f $file && $view $file --mask no >> $log_file
fi


## Move/copy picture files to pic folder
echo "Copy *.txt files into ./pic folder."
cp *.txt pic/
echo "Move *.png/pdf/kmz files into ./pic folder."
mv *.png *.pdf *.kmz ./geo/*.kmz pic/
