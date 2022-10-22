#! /bin/sh

##### Customized Processing After 1st smallbaselineApp
project='KyushuT424F610_640AlosA'
#####
trans='geomap_4rlks.trans'
tmpl=$TE/$project.template


#cd $SC/$project/TSSAR

smallbaselineApp=0
plot_load_data=0
modify_network=0
unwrap_error=0
phaseCor=0
ts_sum=0
ref_dates=0
deramp=0
timeseries_view=0
velocity=0
geo=1
geo2=1


##------------ MintPy routine processing --------------------##
if [ $smallbaselineApp -eq 1 ]; then
    smallbaselineApp.py $tmpl
fi


##------------ Plot Load Data ------------------------------##
if [ $plot_load_data -eq 1 ]; then
    ## Average Spatial Coherence, more reliable mask source than non-nan and temporal coherence
    #generate_mask.py -f average_spatial_coherence.h5 -m 0.7 -o Mask_spatialCoh.h5

    view='view4job.py -c gray -m 0 -M 1 --nodisplay -f'
    #$view average_spatial_coherence.h5
    #$view Mask_spatialCoh.h5
    #$view Mask.h5
    #$view temporal_coherence.h5

    msk='Mask_spatialCoh.h5'
    #remove_ramp.py Seeded_LoadedData.h5 quadratic $msk

    view='view4job.py -r 3 -p 12 --nodisplay'
    #$view -m  0 -M 1 -c gray     -f Coherence.h5
    #$view -m -7 -M 7 --mask $msk -f Seeded_LoadedData.h5
    #$view -m -7 -M 7 --mask $msk -f Seeded_LoadedData_quadratic.h5

    #unwrap_error.py -f Seeded_LoadedData.h5 -m $msk -o Seeded_LoadedData_unwCorTri.h5
    #remove_plane.py Seeded_LoadedData_unwCorTri.h5 quadratic $msk
fi


## Based on plotted Coherence and LoadedData image, select low coherent ifgs and put them in template file.
##------------ Modify Network and Re-run smallbaselineApp ----------##
if [ $modify_network -eq 1  ]; then
    modify_network.py -f LoadedData.h5 -C Coherence.h5 -T $tmpl
    mean_spatial_coherence.py Modified_Coherence.h5 Modified_average_spatial_coherence.h5
#    generate_mask.py -f Modified_average_spatial_coherence.h5 -m 0.75 -o Mask_spatialCoh_Modified.h5

    view='view4job.py -m 0 -M 1 -c gray --nodisplay -f'
    $view Modified_Mask.h5
    $view Modified_average_spatial_coherence.h5
    #plot_network.py -f Modified_LoadedData.h5 --nodisplay
fi



##------------ PU Error correction -------------------------##
##### PU Correction using bonding points
igram='Seeded_Modified_LoadedData'
if [ $unwrap_error -eq 1 ]; then

  ## Generate Individual Mask for Unwrapping Error Correction
  ue_generate_mask=0
  if [ $ue_generate_mask -eq 1 ]; then
    #reference_point.py -f Modified_LoadedData.h5 -t $tmpl
    #cp temporal_coherence.h5 temporal_coherence_1.h5
    #generate_mask.py -f temporal_coherence_1.h5 -m 0.7 -o mask_1.h5

    ##### Calculate with Subset, faster
    ## width = 2832, length = 6728
    sub='-y 4500:6700 -x 0:1500'
    sub2='-y -4500:2228 -x 0:2832'
    subset.py -f ${igram}.h5 $sub
    subset.py -f Mask_spatialCoh.h5 $sub

    reference_point.py       -f subset_${igram}.h5 -x 550 -y 550 -M subset_Mask_spatialCoh.h5 -o subset_${igram}_2.h5
    igram_inversion.py -f subset_${igram}_2.h5 -o subset_timeseries_2.h5
    temporal_coherence.py subset_${igram}_2.h5    subset_timeseries_2.h5 subset_temporal_coherence_3.h5
    subset.py -f subset_temporal_coherence_3.h5 $sub2 -o temporal_coherence_3.h5
    generate_mask.py -f temporal_coherence_3.h5 -m 0.7 -o mask_3.h5
    image_math.py mask_3.h5 '*' 3 mask_3_3.h5


    #reference_point.py       -f ${igram}.h5 -x 840 -y 5390 -o ${igram}_2.h5
    #igram_inversion.py -f ${igram}_2.h5 -o timeseries_2.h5
    #temporal_coherence.py ${igram}_2.h5    timeseries_2.h5 temporal_coherence_2.h5
    #generate_mask.py -f temporal_coherence_2.h5 -m 0.7 -o mask_2.h5
    #image_math.py mask_2.h5 '*' 2 mask_2_2.h5

    reference_point.py       -f ${igram}.h5 -x 550 -y 5050 -o ${igram}_2.h5
    igram_inversion.py -f ${igram}_2.h5 -o timeseries_2.h5
    temporal_coherence.py ${igram}_2.h5    timeseries_2.h5 temporal_coherence_3.h5
    generate_mask.py -f temporal_coherence_3.h5 -m 0.7 -o mask_3.h5
    image_math.py mask_3.h5 '*' 3 mask_3_3.h5

    reference_point.py       -f ${igram}.h5 -x 630 -y 5860 -o ${igram}_2.h5
    igram_inversion.py -f ${igram}_2.h5 -o timeseries_2.h5
    temporal_coherence.py ${igram}_2.h5    timeseries_2.h5 temporal_coherence_4.h5
    generate_mask.py -f temporal_coherence_4.h5 -m 0.7 -o mask_4.h5
    image_math.py mask_4.h5 '*' 4 mask_4_4.h5

  fi

    #cd mask
    #generate_mask.py -f temporal_coherence_2.h5 -m 0.7 -x 660:930 -y 5200:5500 -o mask_2.h5
    #image_math.py mask_2.h5 '*' 2 mask_2_2.h5

    #diff.py mask_1.h5 mask_2.h5 mask_1_2.h5
    #generate_mask.py -f mask_1_2.h5 -m 0.5 -x 910:944 -y 5457:5496 -o mask_small_1.h5
    #generate_mask.py -f mask_1_2.h5 -m 0.5 -x 877:985 -y 5263:5357 -o mask_small_2.h5
    #diff.py mask_1_2.h5 mask_small_1.h5 mask_1_3.h5
    #diff.py mask_1_3.h5 mask_small_2.h5 mask_1_4.h5

    #generate_mask.py -f mask_1_4.h5 -m 0.5 -o mask_1_1.h5
    #cd ..

    #add.py -f mask_1_1.h5,mask_2_2.h5,mask_3_3.h5,mask_4_4.h5 -o Mask_all.h5
    msk='Mask_all.h5'
    view='view.py --nodisplay -f'
    #$view mask_1.h5
    #$view mask_2.h5
    #$view $msk

    #reference_point.py -f Modified_LoadedData.h5 -t $tmpl -M $msk

    unwrap_error.py -f ${igram}.h5 -m $msk -t $tmpl --ramp quadratic --save-detrend
    remove_plane.py ${igram}_unwCor.h5 quadratic $msk

    view='view4job.py -r 6 -p 19 -m -7 -M 7 --mask Mask_all.h5 --noaxis --nodisplay -f'
    #$view ${igram}.h5
    $view ${igram}_unwCor.h5
    $view ${igram}_unwCor_quadratic.h5

    #### Mask with Bonding Points
    #point_yx='5430,955,5429,865,5293,834,5178,832,5403,780,5714,612'
    #line_yx='5430,955,5429,865;5293,834,5178,832;5403,780,5714,612'
    #view.py --point=${point_yx} --line=${line_yx} --nodisplay -f $msk -o bonding_points.png

    #smallbaselineApp.py $tmpl

fi


##### Display Iono Pairs
view='view.py -r 3 -p 5 -m -7 -M 7 -f Seeded_Modified_LoadedData_unwCor_plane_masked.h5 --nodisplay'
#$view -d 100606 -o Seeded_Modified_LoadedData_unwCor_plane_masked_100606.png
#$view -d 060526 -o Seeded_Modified_LoadedData_unwCor_plane_masked_060526.png

##------------ Mask ----------------------------------------##
#igram_inversion.py    Seeded_Modified_LoadedData_unwCor.h5
#temporal_coherence.py Seeded_Modified_LoadedData_unwCor.h5 timeseries.h5
#generate_mask.py -f temporal_coherence.h5 -m 0.7 -x 180:2560 -o Mask_tempCoh_dis.h5
msk='Mask_tempCoh_dis.h5'
view='view.py -m 0 -M 1 -c gray --nodisplay -f'
#$view $msk
#$view temporal_coherence.h5


##------------ Phase Correction ----------------------------##
if [ $phaseCor -eq 1 ]; then
    remove_plane.py -f timeseries_ECMWF.h5 -t $tmpl
    dem_error.py timeseries_ECMWF_plane.h5 Seeded_Modified_LoadedData_unwCor.h5
    masking.py -f DEM_error.h5 -m $msk
    view.py -m -25 -M 25 --nodisplay -f DEM_error_masked.h5
fi
ts='timeseries_ECMWF_plane_demCor.h5'



##------------ Reference & Exclude Dates -------------------##
## Reference date and drop turbulence dates
if [ $ts_sum -eq 1 ]; then
    #sum_epochs.py timeseries_ECMWF_demCor_plane.h5 ts_sum_0.h5
    view.py -r 3  -p 10 -m 0 -M 1 --mask $msk --noaxis --nodisplay -f ts_sum_0.h5
fi

#ex_date='20060526,20060826,20070111,20091019,20100606,20100906'
ex_date='20060624,20070812,20080629,20090702,20090817,20091002,20100705'
ref_date='20090214'

if [ $ref_dates -eq 1 ]; then
    #reference_date.py timeseries_ECMWF_demCor.h5 $ref_date
    #remove_plane.py -f timeseries_ECMWF_demCor_ref${ref_date}.h5 -s quadratic -m $msk
    ts='timeseries_ECMWF_demCor_ref'${ref_date}'_quadratic.h5'
    #sum_epochs.py $ts ts_sum.h5
    view4job.py -r 3 -p 10 -m 0 -M 1 --noaxis --nodisplay --mask $msk -f ts_sum.h5

fi

if [ $deramp -eq 1 ]; then
    ysub='0,4000,3800,6728'
    remove_plane.py -f timeseries_ECMWF_demCor_ref${ref_date}.h5 -s quadratic -m $msk -y $ysub --save-mask
    ts='timeseries_ECMWF_demCor_ref'${ref_date}'_quadratic.h5'

    #sum_epochs.py $ts ts_sum.h5
    #view4job.py -r 2 -p 11 -m 0 -M 1 --noaxis --nodisplay --mask $msk -f ts_sum.h5
    #diff.py timeseries_ECMWF_demCor_ref${ref_date}.h5 ${ts} 2quadratic.h5

    view='view.py -r 4 -p 8 --noaxis --nodisplay --mask '$msk
    #$view       -m  0  -M 1  -f ts_sum.h5
    #$view -u cm -m -25 -M 25 -f quadratic.h5
    #$view -u cm -m -5  -M 5  -f timeseries_ECMWF_demCor_quadratic.h5

    #reference_date.py timeseries_ECMWF_demCor_quadratic.h5 $ref_date
    #$view -u cm -m -5 -M 5 -f timeseries_ECMWF_demCor_quadratic_ref${ref_date}.h5
    #$view -u cm -m -5 -M 5 -f timeseries_ECMWF_demCor_quadratic_ref${ref_date}.h5 -E ${ex_date}
fi


ts='timeseries_ECMWF_demCor_ref'${ref_date}'_quadratic.h5'
##------------ Time Series Correction View -----------------##
## MERRA
#tropo_pyaps.py -f timeseries.h5 -d radar_4rlks.hgt -s MERRA -h 12:00 -i incidence_angle.h5
if [ $timeseries_view -eq 1 ]; then
    view='view4job.py -r 3 -p 10 -u cm -m -25 -M 25 --mask '$msk' --noaxis --nodisplay -f'
    #$view timeseries.h5
    #$view timeseries_ECMWF.h5
    #$view timeseries_ECMWF_demCor.h5
    #$view timeseries_ECMWF_demCor_ref${ref_date}.h5
    #$view ${ts} -o ${ts}_0.png

    view='view4job.py -r 3 -p 10 -u cm -m  -5 -M  5 --mask '$msk' --noaxis --nodisplay -f'
    #$view ${ts}
    $view ${ts} -E $ex_date

    ##### error
    #diff.py timeseries_ECMWF_demCor_plane.h5 timeseries_ECMWF_demCor.h5 plane.h5
    #masking.py -m $msk -f ECMWF.h5
    #masking.py -m $msk -f MERRA.h5
    #masking.py -m $msk -f DEM_error.h5
    #masking.py -m $msk -f plane.h5
    #view.py -r 4 -p 8 -u cm  -m -5  -M 5   --noaxis --nodisplay -f ECMWF_masked.h5
    #view.py -r 4 -p 8 -c jet -m -25 -M 25           --nodisplay -f DEM_error_masked.h5
    #view.py -r 4 -p 8 -u cm  -m -25 -M 25  --noaxis --nodisplay -f plane_masked.h5
fi


##------------ Velocity ---------------------------------------##
if [ $velocity -eq 1 ]; then
    timeseries2velocity.py -f $ts -E $ex_date
    view.py -u cm/yr -m -2 -M 2 --mask $msk --nodisplay -f velocity_ex.h5
fi

if [ $geo -eq 1 ]; then
    ##### Geocoding
    #geocode.py $msk $trans
    #geocode.py velocity_ex.h5 $trans
    #geocode.py ${ts} $trans
    geocode.py Modified_Coherence.h5 $trans

    ## Find subset.lalo based on geo_velocity.h5

    ##### Seeding
    #reference_point.py -t $tmpl -f geo_velocity_ex.h5
    #reference_point.py -t $tmpl -f geo_${ts}
fi

if [ $geo2 -eq 1 ]; then
    ##### Subset
    #subset.py -t $tmpl -f Seeded_geo_*.h5,geo_$msk,gsi10m_30m.dem --parallel
    subset.py -t $tmpl -f geo_Modified_Coherence.h5

    ##### masking
    #masking.py -m subset_geo_$msk -f subset_Seeded_geo_velocity_ex.h5
    #masking.py -m subset_geo_$msk -f subset_Seeded_geo_${ts}

    #multi_looking.py subset_gsi10m_30m.dem 5 5 subset_gsi10m_150m.dem
    #view.py -u cm/yr -m -2 -M 2 --mask subset_geo_$msk -D subset_gsi10m_150m.dem --save -f subset_Seeded_geo_velocity_ex.h5
    #save_kml.py -m -0.02 -M 0.02 -f subset_Seeded_geo_velocity_ex_masked.h5
fi
