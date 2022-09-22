#! /bin/sh
# This is an example for creating your own processing recipe using MintPy
# Created by Zhang Yunjun, 2016

##### Initial Setting #####
satellite='AlosA'
track='T424'
frame='F610_640'

sub='-l 31.14:31.30 -L 130.48:130.68'
seed='-l 31.2352 -L 130.5029'
row_col='-r 5 -p 6'

##### Process Flag ########
subset=0
timeseries=0
velocity=0
point_ts=1
key_ifg=0


##------------- Data Source ---------------------------------------##
projectDir=$KYUSHU_PROJECT/${satellite}${track}${frame}/TSSAR
#projectDir=$SC/Kyushu${track}${frame}${satellite}/TSSAR/Product
workDir=$KYUSHU_PROJECT/Volcanoes/Ata/${satellite}${track}

msk=$projectDir'/subset_geo_Mask_tempCoh_dis.h5'
ts=$projectDir'/subset_Seeded_geo_timeseries_ECMWF_demCor_refDate_quadratic.h5'
vel=$projectDir'/subset_Seeded_geo_velocity_ex_masked.h5'
dem=$projectDir'/subset_gsi10m_30m.dem'
coh=$projectDir'/subset_geo_temporal_coherence.h5'

cd $workDir

##--------------- Subset -------------------------------------------##
if [ $subset -eq 1 ]; then

    subset.py --outfill-nan  $sub -o ts.h5          -f ${ts}
    subset.py --outfill-zero $sub -o gsi10m_30m.dem -f ${dem}
    subset.py --outfill-zero $sub -o Mask.h5        -f ${msk}
    subset.py --outfill-zero $sub -o tempCoh.h5     -f ${coh}
    generate_mask.py -f tempCoh.h5 -m 0.7 -o Mask_tempCoh.h5

    view.py -c gray --nodisplay -f Mask.h5
    view.py -c gray --nodisplay -f tempCoh.h5
    view.py -c gray --nodisplay -f Mask_tempCoh.h5
    view.py         --nodisplay -f gsi10m_30m.dem

    msk='Mask_tempCoh.h5'
    view='view.py '${row_col}' -u cm -m -5 -M 5 -D gsi10m_30m.dem --dem-noshade --mask '$msk' --noaxis --nodisplay'
    $view -f ts.h5

    #generate_mask.py -f velocity_ex_masked.h5 -m 0.01 -o mask_vel_1.h5
    #diff.py Mask_tempCoh.h5 mask_vel_1.h5 Mask_tempCoh_dis.h5
fi
msk='Mask_tempCoh_dis.h5'
dem='gsi10m_30m.dem'


##--------------- Timeseries ----------------------------------------##
if [ $timeseries -eq 1 ]; then

    view='view.py '${row_col}' -u cm -m -5 -M 5 -D gsi10m_30m.dem --mask '$msk' --noaxis --nodisplay'

    reference_point.py $seed -f ts.h5 -m $msk
    $view -f Seeded_ts.h5

    sum_epochs.py Seeded_ts.h5
    $view -u 1 -m 0 -M 1 -f sum_Seeded_ts.h5 -c gray

    #mean_spatial.py -f sum_Seeded_ts.h5 -m $msk
    cir_par='31.1800,130.5290,60;31.1892,130.6147,40'
    mean_spatial.py -f sum_Seeded_ts.h5 -m $msk --circle $cir_par

    reference_date.py Seeded_ts.h5 reference_date.txt
    masking.py -f Seeded_ts_refDate.h5 -m $msk
    $view -f Seeded_ts_refDate.h5
    $view -f Seeded_ts_refDate.h5 -E drop_date.txt

    #$view -f Seeded_ts_refDate.h5 --ref-epoch 20110420
fi

ts='Seeded_ts_refDate.h5'
##--------------- Velocity -------------------------------------------##
if [ $velocity -eq 1 ]; then
    timeseries2velocity.py -f $ts -E drop_date.txt
    view='view.py -u cm -m -2 -M 2 -D '$dem' --nodisplay --mask '${msk}' -f '
    $view velocity_ex.h5
    masking.py -f velocity_ex.h5 -m $msk
    save_kmz.py -f velocity_ex_masked.h5 -m -0.02 -M 0.02
fi


##---------------- Point Time Series Displacement --------------------##
if [ $point_ts -eq 1 ]; then

    view='view.py --displacement -u cm -m -5 -M 5 --mask '$msk' -D '$dem' --nodisplay -f'
    ##--- 1.1 Pre-eruptive Inflation
    #save_unw.py Seeded_ts.h5 20100102 20101120
    #$view 100102_101120.unw
    #masking.py -f 100102_101120.unw -m Mask.h5
    #save_kmz.py -f 100102_101120_masked.unw -m -0.05 -M 0.05 --displacement --cbar-label 'LOS displacement'

    ##--- 2. Point TS
    #masking.py -f Seeded_ts_refDate.h5 -m $msk
    ts='Seeded_ts_refDate_masked.h5'
    vel='velocity_ex_masked.h5'
    tsview='tsview.py -f '${ts}' -E drop_date.txt -F '$ts' -r 3 -D '$dem' -v '$vel' -a -0.015 -b 0.015 --rect-color crimson'

    #$tsview -l -10 -h 4 --lalo 31.1783,130.5364 --LALO 31.1853,130.5231      --nodisplay -t 20060623 -T 20110421
    #$tsview -l -7  -h 4 --lalo 31.1923,130.6184 --LALO 31.1959,130.5971 -r 2 --nodisplay -t 20060623 -T 20110421
    $tsview -l -7  -h 4 --lalo 31.1923,130.6184 -r 2 --nodisplay -t 20060623 -T 20110421

    point_lalo='31.1923,130.6184'
    view='view.py -u cm -m -2 -M 2 --mask '$msk' -D '$dem' --nodisplay'
    sub='-l 31.155:31.21 -L 130.505:130.635'
    lineFile='../transect_lonlat.xy'
    $view -f velocity_ex.h5 --point-lalo $point_lalo -o velocity_ex_poi.png $sub --line-lalo $lineFile
fi


##----------------- Velocity and Interferograms ----------------------##
if [ $key_ifg -eq 1 ]; then

    view='view.py --mask '$msk' --displacement -u cm -m -5 -M 5 -D '$dem' --nodisplay'
    save_kmz='save_kmz.py -m -0.05 -M 0.05 --displacement'

    #### Post-eruptive
    d1='110305'
    d2='110420'
    save_unw.py Seeded_ts.h5 $d1 $d2
    #$view -f 091130_110118.unw
    masking.py -f ${d1}_${d2}.unw -m $msk
    $save_kmz -f ${d1}_${d2}_masked.unw --cbar-label 'LOS displacement'

    ts='Seeded_ts_refDate_masked.h5'
    vel='091130_110118_masked.unw'
    tsview='tsview.py -f '${ts}' -E drop_date.txt -F '$ts' -r 3 -D '$dem' -v '$vel' --displacement -a -0.05 -b 0.05 --rect-color black'
    #$tsview --lalo 31.9076,130.8284 -l -3 -h 3 --nodisplay --zoom-y 200:900 --zoom-x 700:1600
fi
