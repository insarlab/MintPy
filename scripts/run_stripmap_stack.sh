#! /bin/sh

# clean folders before re-run
if [ -d "Igrams" ]; then
    echo "clean obsolete files/folders before rerunning"
    rm -r baselines/ configs/ coregSLC/ geom_reference/ Igrams/ merged/ offsets/ refineSecondaryTiming/ run_files/ SLC/
    rm run_unPackALOS
    cd download
    rm -rf 20*
    mv ARCHIVED_FILES/* .
    cd ..
fi

# prepare SAR data
prepRawALOS.py -i download/ -o SLC -t '' --dual2single
chmod 755 run_unPackALOS
./run_unPackALOS

# stack processing
stackStripMap.py -s SLC/ -d DEM/gsi*.dem -t 1800 -b 1800 -a 20 -r 8 -u snaphu -W interferogram -m 20080212 -f 0.5 --applyWaterMask
submit_run_files4stripmap_stack.py

# mintpy
cd mintpy
smallbaselineApp.py ../KirishimaAlosAT424F620_630.txt
