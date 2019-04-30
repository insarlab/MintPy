#!/bin/sh
cd /Users/yunjunz/insarlab/Galapagos/GalapagosSenDT128_30CONN/PYSAR
#reference_point.py INPUTS/ifgramStack.h5 -l -0.81 -L -91.14 --lookup INPUTS/geometryRadar.h5
for num in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20
do
    out_dir=./${num}CONN
    if [ ! -d $out_dir ]; then
        mkdir $out_dir
    fi

    #modify_network.py INPUTS/ifgramStack.h5 --max-conn-num $num -t template.txt
    #plot_network.py INPUTS/ifgramStack.h5 --nodrop --mask None -t template.txt --nodisplay
    #ifgram_inversion.py INPUTS/ifgramStack.h5 -w var -i unwrapPhase
    #mv timeseries.h5 temporalCoherence.h5 numInvIfgram.h5 $out_dir/

    #view.py $out_dir/temporalCoherence.h5 -o ./PIC/tcoh_${num}.png --nodisplay
    #mv CoherenceMatrix.pdf ./PIC/CoherenceMatrix_${num}.pdf
done


