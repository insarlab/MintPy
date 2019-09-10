#! /bin/sh
# generate multilooked geometry files for ISCE/stripmapStack product
# Author: Zhang Yunjun, 2019-03-10


if [ ! $# -eq 2 ]; then
    echo 'No input argument for alks and rlks! Example:'
    echo './run_multilook_geometry 20 4'
    exit 1
fi

alks=$1
rlks=$2
echo 'number of multilooks in azimuth direction: '${alks}
echo 'number of multilooks in range   direction: '${rlks}
in_dir='./merged/geom_master'
out_dir='./geom_master'

#geom_master folder
if [ ! -d ${out_dir} ]; then
    echo 'create '${out_dir}
    mkdir ${out_dir}
fi

#multilook geometry files
for fbase in hgt incLocal lat lon los shadowMask
do
    fname=${fbase}'.rdr'
    echo 'multilook '${fname}
    looks.py -a ${alks} -r ${rlks} -i ${in_dir}/${fname} -o ${out_dir}/${fname}
    cp ${in_dir}/${fname}.xml ${out_dir}/${fname}.full.xml
    cp ${in_dir}/${fname}.vrt ${out_dir}/${fname}.full.vrt
done

