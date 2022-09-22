#! /bin/bash
# Extract input datasets in mintpy/inputs for an area of interest (AOI) from the one of the whole frame
# Parameters: SNWE : string for the AOI
#             step : string/number, output resolution in degree
#             inputs_dir_src : source mintpy/inputs folder for the whole track in radar coordinates
# Returns:    inputs_dir_dst : destination mintpy/inputs folder for the AOI in geo coordinates
# Author: Zhang Yunjun, 2019-05-18

# setup input/output directories
inputs_dir_src=$HOME'/insarlab/Kirishima/Alos2DT23F2970/mintpy/inputs'        #input dir
inputs_dir_dst=$HOME'/insarlab/Kirishima/ShinmoedakeAlos2DT23/mintpy/inputs'  #output dir

# setup AOI
SNWE="31.88 31.94 130.85 130.91"
step="0.000185185"   #degrees
# degrees --> meters on equator
# 0.000925926 --> 100 m
# 0.000555556 --> 60 m
# 0.000277778 --> 30 m
# 0.000185185 --> 20 m
# 0.000092593 --> 10 m


###########################################################################
# create destination directory
if [ ! -d $inputs_dir_dst ]; then
    echo "create "$inputs_dir_dst
    mkdir -p $inputs_dir_dst
fi

# geocode - prepare script options
lut_file="$inputs_dir_src/geometryRadar.h5"
geocode_opt="-l $lut_file --bbox $SNWE --lat-step -$step --lon-step $step"

echo "geocode & subset - geometry file"
src_file=$inputs_dir_src"/geometryRadar.h5"
dst_file=$inputs_dir_dst"/geometryGeo.h5"
echo "geocode.py $src_file -o $dst_file $geocode_opt"
geocode.py $src_file -o $dst_file $geocode_opt

echo "geocode & subset - ifgramStack file"
src_file=$inputs_dir_src"/ifgramStack.h5"
dst_file=$inputs_dir_dst"/ifgramStack.h5"
echo "geocode.py $src_file -o $dst_file $geocode_opt"
geocode.py $src_file -o $dst_file $geocode_opt

# subset - prepare script options
# split SNWE into four variables to be used by subset.py
S="$(cut -d' ' -f1 <<<"$SNWE")"
N="$(cut -d' ' -f2 <<<"$SNWE")"
W="$(cut -d' ' -f3 <<<"$SNWE")"
E="$(cut -d' ' -f4 <<<"$SNWE")"
subset_opt=" --lat $S $N --lon $W $E"

echo "subset - DEM"
src_file=$inputs_dir_src"/../../../DEM/gsi10m.dem.wgs84"   #adjust filename for specific dataset
dst_file=$inputs_dir_dst"/gsi10m.dem.wgs84"                #adjust filename for specific dataset
echo "subset.py $src_file -o $dst_file $subset_opt"
subset.py $src_file -o $dst_file $subset_opt

echo "Done."
