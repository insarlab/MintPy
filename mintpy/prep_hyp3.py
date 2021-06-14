#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Forrest Williams, Mar 2021                       #
############################################################


import os
import sys
import argparse
from datetime import datetime
from mintpy.objects import sensor
from mintpy.utils import readfile, writefile, utils as ut


SPEED_OF_LIGHT = 299792458  # m/s


#########################################################################
EXAMPLE = """example:
  prep_hyp3.py  interferograms/*/*unw_phase_clip.tif
  prep_hyp3.py  interferograms/*/*corr_clip.tif
  prep_hyp3.py  interferograms/*/*dem_clip.tif
  prep_hyp3.py  interferograms/*/*inc_map_clip.tif
  prep_hyp3.py  interferograms/*/*clip.tif
"""

DESCRIPTION = """
  For each interferogram, the unwrapped interferogram, coherence, and metadata the file name is required e.g.:
  1) S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_unw_phase.tif
  2) S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_corr.tif
  3) S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2.txt

  A DEM filename is needed and a incidence angle filename is recommended  e.g.:
  1) S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_dem.tif
  2) S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_inc_map.tif

  This script will read these files, read the geospatial metadata from GDAL,
  find the corresponding HyP3 metadata file (for interferograms and coherence),
  and write to a ROI_PAC .rsc metadata file with the same name as the input file with suffix .rsc,
  e.g. S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_unw_phase.tif.rsc

  Here is an example of how your HyP3 files should look:

  Before loading:
      For each interferogram, 3 files are needed:
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_unw_phase_clip.tif
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_corr_clip.tif
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2.txt
      For the geometry file 2 file are recommended:
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_dem_clip.tif     (required)
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_inc_map_clip.tif (optional but recommended)

  After running prep_hyp3.py:
      For each interferogram:
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_unw_phase_clip.tif
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_unw_phase_clip.tif.rsc
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_corr_clip.tif
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_corr_clip.tif.rsc
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2.txt
      For the input geometry files:
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_dem_clip.tif
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_dem_clip.tif.rsc
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_inc_map_clip.tif
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_inc_map_clip.tif.rsc

  Notes:
    HyP3 currently only supports generation of Sentinel-1 interferograms, so
    some Sentinel-1 metadata is hard-coded. If HyP3 adds processing of interferograms
    from other satellites, changes will be needed. 
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Prepare attributes file for HyP3 InSAR product.\n'+
                                     DESCRIPTION,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs='+', help='HyP3 file(s)')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    inps.file = ut.get_file_list(inps.file, abspath=True)
    return inps


#########################################################################
def add_hyp3_metadata(fname,meta,is_ifg=True):
    '''Read/extract attribute data from HyP3 metadata file and add to metadata dictionary
    Inputs:
        *unw_phase.tif, *corr.tif file name, *dem.tif, *inc_map.tif, e.g.
            S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_unw_phase_clip.tif
            S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_corr_clip.tif
            S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_dem_clip.tif
        Metadata dictionary (meta)
    Output:
        Metadata dictionary (meta)
    '''

    # determine interferogram pair info and hyp3 metadata file name
    sat, date1_string, date2_string, pol, res, soft, proc, ids, *_ = os.path.basename(fname).split('_')
    job_id = '_'.join([sat, date1_string, date2_string, pol, res, soft, proc, ids])
    directory = os.path.dirname(fname)
    meta_file = f'{os.path.join(directory,job_id)}.txt'

    # open and read hyp3 metadata
    hyp3_meta = {}
    with open(meta_file, 'r') as f:
        for line in f:
            key, value = line.strip().split(': ')
            hyp3_meta[key] = value

    # add universal hyp3 metadata
    meta['PROCESSOR'] = 'hyp3'
    meta['HEADING'] = hyp3_meta['Heading']
    meta['ALOOKS'] = hyp3_meta['Azimuth looks']
    meta['RLOOKS'] = hyp3_meta['Range looks']
    meta['P_BASELINE_TOP_HDR'] = hyp3_meta['Baseline']
    meta['P_BASELINE_BOTTOM_HDR'] = hyp3_meta['Baseline']
    meta['EARTH_RADIUS'] = hyp3_meta['Earth radius at nadir']
    meta['HEIGHT'] = hyp3_meta['Spacecraft height']

    # add LAT/LON_REF1/2/3/4 based on whether satellite ascending or descending
    meta['ORBIT_DIRECTION'] = 'ASCENDING' if float(meta['HEADING']) > -90 else 'DESCENDING'
    N = float(meta['Y_FIRST'])
    W = float(meta['X_FIRST'])
    S = N + float(meta['Y_STEP']) * int(meta['LENGTH'])
    E = W + float(meta['X_STEP']) * int(meta['WIDTH'])

    if meta['ORBIT_DIRECTION'] == 'ASCENDING':
        meta['LAT_REF1'] = str(S)
        meta['LAT_REF2'] = str(S)
        meta['LAT_REF3'] = str(N)
        meta['LAT_REF4'] = str(N)
        meta['LON_REF1'] = str(W)
        meta['LON_REF2'] = str(E)
        meta['LON_REF3'] = str(W)
        meta['LON_REF4'] = str(E)
    else:
        meta['LAT_REF1'] = str(N)
        meta['LAT_REF2'] = str(N)
        meta['LAT_REF3'] = str(S)
        meta['LAT_REF4'] = str(S)
        meta['LON_REF1'] = str(E)
        meta['LON_REF2'] = str(W)
        meta['LON_REF3'] = str(E)
        meta['LON_REF4'] = str(W)

    # note: HyP3 currently only supports Sentinel-1 data, so Sentinel-1
    #       configuration is hard-coded.
    meta['PLATFORM'] = 'Sen'
    meta['ANTENNA_SIDE'] = -1
    meta['WAVELENGTH'] = SPEED_OF_LIGHT / sensor.SEN['carrier_frequency']

    # note: HyP3 (incidence) angle datasets are in the unit of radian
    # which is different from the isce-2 convention of degree
    if 'inc_map' in os.path.basename(fname):
        meta['UNIT'] = 'radian'

    # add metadata that is only relevant to interferogram files
    if is_ifg:
        date1 = datetime.strptime(date1_string,'%Y%m%dT%H%M%S')
        date2 = datetime.strptime(date2_string,'%Y%m%dT%H%M%S')
        date_avg = date1 + (date2 - date1) / 2
        date_avg_seconds = (date_avg - date_avg.replace(hour=0, minute=0, second=0, microsecond=0)).total_seconds()

        meta['CENTER_LINE_UTC'] = date_avg_seconds
        meta['DATE12'] = f'{date1.strftime("%y%m%d")}-{date2.strftime("%y%m%d")}'

    return(meta)


#########################################################################
def main(iargs=None):
    # read in arguments
    inps = cmd_line_parse(iargs)

    # for each filename, generate metadata rsc file
    for fname in inps.file:
        is_ifg = any([x in fname for x in ['unw_phase','corr']])
        meta = readfile.read_gdal_vrt(fname)
        meta = add_hyp3_metadata(fname, meta, is_ifg=is_ifg)

        # write
        rsc_file = fname+'.rsc'
        writefile.write_roipac_rsc(meta, out_file=rsc_file)

    return


###################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
