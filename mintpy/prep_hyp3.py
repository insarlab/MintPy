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

try:
    from osgeo import gdal
except ImportError:
    raise ImportError('Can not import gdal!')

from mintpy.objects import sensor
from mintpy.utils import writefile, utils as ut


SPEED_OF_LIGHT = 299792458  # m/s

EXAMPLE = """example:
  prep_hyp3.py  interferograms/*/*unw_phase.tif
  prep_hyp3.py  interferograms/*/*corr.tif
  prep_hyp3.py  interferograms/*/*unw_phase_clip.tif
  prep_hyp3.py  interferograms/*/*corr_clip.tif
  prep_hyp3.py  interferograms/*/*clip.tif
  prep_hyp3.py  dem.tif
  prep_hyp3.py  dem_clip.tif
"""

DESCRIPTION = """
  For each interferogram, the unwrapped interferogram, coherence, and metadata the file name is required e.g.:
  1) S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_unw_phase.tif
  2) S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_corr.tif
  3) S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2.txt

  A DEM is also needed e.g.:
  1) dem.tif

  This script will read these files, read the geospatial metadata from GDAL,
  find the corresponding HyP3 metadata file (for interferograms and coherence),
  and write to a ROI_PAC .rsc metadata file with the same name as the input file with suffix .rsc,
  e.g. S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_unw_phase.tif.rsc

  Here is an example of how your HyP3 files should look:

  Before loading:
      For each interferogram, 3 files are needed:
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_unw_phase.tif
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_corr.tif
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2.txt
      For the geometry file 1 file is needed:
          dem.tif (a DEM with any name)

  After running prep_hyp3.py:
      For each interferogram:
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_unw_phase.tif
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_unw_phase.tif.rsc
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_corr.tif
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_corr.tif.rsc
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2.txt
      For the input DEM geometry file:
          dem.tif
          dem.tif.rsc

  Notes:
    HyP3 currently only supports generation of Sentinel-1 interferograms, so
    some Sentinel-1 metadata is hard-coded. If HyP3 adds processing of interferograms
    from other satellites, changes will be needed. 
"""

#########################################################################

# CMD parsing
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

# extract geospatial and add to meta
def add_geospatial_metadata(fname, meta):
    '''Read/extract attribute data from HyP3 gdal geotiff metadata and add to metadata dictionary
    Inputs:
        *unw_phase.tif, *corr.tif file name or *dem.tif file name, e.g.
            S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_unw_phase.tif
            S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_corr.tif
            dem.tif
        Metadata dictionary (meta)
    Output:
        Metadata dictionary (meta)
    '''

    # geospatial Data
    ds = gdal.Open(fname, gdal.GA_ReadOnly)
    transform = ds.GetGeoTransform()

    x_step = abs(transform[1])
    y_step = abs(transform[5]) * -1.

    W = transform[0] - x_step / 2.
    N = transform[3] - y_step / 2.
    E = W + x_step * ds.RasterXSize
    S = N + y_step * ds.RasterYSize

    meta['FILE_LENGTH'] = ds.RasterYSize
    meta['LENGTH'] = ds.RasterYSize
    meta['WIDTH'] = ds.RasterXSize

    meta['NORTH'] = N
    meta['SOUTH'] = S
    meta['EAST'] = E
    meta['WEST'] = W

    meta['X_FIRST'] = W
    meta['Y_FIRST'] = N
    meta['X_STEP'] = x_step
    meta['Y_STEP'] = y_step
    meta['X_UNIT'] = 'meters'
    meta['Y_UNIT'] = 'meters'

    # Earth radius probably won't be used anywhere for the geocoded data
    meta['EARTH_RADIUS'] = 6337286.638938101

    del ds

    return(meta)

# extract data from HyP3 interferogram metadata
def add_hyp3_metadata(fname,meta,is_ifg=True):
    '''Read/extract attribute data from HyP3 metadata file and add to metadata dictionary
    Inputs:
        *unw_phase.tif or *corr.tif file name, *dem.tif e.g.
            S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_unw_phase.tif
            S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_corr.tif
            dem.tif
        Metadata dictionary (meta)
    Output:
        Metadata dictionary (meta)
    '''

    if is_ifg:
        # determine interferogram pair info and hyp3 metadata file name
        sat, date1_string, date2_string, pol, res, soft, proc, ids, *_ = os.path.basename(fname).split('_')

        job_id = '_'.join([sat, date1_string, date2_string, pol, res, soft, proc, ids])
        directory = os.path.dirname(fname)
        meta_file = f'{os.path.join(directory,job_id)}.txt'

        date1 = datetime.strptime(date1_string,'%Y%m%dT%H%M%S')
        date2 = datetime.strptime(date2_string,'%Y%m%dT%H%M%S')
        date_avg = date1 + (date2 - date1) / 2
        date_avg_seconds = (date_avg - date_avg.replace(hour=0, minute=0, second=0, microsecond=0)).total_seconds()

        # open and read hyp3 metadata
        hyp3_meta = {}
        with open(meta_file, 'r') as f:
            for line in f:
                key, value = line.strip().split(': ')
                hyp3_meta[key] = value

        # add relevant data to meta object
        meta['CENTER_LINE_UTC'] = date_avg_seconds
        meta['DATE12'] = f'{date1.strftime("%y%m%d")}-{date2.strftime("%y%m%d")}'
        meta['HEADING'] = hyp3_meta['Heading']
        meta['ALOOKS'] = hyp3_meta['Azimuth looks']
        meta['RLOOKS'] = hyp3_meta['Range looks']
        meta['P_BASELINE_TOP_HDR'] = hyp3_meta['Baseline']
        meta['P_BASELINE_BOTTOM_HDR'] = hyp3_meta['Baseline']

    else:
        meta_file = os.path.splitext(fname)[0]+'.txt'

        # open and read hyp3 metadata
        hyp3_meta = {}
        with open(meta_file, 'r') as f:
            for line in f:
                key, value = line.strip().split(': ')
                hyp3_meta[key] = value

        # add relevant data to meta object
        meta['HEADING'] = hyp3_meta['Heading']

    return(meta)

# add sentinel-1 metadata and remove unnessecary metadata
def add_sentinel1_metadata(meta):
    '''Add Sentinel-1 attribute data to metadata dictionary
    Inputs:
        Metadata dictionary (meta)
    Output:
        Metadata dictionary (meta)
    '''

    # note: HyP3 currently only supports Sentinel-1 data, so Sentinel-1
    #       configuration is hard-coded.

    meta['PROCESSOR'] = 'hyp3'
    meta['PLATFORM'] = 'Sen'
    meta['ANTENNA_SIDE'] = -1
    meta['WAVELENGTH'] = SPEED_OF_LIGHT / sensor.SEN['carrier_frequency']
    meta['HEIGHT'] = 693000.0 # nominal altitude of Sentinel1 orbit

    # add LAT/LON_REF1/2/3/4 based on whether satellite ascending or descending
    meta['ORBIT_DIRECTION'] = 'ASCENDING' if float(meta['HEADING']) > -90 else 'DESCENDING'

    if meta['ORBIT_DIRECTION'] == 'ASCENDING':
        meta['LAT_REF1'] = str(meta['SOUTH'])
        meta['LAT_REF2'] = str(meta['SOUTH'])
        meta['LAT_REF3'] = str(meta['NORTH'])
        meta['LAT_REF4'] = str(meta['NORTH'])
        meta['LON_REF1'] = str(meta['WEST'])
        meta['LON_REF2'] = str(meta['EAST'])
        meta['LON_REF3'] = str(meta['WEST'])
        meta['LON_REF4'] = str(meta['EAST'])
    else:
        meta['LAT_REF1'] = str(meta['NORTH'])
        meta['LAT_REF2'] = str(meta['NORTH'])
        meta['LAT_REF3'] = str(meta['SOUTH'])
        meta['LAT_REF4'] = str(meta['SOUTH'])
        meta['LON_REF1'] = str(meta['EAST'])
        meta['LON_REF2'] = str(meta['WEST'])
        meta['LON_REF3'] = str(meta['EAST'])
        meta['LON_REF4'] = str(meta['WEST'])

    del meta['NORTH'], meta['SOUTH'], meta['EAST'], meta['WEST']

    return(meta)


#########################################################################

def main(iargs=None):
    # read in arguments
    inps = cmd_line_parse(iargs)

    # for each filename, generate metadata rsc file
    for fname in inps.file:
        meta = {}
        meta = add_geospatial_metadata(fname,meta)

        if any([x in fname for x in ['unw_phase','corr']]):
            meta = add_hyp3_metadata(fname, meta, is_ifg=True)
        else:
            meta = add_hyp3_metadata(fname, meta, is_ifg=False)

        meta = add_sentinel1_metadata(meta)
        rsc_file = fname+'.rsc'
        writefile.write_roipac_rsc(meta, out_file=rsc_file)

    return


###################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
