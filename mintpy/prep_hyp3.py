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

# Extract geospatial and sensor metadata
def general_metadata(fname):
    '''Read/extract attribute data from HyP3 gdal geotiff metadata file and add to metadata dictionary
    Inputs:
        *unw_phase.tif or *corr.tif file name or *dem.tif file name, e.g.
            S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_unw_phase.tif
            S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_corr.tif
            dem.tif
    Output:
        Metadata dictionary (meta)
    '''

    meta = {}

    # Geospatial Data
    ds = gdal.Open(fname, gdal.GA_ReadOnly)
    meta['FILE_LENGTH'] = ds.RasterYSize
    meta['LENGTH'] = ds.RasterYSize
    meta['WIDTH'] = ds.RasterXSize
    geoTrans = ds.GetGeoTransform()
    del ds

    x0 = geoTrans[0]
    y0 = geoTrans[3]
    x_step = geoTrans[1]
    y_step = geoTrans[5]
    x1 = x0 + x_step * meta['WIDTH']
    y1 = y0 + y_step * meta['LENGTH']

    meta['X_FIRST'] = f'{x0:.9f}'
    meta['Y_FIRST'] = f'{y0:.9f}'
    meta['X_STEP'] = f'{x_step:.9f}'
    meta['Y_STEP'] = f'{y_step:.9f}'
    meta['X_UNIT'] = 'METERS'
    meta['Y_UNIT'] = 'METERS'

    meta["LON_REF1"] = x0
    meta["LON_REF2"] = x1
    meta["LON_REF3"] = x0
    meta["LON_REF4"] = x1

    meta["LAT_REF1"] = y0
    meta["LAT_REF2"] = y0
    meta["LAT_REF3"] = y1
    meta["LAT_REF4"] = y1

    # Satellite platform data
    #   Note: HyP3 currently only supports Sentinel-1 data, so Sentinel-1
    #         configuration is hard-coded.
    meta['PROCESSOR'] = 'hyp3'
    meta['PLATFORM'] = 'Sen'
    meta['ANTENNA_SIDE'] = -1
    meta['WAVELENGTH'] = SPEED_OF_LIGHT / sensor.SEN['carrier_frequency']
    meta['HEIGHT'] = 693000.0 # nominal altitude of Sentinel1 orbit

    # Earth radius probably won't be used anywhere for the geocoded data
    meta['EARTH_RADIUS'] = 6337286.638938101

    return(meta)


# Extract data from HyP3 interferogram metadata
def add_ifg_metadata(fname,meta):
    '''Read/extract attribute data from HyP3 metadata file and add to metadata dictionary
    Inputs:
        *unw_phase.tif or *corr.tif file name, e.g.
            S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_unw_phase.tif
            S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_corr.tif
        Metadata dictionary (meta)
    Output:
        Metadata dictionary (meta)
    '''
    directory = os.path.dirname(fname)

    sat, date1_string, date2_string, pol, res, soft, proc, ids, *_ = os.path.basename(fname).split('_')

    job_id = '_'.join([sat, date1_string, date2_string, pol, res, soft, proc, ids])
    date1 = datetime.strptime(date1_string,'%Y%m%dT%H%M%S')
    date2 = datetime.strptime(date2_string,'%Y%m%dT%H%M%S')
    date_avg = date1 + (date2 - date1) / 2

    ifg_meta_file = f'{os.path.join(directory,job_id)}.txt'
    ifg_meta = {}
    with open(ifg_meta_file, 'r') as f:
        for line in f:
            key, value = line.strip().split(': ')
            ifg_meta[key] = value

    # Add relevant data to meta object
    meta['CENTER_LINE_UTC'] = date_avg.strftime('%y%m%d')
    meta['DATE12'] = f'{date1.strftime("%y%m%d")}-{date2.strftime("%y%m%d")}'
    meta['HEADING'] = ifg_meta['Heading']
    meta['ORBIT_DIRECTION'] = 'ASCENDING' if float(meta['HEADING']) > -90 else 'DESCENDING'
    meta['ALOOKS'] = ifg_meta['Azimuth looks']
    meta['RLOOKS'] = ifg_meta['Range looks']
    meta['P_BASELINE_TOP_HDR'] = ifg_meta['Baseline']
    meta['P_BASELINE_BOTTOM_HDR'] = ifg_meta['Baseline']

    return(meta)


#########################################################################

def main(iargs=None):
    # read in arguments
    inps = cmd_line_parse(iargs)

    # for each filename, generate metadata rsc file
    for fname in inps.file:
        meta = general_metadata(fname)

        # add interferogram-specific metadata
        if any([x in fname for x in ['unw_phase','corr']]):
            meta = add_ifg_metadata(fname, meta)

        rsc_file = fname+'.rsc'
        writefile.write_roipac_rsc(meta, out_file=rsc_file)

    return


###################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
