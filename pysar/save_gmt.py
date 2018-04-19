#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Dec 2015: add support for ROI_PAC product


import os
import sys
import argparse

import h5py
import numpy as np
import scipy.io.netcdf as netcdf

import pysar.utils.readfile as readfile
import pysar.view as pv
from pysar.utils.readfile import multi_group_hdf5_file, multi_dataset_hdf5_file, single_dataset_hdf5_file


####################################################################################
def write_gmt_simple(lons, lats, z, fname, title='default', name='z', scale=1.0, offset=0, units='meters'):
    '''Writes a simple GMT grd file with one array.
    This is based on the gdal2grd.py script found at:
        http://http://www.vso.cape.com/~nhv/files/python/gdal/gdal2grd.py
    
    .. Args:
        
        * lons     -> 1D Array of lon values
        * lats     -> 1D Array of lat values
        * z        -> 2D slice to be saved
        * fname    -> Output file name
        
    .. Kwargs:
        
        * title    -> Title for the grd file
        * name     -> Name of the field in the grd file
        * scale    -> Scale value in the grd file
        * offset   -> Offset value in the grd file
        
    .. Returns:
        
        * None'''
    fid = netcdf.netcdf_file(fname,'w')

    ####Create a dimension variable
    fid.createDimension('side',2)
    fid.createDimension('xysize',np.prod(z.shape))

    ####Range variables
    fid.createVariable('x_range','d',('side',))
    fid.variables['x_range'].units = 'degrees'

    fid.createVariable('y_range','d',('side',))
    fid.variables['y_range'].units = 'degrees'

    fid.createVariable('z_range','d',('side',))
    fid.variables['z_range'].units = units

    #####Spacing
    fid.createVariable('spacing','d',('side',))
    fid.createVariable('dimension','i4',('side',))

    fid.createVariable('z','f',('xysize',))
    fid.variables['z'].long_name = name
    fid.variables['z'].scale_factor = scale
    fid.variables['z'].add_offset = offset
    fid.variables['z'].node_offset=0

    fid.title = title
    fid.source = 'PySAR v1.2'

    #####Filling in the actual data
    fid.variables['x_range'][0] = lons[0]
    fid.variables['x_range'][1] = lons[-1]
    fid.variables['spacing'][0] = lons[1]-lons[0]

    fid.variables['y_range'][0] = lats[0]
    fid.variables['y_range'][1] = lats[-1]
    fid.variables['spacing'][1] = lats[1]-lats[0]

    #####Range
    zmin = np.nanmin(z)
    zmax = np.nanmax(z)

    fid.variables['z_range'][0] = zmin
    fid.variables['z_range'][1] = zmax

    fid.variables['dimension'][:] = z.shape[::-1]
    fid.variables['z'][:] = np.flipud(z).flatten()
    fid.close()

    ############################################################
    # Program is part of GIAnT v1.0                            #
    # Copyright 2012, by the California Institute of Technology#
    # Contact: earthdef@gps.caltech.edu                        #
    ############################################################


def get_geo_lat_lon(atr):
    X_FIRST = float(atr['X_FIRST'])
    Y_FIRST = float(atr['Y_FIRST'])
    X_STEP = float(atr['X_STEP'])
    Y_STEP = float(atr['Y_STEP'])
    W = float(atr['WIDTH'])
    L = float(atr['LENGTH'])
    Y_END = Y_FIRST + L*Y_STEP
    X_END = X_FIRST + W*X_STEP

    X = np.linspace(X_FIRST,X_END,W)
    Y = np.linspace(Y_FIRST,Y_END,L)
    #XI,YI = np.meshgrid(X,Y)

    return Y,X


def write_grd_file(data, atr, fname_out=None):
    '''Write GMT .grd file for input data matrix, using giant._gmt module.
    Inputs:
        data - 2D np.array in int/float, data matrix to write
        atr  - dict, attributes of input data matrix
        fname_out - string, output file name
    Output:
        fname_out - string, output file name
    '''
    if not fname_out:
        fname_out = os.path.splitext(atr['FILE_PATH'])[0]+'.grd'

    # Get 1D array of lats and lons
    lats, lons = get_geo_lat_lon(atr)

    # writing
    print('writing >>> '+fname_out)
    write_gmt_simple(lons, np.flipud(lats), np.flipud(data), fname_out,\
                         title='default', name=atr['FILE_TYPE'], scale=1.0, offset=0, units=atr['UNIT'])
    return fname_out



####################################################################################
EXAMPLE='''example:
  save_gmt.py  geo_velocity.h5
  save_gmt.py  geo_timeseries.h5     20071031
  save_gmt.py  geo_timeseries.h5
  save_gmt.py  geo_filt_100608-101024-sim_HDR_16rlks_c10.unw
  save_gmt.py  gsi10m.dem
'''

def create_parser():
    parser = argparse.ArgumentParser(description='Export geocoded file to GMT grd file',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='file to be converted, in geo coordinate.')
    parser.add_argument('dset', nargs='?', help='date of timeseries, or date12 of interferograms to be converted')
    parser.add_argument('-o','--output', dest='outfile', help='output file base name. Extension is fixed with .kmz')
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


####################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    ##### 1. Read data
    atr = readfile.read_attribute(inps.file)
    k = atr['FILE_TYPE']
    print('Input file is '+k)

    # Check: file in geo coord
    if 'X_FIRST' not in atr.keys():
        sys.exit('ERROR: Input file is not geocoded.')

    # Check: dset is required for multi_dataset/group files
    if not inps.dset:
        if k in multi_group_hdf5_file:
            print("No date/date12 input.\nIt's required for "+k+" file")
            sys.exit(1)
        elif k in multi_dataset_hdf5_file:
            print('No input date ..., continue to convert the last date of time-series.')
            h5 = h5py.File(inps.file, 'r')
            date_list = sorted(h5[k].keys())
            h5.close()
            inps.dset = date_list[-1]

    # Read data
    data, atr = readfile.read(inps.file, datasetName=inps.dset)

    # Output filename
    if not inps.outfile:
        inps.outfile = pv.auto_figure_title(inps.file, inps.dset, vars(inps))
    inps.outfile = os.path.splitext(inps.outfile)[0]+'.grd'

    ##### 2. Write GMT .grd file
    inps.outfile = write_grd_file(data, atr, inps.outfile)
    print('Done.')
    return inps.outfile


####################################################################################
if __name__ == '__main__':
    main()

