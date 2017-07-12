#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Dec 2015: add support for ROI_PAC product


import os
import sys
import argparse

import _gmt as gmt
import h5py
import numpy as np

import pysar._readfile as readfile
import pysar.view as pview
from pysar._readfile import multi_group_hdf5_file, multi_dataset_hdf5_file, single_dataset_hdf5_file


####################################################################################
def get_geo_lat_lon(atr):
    X_FIRST = float(atr['X_FIRST'])
    Y_FIRST = float(atr['Y_FIRST'])
    X_STEP = float(atr['X_STEP'])
    Y_STEP = float(atr['Y_STEP'])
    W = float(atr['WIDTH'])
    L = float(atr['FILE_LENGTH'])
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
    print 'writing >>> '+fname_out
    gmt.write_gmt_simple(lons, np.flipud(lats), np.flipud(data), fname_out,\
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

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Export geocoded file to GMT grd file',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='file to be converted, in geo coordinate.')
    parser.add_argument('epoch', nargs='?', help='date of timeseries, or date12 of interferograms to be converted')
    parser.add_argument('-o','--output', dest='outfile', help='output file base name. Extension is fixed with .kmz')

    inps = parser.parse_args()
    return inps


####################################################################################
def main(argv):
    inps = cmdLineParse()

    ##### 1. Read data
    atr = readfile.read_attribute(inps.file)
    k = atr['FILE_TYPE']
    print 'Input file is '+k

    # Check: file in geo coord
    if 'X_FIRST' not in atr.keys():
        sys.exit('ERROR: Input file is not geocoded.')

    # Check: epoch is required for multi_dataset/group files
    if not inps.epoch:
        if k in multi_group_hdf5_file:
            print "No date/date12 input.\nIt's required for "+k+" file"
            sys.exit(1)
        elif k in multi_dataset_hdf5_file:
            print 'No input date ..., continue to convert the last date of time-series.'
            h5 = h5py.File(inps.file, 'r')
            date_list = sorted(h5[k].keys())
            h5.close()
            inps.epoch = date_list[-1]

    # Read data
    data, atr = readfile.read(inps.file, (), inps.epoch)

    # Output filename
    if not inps.outfile:
        inps.outfile = pview.auto_figure_title(inps.file, inps.epoch, vars(inps))
    inps.outfile = os.path.splitext(inps.outfile)[0]+'.grd'

    ##### 2. Write GMT .grd file
    inps.outfile = write_grd_file(data, atr, inps.outfile)
    print 'Done.'
    return inps.outfile


####################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])




