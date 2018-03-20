#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Jul 2017: rewrite using pysay module

import os
import sys
import time
import datetime
import argparse

import h5py
import numpy as np

import pysar.utils.datetime as ptime
import pysar.utils.readfile as readfile


############################################################
EXAMPLE='''example:
 temporal_filter.py timeseries_ECMWF_demErr_refDate.h5
 temporal_filter.py timeseries_ECMWF_demErr_refDate.h5 -t 0.3
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Smoothing Timeseries using moving Gaussian window\n'+\
                                     '  https://en.wikipedia.org/wiki/Gaussian_blur',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('timeseries_file', help='timeseries file to be smoothed.')
    parser.add_argument('-t','--time-win', dest='time_win', type=float, default=0.3,\
                        help='time window in years (Sigma of the assmued Gaussian distribution.)')
    parser.add_argument('-o','--outfile', help='Output file name.')

    inps = parser.parse_args()
    return inps


############################################################
def main(argv):
    inps = cmdLineParse()

    # Basic info
    atr = readfile.read_attribute(inps.timeseries_file)
    k = atr['FILE_TYPE']
    if k not in ['timeseries']:
        sys.exit('ERROR: only timeseries file supported, input is '+k+' file!')

    h5 = h5py.File(inps.timeseries_file,'r')
    date_list = sorted(h5[k].keys())
    date_num = len(date_list)
    length = int(atr['LENGTH'])
    width = int(atr['WIDTH'])
    pixel_num = length*width

    tbase = np.array(ptime.date_list2tbase(date_list)[0], np.float32).reshape((date_num,1))
    tbase /= 365.25

    # Read timeseries
    print('loading time-series ...')
    timeseries = np.zeros((date_num, pixel_num))
    prog_bar = ptime.progress_bar(maxValue=date_num)
    for i in range(date_num):
        date = date_list[i]
        d = h5[k].get(date)[:]
        timeseries[i,:] = d.flatten(0)
        prog_bar.update(i+1, suffix=date)
    del d
    h5.close()
    prog_bar.close()

    # Smooth timeseries with moving window in time
    print('smoothing time-series using moving gaussian window with size of %.1f years' % inps.time_win)
    timeseries_filt = np.zeros((date_num, pixel_num))
    prog_bar = ptime.progress_bar(maxValue=date_num)
    for i in range(date_num):
        date = date_list[i]
        # Weight from Gaussian (normal) distribution in time
        t_diff = tbase[i] - tbase
        weight = np.exp(-0.5*(t_diff**2)/(inps.time_win**2))
        weight /= np.sum(weight)
        weightMat = np.tile(weight, (1,pixel_num))
        # Smooth the current acquisition - moving window in time one by one
        timeseries_filt[i,:] = np.sum(timeseries*weightMat, 0)
        prog_bar.update(i+1, suffix=date)
    del weightMat
    del timeseries
    prog_bar.close()

    # Write smoothed timeseries file
    try:    ref_date = atr['REF_DATE']
    except: ref_date = date_list[0]
    ref_date_idx = date_list.index(ref_date)
    print('reference date: '+ref_date)
    print('reference date index: '+str(ref_date_idx))
    ref_data = np.reshape(timeseries_filt[ref_date_idx,:], [length, width])

    if not inps.outfile:
        inps.outfile = os.path.splitext(inps.timeseries_file)[0]+'_smooth.h5'
    print('writing >>> '+inps.outfile)
    print('number of acquisitions: '+str(date_num))

    h5out = h5py.File(inps.outfile, 'w')
    group = h5out.create_group(k)
    prog_bar = ptime.progress_bar(maxValue=date_num)
    for i in range(date_num):
        date = date_list[i]
        data = np.reshape(timeseries_filt[i,:], [length, width])
        dset = group.create_dataset(date, data=data-ref_data, compression='gzip')
        prog_bar.update(i+1, suffix=date)
    for key,value in iter(atr.items()):
        group.attrs[key] = value
    h5out.close()
    prog_bar.close()

    print('Done.')
    return inps.outfile


############################################################
if __name__ == '__main__':
    main(sys.argv[1:])
