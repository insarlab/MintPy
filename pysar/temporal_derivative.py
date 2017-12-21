#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Jul 2017: re-write using pysar modules


import os
import sys
import time
import datetime

import h5py
import numpy as np

import pysar._datetime as ptime
import pysar._readfile as readfile


############################################################################
def usage():
    print('''usage:  temporal_derivative.py  timeseries_file 

Calculate the temporal derivative of time-series displacement.
  Useful to check time-dependent deformation.

example:
  temporal_derivative.py  timeseries.h5 
    ''')
    return


############################################################################
def main(argv):
    try:
        timeseries_file = argv[0]
    except:
        usage() ; sys.exit(1)

    # Basic info
    atr = readfile.read_attribute(timeseries_file)
    k = atr['FILE_TYPE']
    length = int(atr['FILE_LENGTH'])
    width = int(atr['WIDTH'])

    ##### Read time-series
    print("loading time series: " + timeseries_file)
    h5 = h5py.File(timeseries_file)
    date_list = sorted(h5[k].keys())
    date_num = len(date_list)
    pixel_num = length*width

    tbase = np.array(ptime.date_list2tbase(date_list)[0], np.float32)

    prog_bar = ptime.progress_bar(maxValue=date_num)
    timeseries = np.zeros((date_num, pixel_num),np.float32)
    for i in range(date_num):
        date = date_list[i]
        d = h5[k].get(date)[:]
        timeseries[i,:] = d.flatten(0)
        prog_bar.update(i+1, suffix=date)
    prog_bar.close()
    del d
    h5.close()

    ##### Calculate 1st and 2nd temporal derivatives
    print("calculating temporal 1st derivative ... ")
    timeseries_1st = np.zeros((date_num-1,pixel_num),np.float32)
    for i in range(date_num-1):
        timeseries_1st[i][:] = timeseries[i+1][:] - timeseries[i][:]

    print("calculating temporal 2nd derivative")
    timeseries_2nd = np.zeros((date_num-2,pixel_num),np.float32)
    for i in range(date_num-2):
        timeseries_2nd[i][:] = timeseries_1st[i+1][:] - timeseries_1st[i][:]

    ##### Write 1st and 2nd temporal derivatives
    outfile1 = os.path.splitext(timeseries_file)[0]+'_1stDerivative.h5'
    print('writing >>> '+outfile1)
    h5out = h5py.File(outfile1, 'w')
    group = h5out.create_group(k)

    prog_bar = ptime.progress_bar(maxValue=date_num-1)
    for i in range(date_num-1):
        date = date_list[i+1]
        dset = group.create_dataset(date, data=np.reshape(timeseries_1st[i][:],[length,width]), compression='gzip')
        prog_bar.update(i+1, suffix=date)
    for key,value in atr.items():
        group.attrs[key] = value
    prog_bar.close()
    h5out.close()

    outfile2 = os.path.splitext(timeseries_file)[0]+'_2ndDerivative.h5'
    print('writing >>> '+outfile2)
    h5out = h5py.File(outfile2, 'w')
    group = h5out.create_group(k)

    prog_bar = ptime.progress_bar(maxValue=date_num-2)
    for i in range(date_num-2):
        date = date_list[i+2]
        dset = group.create_dataset(date, data=np.reshape(timeseries_2nd[i][:],[length,width]), compression='gzip')
        prog_bar.update(i+1, suffix=date)
    for key,value in atr.items():
        group.attrs[key] = value
    prog_bar.close()
    h5out.close()

    print('Done.')
    return outfile1, outfile2


############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])

