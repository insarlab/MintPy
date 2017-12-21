#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################


import os
import sys
import time
import datetime

import h5py
import numpy as np

import pysar._datetime as ptime
import pysar._readfile as readfile


#####################################################################
def usage():
    print('''usage: sum_epochs.py  timeseries_file  [output_file]

Calculate the Sum of Time Series Displacement per epoch
  For each epoch, referencing it as the master date,
  get a new time series and calculate the temporal
  average displacement.

arguments:
  timeseries_file : string, name/path of timeseries hdf5 file
  output_file     : string, name/path of sum epochs file
                    default: add prefix, 'sum_' to input timeseries file
example:
  sum_epochs.py  timeseries_ECMWF_demErr.h5
  sum_epochs.py  timeseries_ECMWF_demErr_quadratic.h5  sum_timeseries_ECMWF_demErr_quadratic.h5
    ''')
    return


#####################################################################
def main(argv):
    try:
        timeseriesFile = argv[0]
    except:
        usage(); sys.exit(1)

    try:    outname = argv[1]
    except: outname = 'sum_'+timeseriesFile

    ##### Read Timeseries
    atr = readfile.read_attribute(timeseriesFile)
    k = atr['FILE_TYPE']
    print("loading time series: " + timeseriesFile)
    h5timeseries = h5py.File(timeseriesFile)
    dateList = sorted(h5timeseries['timeseries'].keys())
    date_num = len(dateList)
    print('number of acquisitions: %d' % date_num)

    length = int(atr['FILE_LENGTH'])
    width  = int(atr['WIDTH'])
    D = np.zeros((date_num, length*width), np.float32)

    prog_bar = ptime.progress_bar(maxValue=date_num)
    for i in range(date_num):
        date = dateList[i]
        d = h5timeseries['timeseries'].get(date)[:]
        D[i][:] = d.flatten(0)
        prog_bar.update(i+1, suffix=date)
    prog_bar.close()
    h5timeseries.close()

    ##### Calculate Sum
    print('calculating epochs sum ...')
    sumD = np.zeros(D.shape)
    prog_bar.reset()
    for i in range(date_num):
        sumD[j,:] = np.sum(np.abs(D-D[j,:]),0)/date_num
        prog_bar.update(i+1)
    prog_bar.close()

    ## Normalize to 0 and 1
    ## with high atmosphere equal to 0 and no atmosphere equal to 1
    sumD -= np.max(sumD,0)
    sumD *= -1
    sumD /= np.max(sumD,0)
    sumD[np.isnan(sumD)] = 1

    ##### Write sum epochs file
    print('writing to >>> '+outname)
    h5sum = h5py.File(outname,'w')
    group = h5sum.create_group('timeseries')
    prog_bar.reset()
    for i in range(date_num):
        date = dateList[i]
        d = np.reshape(sumD[i][:],[length,width])
        dset = group.create_dataset(date, data=d, compression='gzip')
        prog_bar.update(i+1, suffix=date)
    prog_bar.close()

    for key,value in atr.items():
        group.attrs[key] = value
    h5sum.close()
    print('Done.')


#####################################################################
if __name__ == '__main__':
    main(sys.argv[1:])

