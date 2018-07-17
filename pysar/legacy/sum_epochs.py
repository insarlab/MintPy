#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013-2018, Heresh Fattahi, Zhang Yunjun     #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import os
import sys
import numpy as np
from pysar.objects import timeseries
from pysar.utils import writefile


#####################################################################
USAGE = """
usage: sum_epochs.py  timeseries_file  [output_file]

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
"""


def usage():
    print(USAGE)
    return


#####################################################################
def main(argv):
    try:
        timeseries_file = argv[0]
    except:
        usage()
        sys.exit(1)

    try:
        out_file = argv[1]
    except:
        out_file = 'sum_'+timeseries_file

    # Read Timeseries
    obj = timeseries(timeseries_file)
    obj.open()
    D = obj.read().reshape(obj.numDate, -1)

    # Calculate Sum
    sumD = np.zeros(D.shape)
    for i in range(obj.numDate):
        sumD[i, :] = np.sum(np.abs(D - D[i, :]), axis=0) / obj.numDate
        sys.stdout.write('\rcalculating epochs sum {}/{} ...'.format(i+1, obj.numDate))
        sys.stdout.flush()
    print('')
    del D

    # Normalize to 0 and 1
    # with high atmosphere equal to 0 and no atmosphere equal to 1
    sumD -= np.max(sumD, 0)
    sumD *= -1
    sumD /= np.max(sumD, 0)
    sumD[np.isnan(sumD)] = 1

    # Write sum epochs file
    sumD = np.reshape(sumD, (obj.numDate, obj.length, obj.width))
    atr = dict(obj.metadata)
    atr['UNIT'] = '1'
    writefile.write(sumD, out_file=out_file, metadata=atr, ref_file=timeseries_file)
    print('Done.')


#####################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
