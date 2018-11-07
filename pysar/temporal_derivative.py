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


############################################################################
USAGE = """
usage:  temporal_derivative.py  timeseries_file 

Calculate the temporal derivative of time-series displacement.
  Useful to check time-dependent deformation.

example:
  temporal_derivative.py  timeseries.h5 
"""

def usage():
    print(USAGE)
    return


############################################################################
def main(argv):
    try:
        timeseries_file = argv[0]
    except:
        usage()
        sys.exit(1)

    obj = timeseries(timeseries_file)
    obj.open(print_msg=False)
    print('reading timeseries data from file: {}'.format(timeseries_file))
    ts_data = obj.read(print_msg=False)

    print('calculate the 1st derivative of timeseries data')
    ts_data_1d = np.zeros(ts_data.shape, np.float32)
    ts_data_1d[1:, :, :] = np.diff(ts_data, n=1, axis=0)

    out_file = '{}_1stDerivative.h5'.format(os.path.splitext(timeseries_file)[0])
    obj_out = timeseries(out_file)
    obj_out.write2hdf5(ts_data_1d, refFile=timeseries_file)

    return out_file


############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
