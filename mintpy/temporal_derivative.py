############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, Zhang Yunjun, 2013               #
############################################################


import os
import numpy as np
from mintpy.objects import timeseries


############################################################################
def temporal_derivative(inputfile, outfile):
    # read data
    obj = timeseries(inputfile)
    obj.open(print_msg=False)
    print('reading timeseries data from file: {}'.format(inputfile))
    ts_data = obj.read(print_msg=False)

    # calculation
    print('calculate the 1st derivative of timeseries data')
    ts_data_1d = np.zeros(ts_data.shape, np.float32)
    ts_data_1d[1:, :, :] = np.diff(ts_data, n=1, axis=0)

    # write to file
    if not outfile:
        outfile = '{}_1stDiff.h5'.format(os.path.splitext(inputfile)[0])
    obj_out = timeseries(outfile)
    obj_out.write2hdf5(ts_data_1d, refFile=inputfile)
