############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, Zhang Yunjun, 2013               #
############################################################


import os
import numpy as np
from mintpy.objects import timeseries
from mintpy.utils import ptime


############################################################
def temporal_filter(timeseries_file, outfile, time_win):

    # read timeseries info / data
    obj = timeseries(timeseries_file)
    obj.open()

    tbase = np.array(obj.yearList, np.float32).reshape(-1, 1)
    tbase -= tbase[obj.refIndex]

    ts_data = obj.read().reshape(obj.numDate, -1)

    # Smooth acquisitions / moving window in time one by one
    print('-'*50)
    print('filtering in time Gaussian window with size of {:.1f} years'.format(time_win))
    ts_data_filt = np.zeros(ts_data.shape, np.float32)
    prog_bar = ptime.progressBar(maxValue=obj.numDate)
    for i in range(obj.numDate):
        # Weight from Gaussian (normal) distribution in time
        tbase_diff = tbase[i] - tbase
        weight = np.exp(-0.5 * (tbase_diff**2) / (time_win**2))
        weight /= np.sum(weight)
        # Smooth the current acquisition
        ts_data_filt[i, :] = np.sum(ts_data * weight, axis=0)
        prog_bar.update(i+1, suffix=obj.dateList[i])
    prog_bar.close()
    del ts_data
    ts_data_filt -= ts_data_filt[obj.refIndex, :]
    ts_data_filt = np.reshape(ts_data_filt, (obj.numDate, obj.length, obj.width))

    # write filtered timeseries file
    if not outfile:
        outfile = '{}_tempGaussian.h5'.format(os.path.splitext(timeseries_file)[0])
    obj_out = timeseries(outfile)
    obj_out.write2hdf5(ts_data_filt, refFile=timeseries_file)
