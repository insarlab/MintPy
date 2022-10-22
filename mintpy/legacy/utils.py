############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, 2013                             #
############################################################


import time

import numpy as np


###################################################
def timeseries_inversion_FGLS(h5flat, h5timeseries):
    """Implementation of the SBAS algorithm.

    Usage:
    timeseries_inversion(h5flat,h5timeseries)
      h5flat: hdf5 file with the interferograms
      h5timeseries: hdf5 file with the output from the inversion
    ##################################################
    """

    total = time.time()
    A, B = design_matrix(h5flat)
    tbase, dateList, dateDict, dateDict2 = date_list(h5flat)
    dt = np.diff(tbase)
    B1 = np.array(np.linalg.pinv(B), dtype=np.float32)
    ifgram_list = list(h5flat['interferograms'].keys())
    ifgram_num = len(ifgram_list)
    #dset = h5flat[ifgram_list[0]].get(h5flat[ifgram_list[0]].keys()[0])
    #data = dset[0:dset.shape[0],0:dset.shape[1]]
    dset = h5flat['interferograms'][ifgram_list[0]].get(ifgram_list[0])
    data = dset[0:dset.shape[0], 0:dset.shape[1]]
    pixel_num = np.shape(data)[0]*np.shape(data)[1]
    print('Reading in the interferograms')
    # print ifgram_num,pixel_num
    print('number of interferograms: '+str(ifgram_num))
    print('number of pixels: '+str(pixel_num))
    pixel_num_step = int(pixel_num/10)

    data = np.zeros((ifgram_num, pixel_num), np.float32)
    for ni in range(ifgram_num):
        dset = h5flat['interferograms'][ifgram_list[ni]].get(ifgram_list[ni])
        #dset = h5flat[ifgram_list[ni]].get(h5flat[ifgram_list[ni]].keys()[0])
        d = dset[0:dset.shape[0], 0:dset.shape[1]]
        # print np.shape(d)

    del d
    dataPoint = np.zeros((ifgram_num, 1), np.float32)
    modelDimension = np.shape(B)[1]
    ts_data = np.zeros((date_num, pixel_num), np.float32)
    for ni in range(pixel_num):
        dataPoint = data[:, ni]
        nan_ndx = dataPoint == 0.
        fin_ndx = dataPoint != 0.
        nan_fin = dataPoint.copy()
        nan_fin[nan_ndx] = 1
        if not nan_fin.sum() == len(nan_fin):
            B1tmp = np.dot(B1, np.diag(fin_ndx))
            tmpe_ratea = np.dot(B1tmp, dataPoint)
            zero = np.array([0.], np.float32)
            defo = np.concatenate((zero, np.cumsum([tmpe_ratea*dt])))
            ts_data[:, ni] = defo
        if not np.remainder(ni, pixel_num_step):
            print('Processing point: %8d of %8d, %3d' %
                  (ni, pixel_num, (10*ni/pixel_num_step))+'%')
    del data
    timeseries = np.zeros((date_num, np.shape(dset)[0], np.shape(dset)[1]), np.float32)
    factor = -1*float(h5flat['interferograms'][ifgram_list[0]].attrs['WAVELENGTH'])/(4.*np.pi)
    for ni in range(date_num):
        timeseries[ni] = ts_data[ni].reshape(np.shape(dset)[1], np.shape(dset)[0]).T
        timeseries[ni] = timeseries[ni]*factor
    del ts_data
    timeseriesDict = {}
    for key, value in h5flat['interferograms'][ifgram_list[0]].attrs.items():
        timeseriesDict[key] = value

    dateIndex = {}
    for ni in range(len(dateList)):
        dateIndex[dateList[ni]] = ni
    if not 'timeseries' in h5timeseries:
        group = h5timeseries.create_group('timeseries')
        for key, value in iter(timeseriesDict.items()):
            group.attrs[key] = value

    for date in dateList:
        if not date in h5timeseries['timeseries']:
            dset = group.create_dataset(date, data=timeseries[dateIndex[date]])
    print('Time series inversion took ' + str(time.time()-total) + ' secs')
    return


def timeseries_inversion_L1(h5flat, h5timeseries):
    try:
        from cvxopt import matrix, normal

        from .l1 import l1
    except ImportError:
        raise ImportError('cvxopt should be installed to be able to use the L1 norm minimization.')
        # modified from sbas.py written by scott baker, 2012

    total = time.time()
    A, B = design_matrix(h5flat)
    tbase, dateList, dateDict, dateDict2 = date_list(h5flat)
    dt = np.diff(tbase)
    BL1 = matrix(B)
    B1 = np.array(np.linalg.pinv(B), dtype=np.float32)
    ifgram_list = list(h5flat['interferograms'].keys())
    ifgram_num = len(ifgram_list)
    #dset = h5flat[ifgram_list[0]].get(h5flat[ifgram_list[0]].keys()[0])
    #data = dset[0:dset.shape[0],0:dset.shape[1]]
    dset = h5flat['interferograms'][ifgram_list[0]].get(ifgram_list[0])
    data = dset[0:dset.shape[0], 0:dset.shape[1]]
    pixel_num = np.shape(data)[0]*np.shape(data)[1]
    print('Reading in the interferograms')
    print(ifgram_num, pixel_num)

    #data = np.zeros((ifgram_num,pixel_num),np.float32)
    data = np.zeros((ifgram_num, pixel_num))
    for ni in range(ifgram_num):
        dset = h5flat['interferograms'][ifgram_list[ni]].get(ifgram_list[ni])
        #dset = h5flat[ifgram_list[ni]].get(h5flat[ifgram_list[ni]].keys()[0])
        d = dset[0:dset.shape[0], 0:dset.shape[1]]
        # print np.shape(d)

        data[ni] = d.flatten(1)
    del d
    dataPoint = np.zeros((ifgram_num, 1), np.float32)
    modelDimension = np.shape(B)[1]
    ts_data = np.zeros((date_num, pixel_num), np.float32)
    print(data.shape)
    DataL1 = matrix(data)
    L1ORL2 = np.ones((pixel_num, 1))
    for ni in range(pixel_num):
        print(ni)
        dataPoint = data[:, ni]
        nan_ndx = dataPoint == 0.
        fin_ndx = dataPoint != 0.
        nan_fin = dataPoint.copy()
        nan_fin[nan_ndx] = 1
        if not nan_fin.sum() == len(nan_fin):

            B1tmp = np.dot(B1, np.diag(fin_ndx))
            #tmpe_ratea = np.dot(B1tmp,dataPoint)
            try:
                tmpe_ratea = np.array(l1(BL1, DataL1[:, ni]))
                zero = np.array([0.], np.float32)
                defo = np.concatenate((zero, np.cumsum([tmpe_ratea[:, 0]*dt])))
            except:
                tmpe_ratea = np.dot(B1tmp, dataPoint)
                L1ORL2[ni] = 0
                zero = np.array([0.], np.float32)
                defo = np.concatenate((zero, np.cumsum([tmpe_ratea*dt])))

            ts_data[:, ni] = defo
        if not np.remainder(ni, 10000):
            print('Processing point: %7d of %7d ' % (ni, pixel_num))
    del data
    timeseries = np.zeros((date_num, np.shape(dset)[0], np.shape(dset)[1]), np.float32)
    factor = -1*float(h5flat['interferograms'][ifgram_list[0]].attrs['WAVELENGTH'])/(4.*np.pi)
    for ni in range(date_num):
        timeseries[ni] = ts_data[ni].reshape(np.shape(dset)[1], np.shape(dset)[0]).T
        timeseries[ni] = timeseries[ni]*factor
    del ts_data
    L1ORL2 = np.reshape(L1ORL2, (np.shape(dset)[1], np.shape(dset)[0])).T

    timeseriesDict = {}
    for key, value in h5flat['interferograms'][ifgram_list[0]].attrs.items():
        timeseriesDict[key] = value

    dateIndex = {}
    for ni in range(len(dateList)):
        dateIndex[dateList[ni]] = ni
    if not 'timeseries' in h5timeseries:
        group = h5timeseries.create_group('timeseries')
        for key, value in iter(timeseriesDict.items()):
            group.attrs[key] = value

    for date in dateList:
        if not date in h5timeseries['timeseries']:
            dset = group.create_dataset(date, data=timeseries[dateIndex[date]])
    print('Time series inversion took ' + str(time.time()-total) + ' secs')
    L1orL2h5 = h5py.File('L1orL2.h5', 'w')
    gr = L1orL2h5.create_group('mask')
    dset = gr.create_dataset('mask', data=L1ORL2, compression='gzip')
    L1orL2h5.close()
    return
