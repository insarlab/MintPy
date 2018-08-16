#####################################################################
# Program is part of PySAR                                          #
# Copyright(c) 2011-2018, Heresh Fattahi, Zhang Yunjun, Scott Baker #
# Author:  Heresh Fattahi, Zhang Yunjun, Scott Baker                #
#####################################################################
# Recommend import:
#     from pysar.utils import deramp


import os
import time
import h5py
import numpy as np
from pysar.utils import ptime, readfile, writefile
from pysar.objects import timeseries, ifgramStack


##################################################################
def remove_data_surface(data, mask, surf_type='linear'):
    '''Remove surface from input data matrix based on pixel marked by mask'''
    mask[np.isnan(data)] = 0
    mask = mask.flatten(1)
    z = data.flatten(1)
    ndx = mask != 0
    x = list(range(0, np.shape(data)[1]))
    y = list(range(0, np.shape(data)[0]))
    x1, y1 = np.meshgrid(x, y)
    points = np.vstack((y1.flatten(1), x1.flatten(1))).T
    if surf_type == 'quadratic':
        G = np.array([points[:, 0]**2, points[:, 1]**2, points[:, 0], points[:, 1], points[:, 0]*points[:, 1],
                      np.ones(np.shape(points)[0])], np.float32).T
    elif surf_type == 'linear':
        G = np.array([points[:, 0], points[:, 1],
                      np.ones(np.shape(points)[0])], np.float32).T
    elif surf_type == 'quadratic_range':
        G = np.array([points[:, 1]**2, points[:, 1],
                      np.ones(np.shape(points)[0])], np.float32).T
    elif surf_type == 'quadratic_azimuth':
        G = np.array([points[:, 0]**2, points[:, 0],
                      np.ones(np.shape(points)[0])], np.float32).T
    elif surf_type == 'linear_range':
        G = np.array([points[:, 1],
                      np.ones(np.shape(points)[0])], np.float32).T
    elif surf_type == 'linear_azimuth':
        G = np.array([points[:, 0],
                      np.ones(np.shape(points)[0])], np.float32).T

    z = z[ndx]
    G = G[ndx]
    G1 = np.linalg.pinv(G)
    plane = np.dot(G1, z)

    if surf_type == 'quadratic':
        zplane = plane[0]*y1**2 + plane[1]*x1**2 + plane[2]*y1 + plane[3]*x1 + plane[4]*y1*x1 + plane[5]
    elif surf_type == 'linear':
        zplane = plane[0]*y1 + plane[1]*x1 + plane[2]
    elif surf_type == 'quadratic_range':
        zplane = plane[0]*x1**2 + plane[1]*x1 + plane[2]
    elif surf_type == 'quadratic_azimuth':
        zplane = plane[0]*y1**2 + plane[1]*y1 + plane[2]
    elif surf_type == 'linear_range':
        zplane = plane[0]*x1 + plane[1]
    elif surf_type == 'linear_azimuth':
        zplane = plane[0]*y1 + plane[1]

    data_n = data - zplane
    data_n[data == 0.] = 0.       # Do not change zero phase value
    data_n = np.array(data_n, data.dtype)
    zplane = np.array(zplane, data.dtype)

    return data_n, zplane


##################################################################
def remove_data_multiple_surface(data, mask, surf_type, ysub):
    ## ysub = [0,2400,2000,6800]
    dataOut = np.zeros(data.shape, data.dtype)
    dataOut[:] = np.nan

    surfaceNum = len(ysub)/2
    # 1st mask
    print('removing 1st surface ...')
    i = 0
    mask_i = np.zeros(data.shape, data.dtype)
    mask_i[ysub[2*i]:ysub[2*i+1], :] = mask[ysub[2*i]:ysub[2*i+1], :]

    dataOut_i, ramp_i = remove_data_surface(data, mask_i, surf_type)
    dataOut[ysub[2*i]:ysub[2*i+1], :] = dataOut_i[ysub[2*i]:ysub[2*i+1], :]

    # 2 - last masks
    for i in range(1, surfaceNum):
        print('removing '+str(i+1)+'th surface ...')
        mask_i = np.zeros(data.shape, data.dtype)
        mask_i[ysub[2*i]:ysub[2*i+1], :] = mask[ysub[2*i]:ysub[2*i+1], :]

        dataOut_i, ramp_i = remove_data_surface(data, mask_i, surf_type)

        if ysub[2*i] < ysub[2*i-1]:
            dataOut[ysub[2*i]:ysub[2*i-1], :] += dataOut_i[ysub[2*i]:ysub[2*i-1], :]
            dataOut[ysub[2*i]:ysub[2*i-1], :] /= 2
            dataOut[ysub[2*i-1]:ysub[2*i+1], :] = dataOut_i[ysub[2*i-1]:ysub[2*i+1], :]
        else:
            dataOut[ysub[2*i]:ysub[2*i+1], :] = dataOut_i[ysub[2*i]:ysub[2*i+1], :]

    return dataOut


##################################################################
def remove_surface(fname, surf_type, mask_file=None, out_file=None, ysub=None):
    start_time = time.time()
    atr = readfile.read_attribute(fname)
    k = atr['FILE_TYPE']
    print('remove {} ramp from file: '.format(surf_type, fname))

    if not out_file:
        out_file = '{base}_ramp{ext}'.format(base=os.path.splitext(fname)[0], ramp=surf_type,
                                               ext=os.path.splitext(fname)[1])

    if os.path.isfile(mask_file):
        mask = readfile.read(mask_file, datasetName='mask')[0]
        print('read mask file: '+mask_file)
    else:
        mask = np.ones((int(atr['LENGTH']), int(atr['WIDTH'])))
        print('use mask of the whole area')

    if k == 'timeseries':
        obj = timeseries(fname)
        data = obj.read()
        numDate = data.shape[0]
        print('estimating phase ramp ...')
        prog_bar = ptime.progressBar(maxValue=numDate)
        for i in range(numDate):
            if not ysub:
                data[i, :, :] = remove_data_surface(np.squeeze(data[i, :, :]),
                                                    mask, surf_type)[0]
            else:
                data[i, :, :] = remove_data_multiple_surface(np.squeeze(data),
                                                             mask, surf_type, ysub)
            prog_bar.update(i+1, suffix=obj.dateList[i])
        prog_bar.close()
        objOut = timeseries(out_file)
        objOut.write2hdf5(data=data, refFile=fname)

    elif k == 'ifgramStack':
        obj = ifgramStack(fname)
        obj.open(print_msg=False)
        with h5py.File(fname, 'a') as f:
            ds = f['unwrapPhase']

            dsName = 'unwrapPhase_{}'.format(surf_type)
            if dsName in f.keys():
                dsOut = f[dsName]
                print('access HDF5 dataset /{}'.format(dsName))
            else:
                dsOut = f.create_dataset(dsName, shape=(obj.numIfgram, obj.length, obj.width),
                                      dtype=np.float32, chunks=True, compression=None)
                print('create HDF5 dataset /{}'.format(dsName))

            prog_bar = ptime.progressBar(maxValue=obj.numIfgram)
            for i in range(obj.numIfgram):
                data = ds[i, :, :]
                mask_n = np.array(mask, np.bool_)
                mask_n[data == 0.] = 0

                if not ysub:
                    data_n = remove_data_surface(data, mask, surf_type)[0]
                else:
                    data_n = remove_data_multiple_surface(data, mask, surf_type, ysub)
                dsOut[i, :, :] = data_n
                prog_bar.update(i+1, suffix='{}/{}'.format(i+1, obj.numIfgram))
            prog_bar.close()

    # Single Dataset File
    else:
        data, atr = readfile.read(fname)
        data_n, ramp = remove_data_surface(data, mask, surf_type)
        writefile.write(data_n, out_file=out_file, metadata=atr)

    m, s = divmod(time.time()-start_time, 60)
    print('\ntime used: {:02.0f} mins {:02.1f} secs'.format(m, s))
    return out_file
