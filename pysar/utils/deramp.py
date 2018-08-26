#####################################################################
# Program is part of PySAR                                          #
# Copyright(c) 2011-2018, Zhang Yunjun, Heresh Fattahi, Scott Baker #
# Author:  Zhang Yunjun, Heresh Fattahi, Scott Baker                #
#####################################################################
# Recommend import:
#     from pysar.utils import deramp


import os
import time
import h5py
import numpy as np
from scipy import linalg
from pysar.utils import ptime, readfile, writefile
from pysar.objects import timeseries, ifgramStack


##################################################################
def deramp_data(data, mask_in, ramp_type='linear', metadata=None):
    '''Remove ramp from input data matrix based on pixel marked by mask
    Ignore data with nan or zero value.
    Parameters: data      : 2D / 3D np.ndarray, data to be derampped
                mask_in   : 2D np.ndarray, mask of pixels used for ramp estimation
                metadata  : dict, containing reference pixel info, REF_Y/X
                ramp_type : str, name of ramp to be estimated.
    Returns:    data_out  : 2D / 3D np.ndarray, data after deramping
                ramp      : 2D / 3D np.ndarray, estimated ramp
    '''
    dshape = data.shape
    length, width = dshape[-2:]
    num_pixel = length * width

    # prepare input data
    if len(dshape) == 3:
        data = np.moveaxis(data, 0, -1)        #reshape to (length, width, numDate)
        data = data.reshape(num_pixel, -1)
        dmean = np.mean(data, axis=-1).flatten()
    else:
        data = data.reshape(-1, 1)
        dmean = np.array(data).flatten()

    # default mask
    if mask_in is None:
        mask_in = np.ones((length, width), dtype=np.float32)

    # design matrix
    xx, yy = np.meshgrid(np.arange(0, width),
                         np.arange(0, length))
    xx = xx.reshape(-1, 1)
    yy = yy.reshape(-1, 1)
    ones = np.ones(xx.shape, dtype=np.float32)
    if ramp_type == 'linear':
        G = np.hstack((yy, xx, ones))
    elif ramp_type == 'quadratic':
        G = np.hstack((yy**2, xx**2, yy*xx, yy, xx, ones))
    elif ramp_type == 'linear_range':
        G = np.hstack((xx, ones))
    elif ramp_type == 'linear_azimuth':
        G = np.hstack((yy, ones))
    elif ramp_type == 'quadratic_range':
        G = np.hstack((xx**2, xx, ones))
    elif ramp_type == 'quadratic_azimuth':
        G = np.hstack((yy**2, yy, ones))
    else:
        raise ValueError('un-recognized ramp type: {}'.format(ramp_type))

    # ignore pixels with NaN or zero data value
    mask = (mask_in != 0).flatten()
    mask[np.isnan(dmean)] = 0
    mask[dmean == 0] = 0

    # estimate ramp
    X = linalg.lstsq(G[mask, :], data[mask, :], cond=1e-8)[0]
    ramp = np.dot(G, X)

    # reference in space if metadata
    if metadata:
        ref_y, ref_x = int(metadata['REF_Y']), int(metadata['REF_X'])
        ref_idx = ref_y * width + ref_x
        ramp -= ramp[ref_idx, :]

    # do not change pixel with original zero value
    ramp[data == 0] = 0
    ramp = np.array(ramp, dtype=data.dtype)

    data_out = data - ramp
    if len(dshape) == 3:
        ramp = np.moveaxis(ramp, -1, 0)
        data_out = np.moveaxis(data_out, -1, 0)
    ramp = ramp.reshape(dshape)
    data_out = data_out.reshape(dshape)
    return data_out, ramp


##################################################################
def remove_data_multiple_surface(data, mask, ramp_type, ysub):
    ## ysub = [0,2400,2000,6800]
    dataOut = np.zeros(data.shape, data.dtype)
    dataOut[:] = np.nan

    surfaceNum = len(ysub)/2
    # 1st mask
    print('removing 1st surface ...')
    i = 0
    mask_i = np.zeros(data.shape, data.dtype)
    mask_i[ysub[2*i]:ysub[2*i+1], :] = mask[ysub[2*i]:ysub[2*i+1], :]

    dataOut_i, ramp_i = remove_data_surface(data, mask_i, ramp_type)
    dataOut[ysub[2*i]:ysub[2*i+1], :] = dataOut_i[ysub[2*i]:ysub[2*i+1], :]

    # 2 - last masks
    for i in range(1, surfaceNum):
        print('removing '+str(i+1)+'th surface ...')
        mask_i = np.zeros(data.shape, data.dtype)
        mask_i[ysub[2*i]:ysub[2*i+1], :] = mask[ysub[2*i]:ysub[2*i+1], :]

        dataOut_i, ramp_i = remove_data_surface(data, mask_i, ramp_type)

        if ysub[2*i] < ysub[2*i-1]:
            dataOut[ysub[2*i]:ysub[2*i-1], :] += dataOut_i[ysub[2*i]:ysub[2*i-1], :]
            dataOut[ysub[2*i]:ysub[2*i-1], :] /= 2
            dataOut[ysub[2*i-1]:ysub[2*i+1], :] = dataOut_i[ysub[2*i-1]:ysub[2*i+1], :]
        else:
            dataOut[ysub[2*i]:ysub[2*i+1], :] = dataOut_i[ysub[2*i]:ysub[2*i+1], :]

    return dataOut


##################################################################
def deramp_file(fname, ramp_type, mask_file=None, out_file=None, datasetName=None):
    print('remove {} ramp from file: {}'.format(ramp_type, fname))
    if not out_file:
        fbase, fext = os.path.splitext(fname)
        out_file = '{}_ramp{}'.format(fbase, fext)

    start_time = time.time()
    atr = readfile.read_attribute(fname)

    # mask
    if os.path.isfile(mask_file):
        mask = readfile.read(mask_file, datasetName='mask')[0]
        print('read mask file: '+mask_file)
    else:
        mask = np.ones((int(atr['LENGTH']), int(atr['WIDTH'])))
        print('use mask of the whole area')

    # deramping
    k = atr['FILE_TYPE']
    if k == 'timeseries':
        print('reading data ...')
        data = readfile.read(fname)[0]
        print('estimating phase ramp ...')
        data = deramp_data(data, mask, ramp_type=ramp_type, metadata=atr)[0]
        writefile.write(data, out_file, ref_file=fname)

    elif k == 'ifgramStack':
        obj = ifgramStack(fname)
        obj.open(print_msg=False)
        if not datasetName:
            datasetName = 'unwrapPhase'
        with h5py.File(fname, 'a') as f:
            ds = f[datasetName]
            dsNameOut = '{}_{}'.format(datasetName, ramp_type)
            if dsNameOut in f.keys():
                dsOut = f[dsNameOut]
                print('access HDF5 dataset /{}'.format(dsNameOut))
            else:
                dsOut = f.create_dataset(dsNameOut, shape=(obj.numIfgram, obj.length, obj.width),
                                         dtype=np.float32, chunks=True, compression=None)
                print('create HDF5 dataset /{}'.format(dsNameOut))

            prog_bar = ptime.progressBar(maxValue=obj.numIfgram)
            for i in range(obj.numIfgram):
                data = ds[i, :, :]
                data = deramp_data(data, mask, ramp_type=ramp_type, metadata=atr)[0]
                dsOut[i, :, :] = data
                prog_bar.update(i+1, suffix='{}/{}'.format(i+1, obj.numIfgram))
            prog_bar.close()
            print('finished writing to file: '.format(fname))

    # Single Dataset File
    else:
        data = readfile.read(fname)[0]
        data = deramp_data(data, mask, ramp_type, metadata=atr)[0]
        print('writing >>> {}'.format(out_file))
        writefile.write(data, out_file=out_file, ref_file=fname)

    m, s = divmod(time.time()-start_time, 60)
    print('\ntime used: {:02.0f} mins {:02.1f} secs'.format(m, s))
    return out_file
