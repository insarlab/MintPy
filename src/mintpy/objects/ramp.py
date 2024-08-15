"""Utilities for ramps."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2018                               #
############################################################
# Recommend import:
#   from mintpy.objects import deramp


import numpy as np

# duplicated in mintpy.cli.remove_ramp
RAMP_LIST = [
    'linear',
    'linear_range',
    'linear_azimuth',
    'quadratic',
    'quadratic_range',
    'quadratic_azimuth',
]


def deramp(data, mask_in=None, ramp_type='linear', metadata=None, max_num_sample=1e6, coeff_file=None,
           ignore_zero_value=True):
    '''Remove ramp from input data matrix based on pixel marked by mask
    Ignore data with nan or zero value.
    Parameters: data       : 2D / 3D np.ndarray, data to be derampped
                             If 3D, it's in size of (num_date, length, width)
                mask_in    : 2D np.ndarray, mask of pixels used for ramp estimation
                ramp_type  : str, name of ramp to be estimated.
                metadata   : dict, containing reference pixel info, REF_Y/X
                max_num_sample : float, max number of pixel sample,
                             above which the uniform sampling is applied to reduce sample size
                coeff_file : str, path to the text file to save the estimated ramp coefficients
                ignore_zero_value : bool, ignore pixels with zero values, default is True
                             Recommend: True for phase data and False for offset data
    Returns:    data_out   : 2D / 3D np.ndarray, data after deramping
                ramp       : 2D / 3D np.ndarray, estimated ramp
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

    ## mask

    # 1. default
    if mask_in is None:
        mask_in = np.ones((length, width), dtype=np.float32)
    mask = (mask_in != 0).flatten()
    del mask_in

    # 2. ignore pixels with NaN and/or zero data value
    mask *= ~np.isnan(dmean)
    if ignore_zero_value:
        mask *= dmean != 0.
    del dmean

    # 3. for big dataset: uniformally sample the data for ramp estimation
    if max_num_sample and np.sum(mask) > max_num_sample:
        step = int(np.ceil(np.sqrt(np.sum(mask) / max_num_sample)))
        if step > 1:
            sample_flag = np.zeros((length, width), dtype=np.bool_)
            sample_flag[int(step/2)::step,
                        int(step/2)::step] = 1
            mask *= sample_flag.flatten()
            del sample_flag

    # design matrix
    xx, yy = np.meshgrid(np.arange(0, width),
                         np.arange(0, length))
    xx = np.array(xx, dtype=np.float32).reshape(-1, 1)
    yy = np.array(yy, dtype=np.float32).reshape(-1, 1)
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
        raise ValueError(f'un-recognized ramp type: {ramp_type}')

    # estimate ramp
    X = np.dot(np.linalg.pinv(G[mask, :], rcond=1e-15), data[mask, :])
    ramp = np.dot(G, X)
    ramp = np.array(ramp, dtype=data.dtype)

    # write estimated coefficient to text file
    if coeff_file is not None:
        with open(coeff_file, 'a') as f:
            for i in range(X.T.shape[0]):
                coeff_str = '    '.join([f'{float(c):16.6e}' for c in X.T[i,:]])
                f.write(f'{coeff_str}\n')

    # reference in space if metadata
    if metadata and all(key in metadata.keys() for key in ['REF_X','REF_Y']):
        ref_y, ref_x = int(metadata['REF_Y']), int(metadata['REF_X'])
        ref_idx = ref_y * width + ref_x
        ramp -= ramp[ref_idx, :]

    # do not change pixel with original zero value
    if ignore_zero_value:
        ramp[data == 0] = 0

    data_out = data - ramp
    if len(dshape) == 3:
        ramp = np.moveaxis(ramp, -1, 0)
        data_out = np.moveaxis(data_out, -1, 0)
    ramp = ramp.reshape(dshape)
    data_out = data_out.reshape(dshape)
    return data_out, ramp
