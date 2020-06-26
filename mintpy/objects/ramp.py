############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2018                               #
############################################################
# Recommend import:
#   from mintpy.objects import deramp


import numpy as np

RAMP_LIST = [
    'linear', 
    'linear_range',
    'linear_azimuth',
    'quadratic',
    'quadratic_range',
    'quadratic_azimuth',
]


def deramp(data, mask_in, ramp_type='linear', metadata=None):
    '''Remove ramp from input data matrix based on pixel marked by mask
    Ignore data with nan or zero value.
    Parameters: data      : 2D / 3D np.ndarray, data to be derampped
                            If 3D, it's in size of (num_date, length, width)
                mask_in   : 2D np.ndarray, mask of pixels used for ramp estimation
                ramp_type : str, name of ramp to be estimated.
                metadata  : dict, containing reference pixel info, REF_Y/X
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

    ## mask
    # default
    if mask_in is None:
        mask_in = np.ones((length, width), dtype=np.float32)
    mask = (mask_in != 0).flatten()
    # ignore pixels with NaN or zero data value
    mask *= np.multiply(~np.isnan(dmean), dmean != 0.)

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
        raise ValueError('un-recognized ramp type: {}'.format(ramp_type))

    # estimate ramp
    X = np.dot(np.linalg.pinv(G[mask, :], rcond=1e-15), data[mask, :])
    ramp = np.dot(G, X)
    del X

    # reference in space if metadata
    if metadata and all(key in metadata.keys() for key in ['REF_X','REF_Y']):
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

