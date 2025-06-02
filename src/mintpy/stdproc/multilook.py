############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Sep 2024                           #
############################################################


import logging
import os
import warnings

import numpy as np

# suppress numpy.RuntimeWarning message
np_logger = logging.getLogger('numpy')
np_logger.setLevel(logging.WARNING)


def multilook_data(data, lks_y=1, lks_x=1, method='mean'):
    """Apply multilooking (spatial averaging/resampling) to a multi-dimensional array.

    Link: https://stackoverflow.com/questions/34689519

    Parameters: data     - 2D / 3D np.array in real or complex
                lks_y    - int, number of multilook in y/azimuth direction
                lks_x    - int, number of multilook in x/range direction
                method   - str, multilook method, mean, median or nearest
    Returns:    out_data - 2D / 3D np.array after multilooking in last two dimension
    """
    # check - method
    method_list = ['mean', 'median', 'nearest']
    if method not in method_list:
        raise ValueError(f'Un-supported multilook method: {method}! Available methods: {method_list}.')

    # check - number of looks: do nothing if no multilook is applied
    lks_y = int(lks_y)
    lks_x = int(lks_x)
    if lks_y * lks_x == 1:
        return data

    shape = np.array(data.shape, dtype=float)
    if len(shape) == 2:
        # prepare - crop data to the exact multiple of the multilook number
        new_shape = np.floor(shape / (lks_y, lks_x)).astype(int) * (lks_y, lks_x)
        crop_data = data[:new_shape[0], :new_shape[1]]

        if method in ['mean', 'median']:
            # approach 1: reshape to higher dimensions, then collapse the extra dimensions
            #   with desired math operation
            # a. reshape to more dimensions
            temp = crop_data.reshape(
                (new_shape[0] // lks_y, lks_y,
                 new_shape[1] // lks_x, lks_x,
                ),
            )

            # b. collapse the extra dimensions with mean / median
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                if method == 'mean':
                    out_data = np.nanmean(temp, axis=(1, 3))
                elif method == 'median':
                    out_data = np.nanmedian(temp, axis=(1, 3))

            # approach 2: indexing + averaging: first in col/x, then in row/y
            # out_len, out_wid = np.floor(shape / (lks_y, lks_x)).astype(int)
            # out_data_col = np.zeros((shape[0], out_wid), dtype=data.dtype)
            # out_data = np.zeros((out_len, out_wid), dtype=data.dtype)

            # for x in range(out_wid):
            #     x0, x1 = x * lks_x, (x + 1) * lks_x
            #     out_data_col[:, x] = np.nanmean(crop_data[:, x0:x1], 1)
            # for y in range(out_len):
            #     y0, y1 = y * lks_y, (y + 1) * lks_y
            #     out_data[y, :] = np.nanmean(out_data_col[y0:y1, :], 0)

        elif method == 'nearest':
            out_data = crop_data[
                int(lks_y/2)::lks_y,
                int(lks_x/2)::lks_x,
            ]

    elif len(shape) == 3:
        # prepare - crop data to the exact multiple of the multilook number
        new_shape = np.floor(shape / (1, lks_y, lks_x)).astype(int) * (1, lks_y, lks_x)
        crop_data = data[
            :new_shape[0],
            :new_shape[1],
            :new_shape[2],
        ]

        if method in ['mean', 'median']:
            # a. reshape to more dimensions
            temp = crop_data.reshape(
                (new_shape[0],
                 new_shape[1] // lks_y, lks_y,
                 new_shape[2] // lks_x, lks_x,
                ),
            )

            # b. collapse the extra dimensions with mean / median
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                if method == 'mean':
                    out_data = np.nanmean(temp, axis=(2, 4))
                elif method == 'median':
                    out_data = np.nanmedian(temp, axis=(2, 4))

        elif method == 'nearest':
            out_data = crop_data[
                :,
                int(lks_y/2)::lks_y,
                int(lks_x/2)::lks_x,
            ]

    else:
        raise ValueError(f'Un-supported data dimension: {shape} --> {len(shape)}!')

    # ensure output data type
    out_data = np.array(out_data, dtype=data.dtype)

    return out_data


def multilook_gdal(in_file, lks_y, lks_x, box=None, out_file=None):
    """Apply multilooking via gdal_translate.

    Parameters: in_file  - str, path to the input GDAL VRT file
                lks_y/x  - int, number of looks in Y/X direction
                box      - tuple(int), area of interest in x0, y0, x1, y1
                out_file - str, path to the output data file
    Returns:    out_file - str, path to the output data file
    """
    print('apply multilooking via gdal.Translate ...')
    from osgeo import gdal

    # default output file name
    if not out_file:
        if in_file.endswith('.rdr.full.vrt'):
            out_file = in_file[:-9]
        elif in_file.endswith('.rdr.vrt'):
            out_file = in_file[:-4] + '.mli'
        else:
            raise ValueError(f'un-recognized ISCE VRT file ({in_file})!')
    print(f'input : {in_file}')
    print(f'output: {out_file}')

    ds = gdal.Open(in_file, gdal.GA_ReadOnly)
    if not box:
        box = (0, 0, ds.RasterXSize, ds.RasterYSize)

    in_wid = box[2] - box[0]
    in_len = box[3] - box[1]

    out_wid = int(in_wid / lks_x)
    out_len = int(in_len / lks_y)

    src_x0 = box[0]
    src_y0 = box[1]
    src_x1 = src_x0 + out_wid * lks_x
    src_y1 = src_y0 + out_len * lks_y

    # link: https://stackoverflow.com/questions/68025043
    options_str = f'-of ENVI -a_nodata 0 -outsize {out_wid} {out_len} '
    options_str += f' -srcwin {src_x0} {src_y0} {src_x1} {src_y1} '
    gdal.Translate(out_file, ds, options=options_str, callback=gdal.TermProgress_nocb)

    # generate GDAL .vrt file
    dso = gdal.Open(out_file, gdal.GA_ReadOnly)
    gdal.Translate(out_file+'.vrt', dso, options=gdal.TranslateOptions(format='VRT'))

    # generate ISCE .xml file
    if not os.path.isfile(out_file+'.xml'):
        from isce.applications.gdal2isce_xml import gdal2isce_xml
        gdal2isce_xml(out_file+'.vrt')

    return out_file
