############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import logging
import os
import warnings

import h5py
import numpy as np

from mintpy.utils import attribute as attr, readfile, writefile

# suppress numpy.RuntimeWarning message
np_logger = logging.getLogger('numpy')
np_logger.setLevel(logging.WARNING)


######################################## Sub Functions ############################################
def multilook_matrix(matrix, lks_y, lks_x):
    """Obsolete multilooking functions"""
    rows, cols = matrix.shape
    lks_x = int(lks_x)
    lks_y = int(lks_y)
    if lks_x == 1 and lks_y == 1:
        return matrix

    rows_mli = int(np.floor(rows/lks_y))
    cols_mli = int(np.floor(cols/lks_x))
    #thr = np.floor(lks_x*lks_y/2)
    matrix_Cmli = np.zeros((rows,    cols_mli))
    matrix_mli = np.zeros((rows_mli, cols_mli))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for c in range(cols_mli):
            matrix_Cmli[:, c] = np.nanmean(matrix[:, (c)*lks_x:(c+1)*lks_x], 1)
        for r in range(rows_mli):
            matrix_mli[r, :] = np.nanmean(matrix_Cmli[(r)*lks_y:(r+1)*lks_y, :], 0)
    del matrix, matrix_Cmli
    return matrix_mli


def multilook_data(data, lks_y=1, lks_x=1, method='mean'):
    """Modified from Praveen on StackOverflow:

    link: https://stackoverflow.com/questions/34689519/how-to-coarser-the-2-d-array-data-resolution

    Parameters: data        : 2D / 3D np.array in real or complex
                lks_y       : int, number of multilook in y/azimuth direction
                lks_x       : int, number of multilook in x/range direction
                method      : str, multilook method, mean, median or nearest
    Returns:    coarse_data : 2D / 3D np.array after multilooking in last two dimension
    """
    method_list = ['mean', 'median', 'nearest']
    if method not in method_list:
        msg = f'un-supported multilook method: {method}. '
        msg += f'Available methods: {method_list}'
        raise ValueError(msg)

    # do nothing if no multilook is applied
    lks_y = int(lks_y)
    lks_x = int(lks_x)
    if lks_y * lks_x == 1:
        return data

    dtype = data.dtype
    shape = np.array(data.shape, dtype=float)
    ysize = int(shape[0] / lks_y)
    xsize = int(shape[1] / lks_x)

    if len(shape) == 2:
        if method in ['mean', 'median']:
            # crop data to the exact multiple of the multilook number
            new_shape = np.floor(shape / (lks_y, lks_x)).astype(int) * (lks_y, lks_x)
            crop_data = data[:new_shape[0],
                             :new_shape[1]]

            # reshape to more dimensions and collapse the extra dimensions with mean
            temp = crop_data.reshape((new_shape[0] // lks_y, lks_y,
                                      new_shape[1] // lks_x, lks_x))

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                if method == 'mean':
                    coarse_data = np.nanmean(temp, axis=(1, 3))
                elif method == 'median':
                    coarse_data = np.nanmedian(temp, axis=(1, 3))

        elif method == 'nearest':
            coarse_data = data[int(lks_y/2)::lks_y,
                               int(lks_x/2)::lks_x]

            # fix size discrepency from average method
            if coarse_data.shape != (ysize, xsize):
                coarse_data = coarse_data[:ysize, :xsize]

    elif len(shape) == 3:
        if method in ['mean', 'median']:
            # crop data to the exact multiple of the multilook number
            new_shape = np.floor(shape / (1, lks_y, lks_x)).astype(int) * (1, lks_y, lks_x)
            crop_data = data[:new_shape[0],
                             :new_shape[1],
                             :new_shape[2]]

            # reshape to more dimensions and collapse the extra dimensions with mean
            temp = crop_data.reshape((new_shape[0],
                                      new_shape[1] // lks_y, lks_y,
                                      new_shape[2] // lks_x, lks_x))

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                if method == 'mean':
                    coarse_data = np.nanmean(temp, axis=(2, 4))
                elif method == 'median':
                    coarse_data = np.nanmedian(temp, axis=(2, 4))

        elif method == 'nearest':
            coarse_data = data[:,
                               int(lks_y/2)::lks_y,
                               int(lks_x/2)::lks_x]

            # fix size discrepency from average method
            if coarse_data.shape[-2:] != (ysize, xsize):
                coarse_data = coarse_data[:, :ysize, :xsize]

    # ensure output data type
    coarse_data = np.array(coarse_data, dtype=dtype)

    return coarse_data


def multilook_file(infile, lks_y, lks_x, outfile=None, method='mean', margin=[0,0,0,0], max_memory=4):
    """ Multilook input file
    Parameters: infile - str, path of input file to be multilooked.
                lks_y  - int, number of looks in y / row direction.
                lks_x  - int, number of looks in x / column direction.
                margin - list of 4 int, number of pixels to be skipped during multilooking.
                         useful for offset product, where the marginal pixels are ignored during
                         cross correlation matching.
                outfile - str, path of output file
    Returns:    outfile - str, path of output file
    """
    lks_y = int(lks_y)
    lks_x = int(lks_x)

    # input file info
    atr = readfile.read_attribute(infile)
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    k = atr['FILE_TYPE']
    print('multilooking {} {} file: {}'.format(atr['PROCESSOR'], k, infile))
    print('number of looks in y / azimuth direction: %d' % lks_y)
    print('number of looks in x / range   direction: %d' % lks_x)
    print(f'multilook method: {method}')

    # margin --> box
    if margin is not [0,0,0,0]:    # top, bottom, left, right
        box = (margin[2], margin[0], width - margin[3], length - margin[1])
        print(f'number of pixels to skip in top/bottom/left/right boundaries: {margin}')
    else:
        box = (0, 0, width, length)

    # output file name
    ext = os.path.splitext(infile)[1]
    if not outfile:
        if os.getcwd() == os.path.dirname(os.path.abspath(infile)):
            suffix = f'_{lks_y}alks_{lks_x}rlks'
            outfile = os.path.splitext(infile)[0] + suffix + ext
        else:
            outfile = os.path.basename(infile)

    # update metadata
    atr = attr.update_attribute4multilook(atr, lks_y, lks_x, box=box)

    if ext in ['.h5', '.he5']:
        writefile.layout_hdf5(outfile, metadata=atr, ref_file=infile)

    # read source data and multilooking
    dsNames = readfile.get_dataset_list(infile)
    maxDigit = max(len(i) for i in dsNames)
    dsDict = dict()
    for dsName in dsNames:
        print('multilooking {d:<{w}} from {f} ...'.format(
            d=dsName, w=maxDigit, f=os.path.basename(infile)))

        # split in Y/row direction for IO for HDF5 only
        if ext in ['.h5', '.he5']:
            # calc step size with memory usage up to 4 GB
            with h5py.File(infile, 'r') as f:
                ds = f[dsName]
                ds_size = np.prod(ds.shape) * 4
            num_step = int(np.ceil(ds_size * 4 / (max_memory * 1024**3)))
            row_step = int(np.rint(length / num_step / 10) * 10)
            row_step = max(row_step, 10)

        else:
            row_step = box[3] - box[1]

        num_step = int(np.ceil((box[3] - box[1]) / (row_step * lks_y)))
        for i in range(num_step):
            r0 = box[1] + row_step * lks_y * i
            r1 = box[1] + row_step * lks_y * (i + 1)
            r1 = min(r1, box[3])
            # IO box
            box_i = (box[0], r0, box[2], r1)
            box_o = (int((box[0] - box[0]) / lks_x),
                     int((r0     - box[1]) / lks_y),
                     int((box[2] - box[0]) / lks_x),
                     int((r1     - box[1]) / lks_y))
            print(f'box: {box_o}')

            # read / multilook
            if method == 'nearest':
                data = readfile.read(infile,
                                     datasetName=dsName,
                                     box=box_i,
                                     xstep=lks_x,
                                     ystep=lks_y,
                                     print_msg=False)[0]

            else:
                data = readfile.read(infile,
                                     datasetName=dsName,
                                     box=box_i,
                                     print_msg=False)[0]

                data = multilook_data(data, lks_y, lks_x, method=method)

            # output block
            if data.ndim == 3:
                block = [0, data.shape[0],
                         box_o[1], box_o[3],
                         box_o[0], box_o[2]]
            else:
                block = [box_o[1], box_o[3],
                         box_o[0], box_o[2]]

            # write
            if ext in ['.h5', '.he5']:
                writefile.write_hdf5_block(outfile,
                                           data=data,
                                           datasetName=dsName,
                                           block=block,
                                           print_msg=False)
            else:
                dsDict[dsName] = data

    # for binary file with 2 bands, always use BIL scheme
    if (len(dsDict.keys()) == 2
            and os.path.splitext(infile)[1] not in ['.h5','.he5']
            and atr.get('INTERLEAVE', 'BIL').upper() != 'BIL'):
        print('the input binary file has 2 bands with band interleave as: {}'.format(atr['INTERLEAVE']))
        print('for the output binary file, change the band interleave to BIL as default.')
        atr['INTERLEAVE'] = 'BIL'

    if ext not in ['.h5', '.he5']:
        atr['BANDS'] = len(dsDict.keys())
        writefile.write(dsDict, out_file=outfile, metadata=atr, ref_file=infile)

        # write extra metadata files for ISCE data files
        if os.path.isfile(infile+'.xml') or os.path.isfile(infile+'.aux.xml'):
            writefile.write_isce_xml(atr, outfile)

    return outfile
