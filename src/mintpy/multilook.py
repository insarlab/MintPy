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

from mintpy.stdproc.multilook import multilook_data, multilook_gdal
from mintpy.utils import attribute as attr, readfile, writefile

# suppress numpy.RuntimeWarning message
np_logger = logging.getLogger('numpy')
np_logger.setLevel(logging.WARNING)


######################################## Sub Functions ############################################
def multilook_file(infile, lks_y, lks_x, outfile=None, method='mean', max_memory=4,
                   search_win=None, xcorr_win=None, margin=0, off_file=None):
    """ Multilook input file.

    Parameters: infile     - str, path of input file to be multilooked.
                lks_y      - int, number of looks in y / row direction.
                lks_x      - int, number of looks in x / column direction.
                outfile    - str, path of output file
                max_memory - float, maximum used memory in GB
                search_win - list(int), ampcor (half) search window in (width, length)
                xcorr_win  - list(int), ampcor cross-correlation window in (width, length)
                margin     - int, ampcor margin
                off_file   - str, path to the ampcor dense offsets file
    Returns:    outfile    - str, path of output file
    """
    lks_y = int(lks_y)
    lks_x = int(lks_x)

    # input file info
    if infile.endswith('.vrt'):
        atr = readfile.read_gdal_vrt(infile)
    else:
        atr = readfile.read_attribute(infile)
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    print(f'multilooking file: {infile}')
    print(f'multilook method: {method}')
    print(f'number of looks in y / x direction: {lks_y} / {lks_x}')

    # box
    if search_win:
        x0 = margin + search_win[0]
        y0 = margin + search_win[1]
        # Note: this is the formular for PyCuAmpcor/cuDenseOffsets.py
        # the one for Ampcor CPU is different from the GPU version
        out_wid = (width - 2*margin - 2*search_win[0] - xcorr_win[0]) // lks_x + 1
        out_len = (length - 2*margin - 2*search_win[1] - xcorr_win[1]) // lks_x + 1
        x1 = x0 + out_wid * lks_x
        y1 = y0 + out_len * lks_y
        box = (x0, y0, x1, y1)

    elif off_file:
        # use the existing ampcor dense offsets file as reference
        # since the formular for Ampcor CPU version is unknown
        # assuming symmetrical margins in top/bottom/left/right
        atr_off = readfile.read_attribute(off_file)
        out_wid, out_len = int(atr_off['WIDTH']), int(atr_off['LENGTH'])
        x0 = (width - out_wid * lks_x) // 2
        y0 = (length - out_len * lks_y) // 2
        x1 = x0 + out_wid * lks_x
        y1 = y0 + out_len * lks_y
        box = (x0, y0, x1, y1)
        if not 0 <= x0 < x1 <= width or not 0 <= y0 < y1 <= length:
            msg = f'Calculated bbox ({box}) exceeds input file coverage (width={width}, length={length})!'
            raise ValueError(msg)

    else:
        box = (0, 0, width, length)

    # use GDAL if input is VRT file
    if infile.endswith('.vrt'):
        outfile = multilook_gdal(infile, lks_y, lks_x, box=box, out_file=outfile)
        return outfile

    # output file name
    if not outfile:
        if os.getcwd() == os.path.dirname(os.path.abspath(infile)):
            fbase, fext = os.path.splitext(infile)
            outfile = f'{fbase}_{lks_y}alks_{lks_x}rlks{fext}'
        else:
            outfile = os.path.basename(infile)
    # use outfile for file extension, to support non-hdf5 infile and hdf5 outfile
    fext = os.path.splitext(outfile)[1]

    # update metadata
    atr = attr.update_attribute4multilook(atr, lks_y, lks_x, box=box)

    if fext in ['.h5', '.he5']:
        writefile.layout_hdf5(outfile, metadata=atr, ref_file=infile)

    # read source data and multilooking
    dsNames = readfile.get_dataset_list(infile)
    maxDigit = max(len(i) for i in dsNames)
    dsDict = dict()
    for dsName in dsNames:
        print('multilooking {d:<{w}} from {f} ...'.format(
            d=dsName, w=maxDigit, f=os.path.basename(infile)))

        # split in Y/row direction for IO for HDF5 only
        if fext in ['.h5', '.he5']:
            # calc step size with memory usage up to 4 GB
            # use outfile as h5py may be be able to handle infile (non-hdf5)
            with h5py.File(outfile, 'r') as f:
                ds_size = np.prod(f[dsName].shape) * 4
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
            box_o = (int((box[0] - box[0]) / lks_x), int((r0 - box[1]) / lks_y),
                     int((box[2] - box[0]) / lks_x), int((r1 - box[1]) / lks_y))
            print(f'box: {box_o}')

            # read / multilook
            kwargs = dict(datasetName=dsName, box=box_i, print_msg=False)
            if method == 'nearest':
                data = readfile.read(infile, xstep=lks_x, ystep=lks_y, **kwargs)[0]
            else:
                data = readfile.read(infile, **kwargs)[0]
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
            if fext in ['.h5', '.he5']:
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
        print(f'the input binary file has 2 bands with band interleave as: {atr["INTERLEAVE"]}')
        print('for the output binary file, change the band interleave to BIL as default.')
        atr['INTERLEAVE'] = 'BIL'

    if fext not in ['.h5', '.he5']:
        atr['BANDS'] = len(dsDict.keys())
        writefile.write(dsDict, out_file=outfile, metadata=atr, ref_file=infile)

        # write extra metadata files for ISCE data files
        if os.path.isfile(infile+'.xml') or os.path.isfile(infile+'.aux.xml'):
            writefile.write_isce_xml(atr, outfile)

    return outfile
