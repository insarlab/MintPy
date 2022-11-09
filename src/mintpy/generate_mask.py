############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os
import time

import h5py
import numpy as np

from mintpy.utils import readfile, utils as ut, writefile


################################################################################################
def run_or_skip(inps):
    print('-'*50)
    print('update mode: ON')
    flag = 'skip'

    # check output file vs input dataset
    if not inps.outfile or not os.path.isfile(inps.outfile):
        flag = 'run'
        print(f'1) output file {inps.outfile} NOT exist.')
    else:
        print(f'1) output file {inps.outfile} already exists.')
        with h5py.File(inps.file, 'r') as f:
            ti = float(f[inps.dset].attrs.get('MODIFICATION_TIME', os.path.getmtime(inps.file)))
        to = os.path.getmtime(inps.outfile)
        if ti > to:
            flag = 'run'
            print(f'2) output file is NOT newer than input dataset: {inps.dset}.')
        else:
            print(f'2) output file is newer than input dataset: {inps.dset}.')

    # result
    print(f'run or skip: {flag}.')
    return flag


################################################################################################
def create_threshold_mask(inps):
    if inps.dset:
        print(f'read {inps.file} {inps.dset}')
    else:
        print('read %s' % (inps.file))
    data, atr = readfile.read(inps.file, datasetName=inps.dset)
    if len(data.shape) > 2:
        raise Exception('Only 2D dataset is supported for threshold method, input is 3D')
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    nanmask = ~np.isnan(data)

    # row/column range for threshold operations
    vmask = np.ones((length, width), dtype=np.bool_)
    if inps.v_subset_x:
        [vx0, vx1] = sorted(inps.v_subset_x)
        vmask[:, :vx0] = 0
        vmask[:, vx1:] = 0
    if inps.v_subset_y:
        [vy0, vy1] = sorted(inps.v_subset_y)
        vmask[:vy0, :] = 0
        vmask[vy1:, :] = 0
    if inps.vroipoly:
        from mintpy.utils import plot_ext
        poly_mask = plot_ext.get_poly_mask(inps.file, datasetName=inps.dset, view_cmd=inps.view_cmd)
        if poly_mask is not None:
            vmask[poly_mask == 0] = 0

    print('create initial mask with the same size as the input file and all = 1')
    mask = np.ones((length, width), dtype=np.bool_)

    # nan value
    if not inps.keep_nan:
        mask *= nanmask
        print('all pixels with nan value = 0')

    if inps.nonzero:
        mask[data == 0.] = 0
        print('exclude pixels with zero value')

    # min threshold
    if inps.vmin is not None:
        mask[vmask] *= ~(data[vmask] < inps.vmin)
        print('exclude pixels with value < %s' % str(inps.vmin))

    # max threshold
    if inps.vmax is not None:
        mask[vmask] *= ~(data[vmask] > inps.vmax)
        print('exclude pixels with value > %s' % str(inps.vmax))

    # remove small pixel clusters
    if inps.minpixels is not None:
        from skimage.morphology import remove_small_objects
        num_pixel = np.sum(mask)
        mask = remove_small_objects(mask, inps.minpixels, connectivity=1)
        print('exclude pixel clusters with size < %d pixels: remove %d pixels' % (inps.minpixels, num_pixel-np.sum(mask)))

    # remove pixels with large velocity STD
    if inps.vstd:
        if atr['FILE_TYPE'] != 'velocity':
            raise ValueError('Input file MUST be a velocity file when using the --vstd option!')
        data_std = readfile.read(inps.file, datasetName='velocityStd')[0]
        mask[nanmask] *= (np.abs(data[nanmask]) > (inps.vstd_num * data_std[nanmask]))
        print(f'exclude pixels according to the formula: |velocity| > {inps.vstd_num} * velocityStd')

    # subset in Y
    if inps.subset_y is not None:
        y0, y1 = sorted(inps.subset_y)
        mask[0:y0, :] = 0
        mask[y1:length, :] = 0
        print('exclude pixels with y OUT of [%d, %d]' % (y0, y1))

    # subset in x
    if inps.subset_x is not None:
        x0, x1 = sorted(inps.subset_x)
        mask[:, 0:x0] = 0
        mask[:, x1:width] = 0
        print('exclude pixels with x OUT of [%d, %d]' % (x0, x1))

    # exclude circular area
    if inps.ex_circle:
        x, y, r = inps.ex_circle
        cmask = ut.get_circular_mask(x, y, r, (length, width))
        mask[cmask == 1] = 0
        print(f'exclude pixels inside of circle defined as (x={x}, y={y}, r={r})')

    # include circular area
    if inps.in_circle:
        x, y, r = inps.in_circle
        cmask = ut.get_circular_mask(x, y, r, (length, width))
        mask[cmask == 0] = 0
        print(f'exclude pixels outside of circle defined as (x={x}, y={y}, r={r})')

    # interactively select polygonal region of interest (ROI)
    if inps.roipoly:
        from mintpy.utils import plot_ext
        poly_mask = plot_ext.get_poly_mask(inps.file, datasetName=inps.dset, view_cmd=inps.view_cmd)
        if poly_mask is not None:
            mask *= poly_mask

    # base mask
    if inps.base_file:
        # read base mask file
        base_data = readfile.read(inps.base_file, datasetName=inps.base_dataset)[0]
        if len(base_data.shape) == 3:
            base_data = np.sum(base_data, axis=0)

        # apply base mask
        mask[base_data == float(inps.base_value)] = 0

        # message
        msg = f'exclude pixels in base file {os.path.basename(inps.base_file)} '
        if inps.base_dataset:
            msg += f'dataset {inps.base_dataset} '
        msg += f'with value == {inps.base_value}'
        print(msg)

    # revert
    if inps.revert:
        temp = np.array(mask, dtype=np.bool_)
        mask[temp == True] = False
        mask[temp == False] = True
        del temp

    # Write mask file
    atr['FILE_TYPE'] = 'mask'
    writefile.write(mask, out_file=inps.outfile, metadata=atr)

    return inps.outfile


def create_mask(inps):
    """Create mask based on non-zero values or threshold."""

    start_time = time.time()
    ftype = readfile.read_attribute(inps.file)['FILE_TYPE']
    print(f'input {ftype} file: {inps.file}')

    # create mask using non-zero
    if inps.nonzero and ftype == 'ifgramStack':
        # update mode
        if inps.update_mode and run_or_skip(inps) == 'skip':
            return

        # run
        inps.outfile = ut.nonzero_mask(
            inps.file,
            out_file=inps.outfile,
            datasetName=inps.dset,
        )

    # create mask using threshold
    else:
        create_threshold_mask(inps)

    # used time
    m, s = divmod(time.time()-start_time, 60)
    print(f'time used: {m:02.0f} mins {s:02.1f} secs.')

    return
