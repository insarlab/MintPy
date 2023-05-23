############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os
import shutil

import numpy as np

from mintpy.utils import readfile, writefile


############################################################
def mask_matrix(data, mask, fill_value=np.nan):
    """mask a 2D matrxi data with mask
    Parameters: data : 2D / 3D np.array
                mask : 2D np.array of bool
                fill_value : number
    Returns:    data : same shape of array as input
    """
    # check data type
    if np.isnan(fill_value) and data.dtype.kind not in ('f','d','c'):
        msg = 'in order to fill the invalid pixels with np.nan'
        msg += f'\n\tconvert input matrix from {data.dtype} to np.float32'
        print(msg)
        data = np.array(data, dtype=np.float32)

    if len(data.shape) == 2:
        data[mask == 0] = fill_value
    elif len(data.shape) == 3:
        data[:, mask == 0] = fill_value
    return data


def update_mask_with_inps(mask, inps, print_msg=True):
    """Update mask matrix from input options: subset_x/y and mask_vmin/vmax."""
    vprint = print if print_msg else lambda *args, **kwargs: None

    if inps.subset_x:
        mask[:, 0:inps.subset_x[0]] = 0
        mask[:, inps.subset_x[1]:] = 0
        vprint(f'mask out area not in x: {inps.subset_x}')

    if inps.subset_y:
        mask[0:inps.subset_y[0], :] = 0
        mask[inps.subset_y[1]:, :] = 0
        vprint(f'mask out area not in y: {inps.subset_y}')

    if inps.mask_vmin:
        mask = mask >= inps.mask_vmin
        vprint(f'mask out pixels with value < {inps.mask_vmin} in mask file')

    if inps.mask_vmax:
        mask = mask <= inps.mask_vmax
        vprint(f'mask out pixels with value > {inps.mask_vmax} in mask file')

    return mask


def mask_file(fname, mask_file, out_file=None, fill_value=np.nan, inps=None):
    """ Mask input fname with mask_file
    Parameters: fname     - str, file to be masked
                mask_file - str, mask file
                out_file  - str, output file name
                inps      - namespace object, including:
                            subset_x/y
                            mask_vmin/vmax
    Returns:    out_file  - str, output file name
    """

    # read mask_file
    mask = readfile.read(mask_file)[0]
    if inps is not None:
        mask = update_mask_with_inps(mask, inps)

    # masking input file
    dsNames = readfile.get_dataset_list(fname)
    maxDigit = max(len(i) for i in dsNames)
    dsDict = {}
    for dsName in dsNames:
        print('masking {d:<{w}} from {f} ...'.format(d=dsName, w=maxDigit, f=fname))
        data = readfile.read(fname, datasetName=dsName, print_msg=False)[0]
        data = mask_matrix(data, mask, fill_value=fill_value)
        dsDict[dsName] = data

    # default output filename
    if not out_file:
        fbase, fext = os.path.splitext(fname)
        out_file = f'{fbase}_msk{fext}'

    writefile.write(dsDict, out_file=out_file, ref_file=fname)
    return out_file


def mask_isce_file(in_file, mask_file, out_file=None):
    if not in_file:
        return

    # read mask_file
    print(f'read mask from {mask_file}')
    mask = readfile.read(mask_file)[0]

    # mask isce file
    atr = readfile.read_attribute(in_file)
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    interleave = atr['INTERLEAVE'].upper()
    num_band = int(atr['BANDS'])

    # default short name for data type from ISCE
    data_type_dict = {
        'byte': 'bool_',
        'float': 'float32',
        'double': 'float64',
        'cfloat': 'complex64',
    }
    data_type = atr['DATA_TYPE'].lower()
    if data_type in data_type_dict.keys():
        data_type = data_type_dict[data_type]

    print(f'read {in_file}')
    print('setting the (phase) value on the masked out pixels to zero')
    fbase, ext = os.path.splitext(in_file)
    if ext == '.unw':
        kwargs = dict(
            shape=(length, width),
            data_type=data_type,
            num_band=num_band,
            interleave=interleave,
        )
        amp = readfile.read_binary(in_file, band=1, **kwargs)
        pha = readfile.read_binary(in_file, band=2, **kwargs)
        pha[mask == 0] = 0
        data = np.hstack((amp, pha)).flatten()
    elif ext == '.int':
        data = np.fromfile(in_file, dtype=data_type, count=length*width).reshape(-1, width)
        data[mask == 0] = 0
    elif ext in ['.cor','.conncomp']:
        data = readfile.read(in_file)[0]
        data[mask == 0] = 0
    else:
        raise ValueError(f'unsupported ISCE file: {in_file}')

    # output filename
    if not out_file:
        if ext in ['.int', '.cor', '.unw']:
            out_file = f'{fbase}_msk{ext}'
        elif in_file.endswith('.unw.conncomp'):
            out_file = '{}_msk.unw.conncomp'.format(in_file.split('.unw.conncomp')[0])
        else:
            raise ValueError(f'unrecognized input file type: {in_file}')

    data.tofile(out_file)
    print(f'finished writing to file {out_file}')

    # prepare ISCE metadata file by
    # 1. copy and rename metadata files
    # 2. update file path inside files
    for ext in ['xml', 'vrt']:
        # copy
        in_meta_file = f'{in_file}.{ext}'
        out_meta_file = f'{out_file}.{ext}'
        shutil.copy2(in_meta_file, out_meta_file)
        print(f'copy {in_meta_file} to {out_meta_file}')
        print('   and update the corresponding filename')

        # update file path
        meta_file = f'{out_file}.{ext}'
        with open(meta_file) as f:
            s = f.read()
        s = s.replace(os.path.basename(in_file),
                      os.path.basename(out_file))
        with open(meta_file, 'w') as f:
            f.write(s)
    return out_file
