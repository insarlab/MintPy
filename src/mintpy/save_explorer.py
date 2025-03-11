############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Mahmud Haghighi, Mar 2025       #
############################################################


import os
from datetime import datetime

import numpy as np

from mintpy.save_gmt import write_grd_file
from mintpy.utils import readfile


####################################################################################
def convert2mm(data, atr):
    if atr['UNIT'] in ['m', 'm/year']:
        return data * 1000
    elif atr['UNIT'] in ['cm', 'cm/year']:
        return data * 10
    elif atr['UNIT'] in ['mm', 'mm/year']:
        return data
    else:
        raise Exception(f"ERROR: unit {atr['UNIT']} is not supported!")


def save_explorer(inps):

    if not os.path.exists(inps.outdir):
        print('output directory does not exist. creating directory: '+inps.outdir)
        os.makedirs(inps.outdir)

    if inps.mask_file:
        mask, atr = readfile.read(inps.mask_file)
    else:
        mask = None

    # export velocity file
    if inps.vel_file:
        data, atr = readfile.read(inps.vel_file)
        data = convert2mm(data, atr)

        if mask is not None:
            data[~mask] = np.nan

        out_file = os.path.join(inps.outdir, os.path.splitext(os.path.basename(inps.vel_file))[0] + '_mm.grd')
        write_grd_file(data, atr, out_file)


    # get slice list
    slice_list = readfile.get_slice_list(inps.file)

    for i, slice_name in enumerate(slice_list): # write each slice to a separate file
        if not slice_name.lower().startswith('timeseries'):
            continue

        data, atr = readfile.read(inps.file, datasetName=slice_name)
        # convert to mm
        data = convert2mm(data, atr)

        if mask is not None:
            data[~mask] = np.nan

        out_file = inps.outdir + '/' + slice_name + '_mm.grd'
        write_grd_file(data, atr, out_file)
        print(f'{i+1}/{len(slice_list)}: {out_file}')

    print('Done.')
    return
