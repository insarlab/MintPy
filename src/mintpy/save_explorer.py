############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Mahmud Haghighi, Mar 2025                        #
############################################################


import os

import numpy as np

from mintpy.objects import timeseries
from mintpy.utils import ptime, readfile
from mintpy.save_gmt import write_grd_file


####################################################################################
def convert2mm(data, atr):
    """Convert displacement (velocity) data to mm or mm/year."""
    if atr['UNIT'] in ['m', 'm/year']:
        return data * 1000
    elif atr['UNIT'] in ['cm', 'cm/year']:
        return data * 10
    elif atr['UNIT'] in ['mm', 'mm/year']:
        return data
    else:
        raise ValueError(f"ERROR: unit {atr['UNIT']} is not supported!")


def save_explorer(inps):

    # create output directory
    inps.outdir = os.path.abspath(inps.outdir)
    print(f'output directory: {inps.outdir}')
    if not os.path.isdir(inps.outdir):
        print('output directory does not exist, creating directory: '+inps.outdir)
        os.makedirs(inps.outdir)

    # grab aux files'name based on input ts_file file name
    inps.ts_file = os.path.abspath(inps.ts_file)
    ts_dir = os.path.dirname(inps.ts_file)
    prefix = 'geo_' if os.path.basename(inps.ts_file).startswith('geo_') else ''
    inps.vel_file = inps.vel_file or os.path.join(ts_dir, f'{prefix}velocity.h5')
    inps.msk_file = inps.msk_file or os.path.join(ts_dir, f'{prefix}maskTempCoh.h5')

    # check aux file existence
    inps.vel_file = os.path.abspath(inps.vel_file) if os.path.isfile(inps.vel_file) else None
    inps.msk_file = os.path.abspath(inps.msk_file) if os.path.isfile(inps.msk_file) else None
    print(f'time series file: {inps.ts_file}')
    print(f'velocity    file: {inps.vel_file}')
    print(f'mask        file: {inps.msk_file}')

    # read mask file
    if inps.msk_file and os.path.isfile(inps.msk_file):
        print(f'read mask data from file: {inps.msk_file}')
        mask = readfile.read(inps.msk_file)[0]
    else:
        mask = None

    # export velocity file
    if inps.vel_file and os.path.isfile(inps.vel_file):
        data, atr = readfile.read(inps.vel_file)
        data = convert2mm(data, atr)  # convert to mm
        if mask is not None:
            data[~mask] = np.nan

        out_base = os.path.splitext(os.path.basename(inps.vel_file))[0]
        out_file = os.path.join(inps.outdir, f'{out_base}_mm.grd')
        write_grd_file(data, atr, out_file, print_msg=True)

    # export time series to a list of grd files
    print('writing time series to a list of timeseries-{YYYYMMDD}_mm.grd files ...')
    date_list = timeseries(inps.ts_file).get_date_list()
    num_date = len(date_list)
    prog_bar = ptime.progressBar(maxValue=num_date)
    for i, date_str in enumerate(date_list):
        prog_bar.update(i+1, suffix=f'{i+1}/{num_date} {date_str}')

        # read
        data, atr = readfile.read(inps.ts_file, datasetName=date_str)
        data = convert2mm(data, atr)  # convert to mm
        if mask is not None:
            data[mask==0] = np.nan

        # write
        out_file = os.path.join(inps.outdir, f'timeseries-{date_str}_mm.grd')
        write_grd_file(data, atr, out_file, print_msg=False)
    prog_bar.close()

    print('Done.')
    return
