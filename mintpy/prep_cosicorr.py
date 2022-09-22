############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Forrest Williams, Mar 2021                       #
############################################################


import datetime as dt
import os

import numpy as np

from mintpy.utils import readfile, utils1 as ut, writefile


#########################################################################
def add_cosicorr_metadata(fname, date12_dict, meta):
    '''Read/extract attribute data from cosicorr metadata file and add to metadata dictionary
    Parameters: fname       - str, offset or SNR file name
                date12_dict - dict, key for file name, value for date12 pairs
                meta        - dict, metadata dictionary
    Returns:    meta        - dict, metadata dictionary
    '''

    # add general attributes
    meta['PROCESSOR'] = 'cosicorr'
    meta['P_BASELINE_TOP_HDR'] = 0.0 #placeholder
    meta['P_BASELINE_BOTTOM_HDR'] = 0.0 #placeholder
    meta['RANGE_PIXEL_SIZE'] = np.abs(meta['X_STEP'])
    meta['AZIMUTH_PIXEL_SIZE'] = np.abs(meta['Y_STEP'])
    meta['RLOOKS'] = 1
    meta['ALOOKS'] = 1

    # Time attributes
    date1_str, date2_str = date12_dict[os.path.basename(fname)].split('-')
    meta['DATE12'] = f'{date1_str}-{date2_str}'
    date1 = dt.datetime.strptime(date1_str, '%Y%m%d')
    date2 = dt.datetime.strptime(date2_str, '%Y%m%d')
    date_avg = date1 + (date2 - date1) / 2
    date_avg_seconds = (date_avg - date_avg.replace(hour=0, minute=0, second=0, microsecond=0)).total_seconds()
    meta['CENTER_LINE_UTC'] = date_avg_seconds

    # add LAT/LON_REF1/2/3/4
    N = float(meta['Y_FIRST'])
    W = float(meta['X_FIRST'])
    S = N + float(meta['Y_STEP']) * int(meta['LENGTH'])
    E = W + float(meta['X_STEP']) * int(meta['WIDTH'])

    meta['LAT_REF1'] = str(S)
    meta['LAT_REF2'] = str(S)
    meta['LAT_REF3'] = str(N)
    meta['LAT_REF4'] = str(N)
    meta['LON_REF1'] = str(W)
    meta['LON_REF2'] = str(E)
    meta['LON_REF3'] = str(W)
    meta['LON_REF4'] = str(E)

    return(meta)


def prep_cosicorr(inps):
    """Prepare the COSI-Corr metadata."""

    # open and read hyp3 metadata
    date12_dict = {}
    with open(inps.meta_file) as f:
        for line in f:
            name, date1, date2 = line.strip().split(' ')
            date12_dict[name] = f'{date1}-{date2}'

    # loop over each filename
    inps.file = ut.get_file_list(inps.file, abspath=True)
    for fname in inps.file:
        # extra metadata
        meta = readfile.read_gdal_vrt(fname)
        meta = add_cosicorr_metadata(fname, date12_dict, meta)

        # write metadata into RSC file
        rsc_file = fname+'.rsc'
        writefile.write_roipac_rsc(meta, out_file=rsc_file)

    return
