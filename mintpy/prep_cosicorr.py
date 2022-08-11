############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Forrest Williams, Mar 2021                       #
############################################################


import os
from datetime import datetime
import numpy as np


#########################################################################
def add_cosicorr_metadata(fname, cosicorr_dates, meta):
    '''Read/extract attribute data from cosicorr metadata file and add to metadata dictionary
    Inputs:
        Offset or SNR file name (fname)
        dictionary of file name and date12 pairs (cosicorr_dates)
        Metadata dictionary (meta)
    Output:
        Metadata dictionary (meta)
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
    date1_string, date2_string = cosicorr_dates[os.path.basename(fname)].split('-')
    meta['DATE12'] = f'{date1_string}-{date2_string}'
    date1 = datetime.strptime(date1_string,'%Y%m%d')
    date2 = datetime.strptime(date2_string,'%Y%m%d')
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
