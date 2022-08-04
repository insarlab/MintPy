#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, Zhang Yunjun, 2013               #
############################################################


import os
import sys
import numpy as np
from mintpy.objects import timeseries
from mintpy.utils import ptime
from mintpy.utils.arg_utils import create_argument_parser


############################################################
REFERENCE="""reference:
  Wikipedia: https://en.wikipedia.org/wiki/Gaussian_blur
"""

EXAMPLE = """example:
 temporal_filter.py timeseries_ERA5_demErr.h5
 temporal_filter.py timeseries_ERA5_demErr.h5 -t 0.1
"""

def create_parser(subparsers=None):
    synopsis = 'Smoothing timeseries in time domain with a moving filter'
    epilog = REFERENCE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('timeseries_file',
                        help='timeseries file to be smoothed.')
    parser.add_argument('-t', '--time-win', dest='time_win', type=float, default=0.1,
                        help='time window in years, default: 0.1 (Sigma of the assmued Gaussian distribution.)')
    parser.add_argument('-o', '--outfile', help='Output file name.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


############################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    # read timeseries info / data
    obj = timeseries(inps.timeseries_file)
    obj.open()

    tbase = np.array(obj.yearList, np.float32).reshape(-1, 1)
    tbase -= tbase[obj.refIndex]

    ts_data = obj.read().reshape(obj.numDate, -1)

    # Smooth acquisitions / moving window in time one by one
    print('-'*50)
    print('filtering in time Gaussian window with size of {:.1f} years'.format(inps.time_win))
    ts_data_filt = np.zeros(ts_data.shape, np.float32)
    prog_bar = ptime.progressBar(maxValue=obj.numDate)
    for i in range(obj.numDate):
        # Weight from Gaussian (normal) distribution in time
        tbase_diff = tbase[i] - tbase
        weight = np.exp(-0.5 * (tbase_diff**2) / (inps.time_win**2))
        weight /= np.sum(weight)
        # Smooth the current acquisition
        ts_data_filt[i, :] = np.sum(ts_data * weight, axis=0)
        prog_bar.update(i+1, suffix=obj.dateList[i])
    prog_bar.close()
    del ts_data
    ts_data_filt -= ts_data_filt[obj.refIndex, :]
    ts_data_filt = np.reshape(ts_data_filt, (obj.numDate, obj.length, obj.width))

    # write filtered timeseries file
    if not inps.outfile:
        inps.outfile = '{}_tempGaussian.h5'.format(os.path.splitext(inps.timeseries_file)[0])
    obj_out = timeseries(inps.outfile)
    obj_out.write2hdf5(ts_data_filt, refFile=inps.timeseries_file)

    return


############################################################
if __name__ == '__main__':
    main(sys.argv[1:])
