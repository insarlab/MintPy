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
from mintpy.utils.arg_utils import create_argument_parser


############################################################################
EXAMPLE = """example:
  temporal_derivative.py  timeseries.h5 
"""

def create_parser(subparsers=None):
    synopsis = 'Calculate the temporal derivative of time-series.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', type=str, help='time-series displacement file.')
    parser.add_argument('-o','--output', dest='outfile', type=str, help='output derivative time-series file.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


############################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    # read data
    obj = timeseries(inps.file)
    obj.open(print_msg=False)
    print('reading timeseries data from file: {}'.format(inps.file))
    ts_data = obj.read(print_msg=False)

    # calculation
    print('calculate the 1st derivative of timeseries data')
    ts_data_1d = np.zeros(ts_data.shape, np.float32)
    ts_data_1d[1:, :, :] = np.diff(ts_data, n=1, axis=0)

    # write to file
    if not inps.outfile:
        inps.outfile = '{}_1stDiff.h5'.format(os.path.splitext(inps.file)[0])
    obj_out = timeseries(inps.outfile)
    obj_out.write2hdf5(ts_data_1d, refFile=inps.file)

    return


############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
