#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Mahmud Haghighi, Mar 2025                        #
############################################################


import sys

from mintpy.utils import readfile
from mintpy.utils.arg_utils import create_argument_parser

####################################################################################
EXAMPLE = """example:
  save_explorer.py geo_timeseries_demErr.h5
  save_explorer.py geo_timeseries_demErr.h5 -v geo_velocity.h5 -m geo_maskTempCoh.h5
  save_explorer.py timeseries_demErr.h5 -v velocity.h5 -m maskTempCoh.h5 -o timeseries
  # check more info on InSAR Explorer documentation: https://luhipi.github.io/insar-explorer
"""


def create_parser(subparsers=None):
    synopsis = 'Convert time series to GRD files for QGIS InSAR Explorer plugin.'
    epilog =  EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('ts_file',
                        help='Time series file to be converted, in geo coordinate.')
    parser.add_argument('-v', '--vel', dest='vel_file',
                        help='velocity file to be converted, in geo coordinate.')
    parser.add_argument('-m', '--mask', dest='msk_file',
                        help='mask file, in geo coordinates.')
    parser.add_argument('-o', '--output', dest='outdir', default='InSAR-Explorer',
                        help='Output directory for GRD files (default: %(default)s).')
    return parser


def cmd_line_parse(iargs=None):
    '''Command line parser.'''
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check: input time series file is geocoded time series
    atr = readfile.read_attribute(inps.ts_file)
    ftype = atr['FILE_TYPE']
    if 'Y_FIRST' not in atr.keys():
        raise Exception('Input file is NOT geocoded! Only geocoded files are supported.')
    if ftype not in ['timeseries']:
        raise Exception(f'Input file ({ftype}) is NOT time series!')

    return inps


####################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.save_explorer import save_explorer

    # run
    save_explorer(inps)


####################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
