#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Mahmud Haghighi, Mar 2025       #
############################################################


import sys

from mintpy.utils import readfile
from mintpy.utils.arg_utils import create_argument_parser

####################################################################################
DESCRIPTION = """More information:
Documentation for InSAR Explorer: https://luhipi.github.io/insar-explorer
Install it on QGIS via Plugins > Manage and Install Plugins....
"""
EXAMPLE = """example:
  save_explorer.py  geo_timeseries_demErr.h5
  save_explorer.py  geo_timeseries_demErr.h5 -v geo_velocity.h5 -m geo_maskTempCoh.h5
  save_explorer.py  geo_timeseries_demErr.h5 -v geo_velocity.h5 -o timeseries -m geo_maskTempCoh.h5
  save_explorer.py  geo_timeseries_demErr_mask.h5 -v geo_velocity_mask.h5 -o timeseries -m geo_maskTempCoh.h5
"""



def create_parser(subparsers=None):
    synopsis = 'Convert time series to GRD files compatible with InSAR Explorer. '
    epilog =  EXAMPLE + '\n' + DESCRIPTION
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file',
                        help='Time series file to be converted, in geo coordinate.')
    parser.add_argument('-v', '--vel', dest='vel_file',
                        help='velocity file to be converted, in geo coordinate.')
    parser.add_argument('-m', '--mask', dest='mask_file',
                        help='mask file, in geo coordinates. Default: no mask applied.')
    parser.add_argument('-o', '--output', dest='outdir',
                        default='InSAR-Explorer',
                        help='Name of the output directory where files will be created. Default: InSAR-Explorer')
    return parser


def cmd_line_parse(iargs=None):
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check
    atr = readfile.read_attribute(inps.file)

    # check: input file coordinate system
    if 'Y_FIRST' not in atr.keys():
        raise Exception('ERROR: input file is not geocoded.')

    ftype = atr['FILE_TYPE']
    if not ftype in ['timeseries']:
        raise Exception(f"NO required timeseries found in file {inps.file}!")

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
