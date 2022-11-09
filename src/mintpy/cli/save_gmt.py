#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Heresh Fattahi, Aug 2022      #
############################################################


import sys

from mintpy.utils.arg_utils import create_argument_parser

####################################################################################
EXAMPLE = """example:
  save_gmt.py  geo_velocity.h5
  save_gmt.py  geo_timeseries.h5  20071031
  save_gmt.py  geo_filt_fine.unw
  save_gmt.py  gsi10m.dem
"""


def create_parser(subparsers=None):
    synopsis = 'Export geocoded file to GMT grd file'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', help='file to be converted, in geo coordinate.')
    parser.add_argument('dset', nargs='?',
                        help='date of timeseries, or date12 of interferograms to be converted')
    parser.add_argument('-o', '--output', dest='outfile',
                        help='output file base name. Extension is fixed with .kmz')
    return parser


def cmd_line_parse(iargs=None):
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # import
    from mintpy.utils import readfile

    # check
    atr = readfile.read_attribute(inps.file)

    # check: input file coordinate system
    if 'Y_FIRST' not in atr.keys():
        raise Exception('ERROR: input file is not geocoded.')

    # check: dset for certain file types (timeseries and ifgramStack)
    ftype = atr['FILE_TYPE']
    if not inps.dset and ftype in ['timeseries', 'ifgramStack']:
        raise Exception(f"NO required dataset is given for {ftype} file!")

    return inps


####################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.save_gmt import save_gmt

    # run
    save_gmt(inps)


####################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
