#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Forrest Williams, Aug 2022    #
############################################################


import sys

from mintpy.utils.arg_utils import create_argument_parser

EXAMPLE = """example:
  save_gdal.py geo/geo_velocity.h5
  save_gdal.py geo/geo_timeseries_ERA5_demErr.h5  -d 20200505_20200517 --of ENVI
  save_gdal.py geo/geo_ifgramStack.h5 -d unwrapPhase-20101120_20110220 --of ISCE
  save_gdal.py geo/geo_ifgramStack.h5 -d   coherence-20101120_20110220 --of ISCE
  save_gdal.py geo_20230225.slc
"""


def create_parser(subparsers=None):
    synopsis = 'Generate GDAL raster from MintPy h5 file.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', help='file to be converted, in geo coordinate.')
    parser.add_argument('-d', '--dset', '--dataset', dest='dset',
                        help='date of timeseries, or date12 of interferograms to be converted')
    parser.add_argument('-o', '--output', dest='outfile',
                        help='output file base name. Extension is fixed by GDAL driver')
    parser.add_argument('--of', '--out-format', '--output-format', dest='out_format', default='GTiff',
                        help='file format as defined by GDAL driver name, e.g. GTiff, ENVI, default: %(default)s\n'
                             'GDAL driver names can be found at https://gdal.org/drivers/raster/index.html')
    return parser


def cmd_line_parse(iargs=None):
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # import
    from mintpy.utils import readfile

    # check: input file coordinate system
    atr = readfile.read_attribute(inps.file)
    if 'X_FIRST' not in atr.keys():
        raise ValueError(f'ERROR: Input file ({inps.file}) is not geocoded!')

    return inps


##############################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.save_gdal import save_gdal

    # run
    save_gdal(inps)


##############################################################################
if __name__ == "__main__":
    main(sys.argv[1:])
