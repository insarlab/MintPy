#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Piyush Agram, Aug 2022        #
############################################################


import os
import sys

from mintpy.utils.arg_utils import create_argument_parser

#########################################################################################
EXAMPLE = """example:
  save_qgis.py timeseries_ERA5_ramp_demErr.h5 -g inputs/geometrygeo.h5
  save_qgis.py timeseries_ERA5_ramp_demErr.h5 -g inputs/geometryRadar.h5
  save_qgis.py timeseries_ERA5_ramp_demErr.h5 -g inputs/geometryRadar.h5 -b 200 150 400 350
  save_qgis.py geo/geo_timeseries_ERA5_ramp_demErr.h5 -g geo/geo_geometryRadar.h5
"""

def create_parser(subparsers=None):
    synopsis = 'Convert to QGIS compatible ps time-series'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('ts_file', type=str, help='time-series HDF5 file')
    parser.add_argument('-g', '--geom', dest='geom_file', type=str, required=True,
                        help='geometry HDF5 file')
    parser.add_argument('-o', '--outshp', dest='shp_file', type=str, help='Output shape file.')

    # bounding box
    parser.add_argument('-b', '--bbox', dest='pix_bbox', type=int, nargs=4, default=None,
                        metavar=('Y0','Y1','X0','X1'), help='bounding box : minLine maxLine minPixel maxPixel')
    parser.add_argument('-B', '--geo-bbox', dest='geo_bbox', type=float, nargs=4, default=None,
                        metavar=('S','N','W','E'), help='bounding box in lat lon: South North West East')

    # other options
    parser.add_argument('--zf', '--zero-first', dest='zero_first', action='store_true',
                        help='Set displacement at the first acquisition to zero.')

    return parser


def cmd_line_parse(iargs=None):
    '''Command line parser.'''
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # default: --outshp option
    if not inps.shp_file:
        fbase = os.path.splitext(inps.ts_file)[0]
        inps.shp_file = fbase + '.shp'

    return inps


#########################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.save_qgis import save_qgis

    # run
    save_qgis(inps)


#########################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
