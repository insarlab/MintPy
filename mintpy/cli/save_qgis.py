############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import os
import sys
from mintpy.utils.arg_utils import create_argument_parser


#########################################################################################
EXAMPLE = """example:
  save_qgis.py timeseries_ERA5_ramp_demErr.h5 -g inputs/geometrygeo.h5
  save_qgis.py timeseries_ERA5_ramp_demErr.h5 -g inputs/geometryRadar.h5
  save_qgis.py geo/geo_timeseries_ERA5_ramp_demErr.h5 -g geo/geo_geometryRadar.h5
  save_qgis.py timeseries_ERA5_ramp_demErr.h5 -g inputs/geometryRadar.h5 -b 200 150 400 350
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

    return parser


def cmd_line_parse(iargs=None):
    '''Command line parser.'''
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # --outshp option
    if not inps.shp_file:
        inps.shp_file = os.path.splitext(inps.ts_file)[0] + '.shp'

    return inps


#########################################################################################
def main(iargs=None):
    from ..save_qgis import read_bounding_box, gather_files, write_shape_file

    # Parse command line
    inps = cmd_line_parse(iargs)

    # Read bounding box
    box = read_bounding_box(pix_box=inps.pix_bbox,
                            geo_box=inps.geo_bbox,
                            geom_file=inps.geom_file)

    # Gather data files
    fDict = gather_files(inps.ts_file, inps.geom_file)

    # Write shape file
    write_shape_file(fDict, inps.shp_file, box=box)


#########################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
