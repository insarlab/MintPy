############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import os
import sys
import glob
from mintpy.utils.arg_utils import create_argument_parser


##################################################################################################
EXAMPLE = """example:
  cd $KIRISHIMA/KirishimaAlosAT424/DEM
  dem_gsi.py -b 31.1 32.8 130.1 131.9
  dem_gsi.py -b 31.1 32.8 130.1 131.9 --grid-dir ~/data/DEM/GSI_DEHM10m
"""

NOTE = """DEHM: Digital Ellipsoidal Height Model
yyxx.dehm with yy and xx indicating the coordinates of the upper left corner of the firt pixel.
where longitude = xx + 100
      latitude  = (yy + 1) / 1.5
"""

def create_parser(subparsers=None):
    synopsis = 'Prepare DEM from GSI (Japan) DEHM grib files.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('-b','--bbox', dest='SNWE', type=float, nargs=4, metavar=('S','N','W','E'), required=True,
                        help='Bounding box in latitude [-90, 90] and longitude [-180, 180].')
    parser.add_argument('-o','--output', dest='outfile', default='gsi10m.dem.wgs84',
                        help='output file name (default: %(default)s).')
    parser.add_argument('-g','--grid-dir', dest='grid_dir', default='$DEMDB/GSI_DEHM10m',
                        help='Directory of DEHM grib files (default: %(default)s).')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    inps.grid_dir = os.path.expanduser(inps.grid_dir)
    inps.grid_dir = os.path.expandvars(inps.grid_dir)
    inps.grid_dir = os.path.abspath(inps.grid_dir)
    if len(glob.glob(os.path.join(inps.grid_dir, '*.dehm'))) == 0:
        raise SystemExit('ERROR: no *.dehm file found in directory: {}'.format(inps.grid_dir))
    return inps


##################################################################################################
def main(iargs=None):
    from ..dem_gsi import write_dem_file, write_rsc_file, write_isce_metadata

    inps = cmd_line_parse(iargs)

    meta = write_dem_file(inps.SNWE,
                          dem_file=inps.outfile,
                          grid_dir=inps.grid_dir)

    # rsc file for roipac
    write_rsc_file(meta, inps.outfile)

    # vrt/xml file for isce
    write_isce_metadata(meta, inps.outfile)


###################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
