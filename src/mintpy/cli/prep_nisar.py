#!/usr/bin/env python
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Sara Mirzaee, Jul 2023                           #
############################################################

import sys

from mintpy.utils.arg_utils import create_argument_parser

############################################################
EXAMPLE = """example:
  prep_nisar.py -i 'interferograms/stitched/*.h5' -d dem.tiff
"""

def create_parser(subparsers=None):
    """Command line parser."""
    synopsis = 'Prepare NISAR GUNW products for MintPy.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument(
        "-i",
        "--input-file-glob",
        dest='input_glob',
        type=str,
        default="./interferograms/NISAR*.h5",
        help="path pattern of NISAR interferograms (default: %(default)s).",
    )

    parser.add_argument(
        "-d",
        "--dem-file",
        dest="dem_file",
        type=str,
        default="./dem.tif",
        help="path to the DEM (default: %(default)s).",
    )

    parser.add_argument(
        "-o",
        "--out-dir",
        dest="out_dir",
        type=str,
        default="./mintpy",
        help="output directory (default: %(default)s).",
    )

    parser.add_argument(
        '--force',
        dest='update_mode',
        action='store_false',
        help='Force to overwrite all .rsc metadata files.'
    )
    parser.add_argument(
        '--sub-lat',
        '--sublat',
        '--subset-lat',
        dest='subset_lat',
        type=float,
        nargs=2,
        metavar=('LATMIN', 'LATMAX'),
        help='subset in latitude'
    )

    parser.add_argument(
        '--sub-lon',
        '--sublon',
        '--subset-lon',
        dest='subset_lon',
        type=float,
        nargs=2,
        metavar=('LONMIN', 'LONMAX'),
        help='subset in longitude'
    )

    return parser


def cmd_line_parse(iargs=None):
    """Create the command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    return inps

####################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.prep_nisar import load_nisar

    # run
    load_nisar(inps)


####################################################################################
if __name__=="__main__":
    main(sys.argv[1:])
