#!/usr/bin/env python
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Sara Mirzaee                                     #
############################################################

import sys
import argparse

from mintpy.utils import arg_utils

############################################################
EXAMPLE = """example:
  python3 ./prep_nisar.py -i 'interferograms/stitched/*.h5' -d 'dem.tiff'  
"""

def _create_parser():
    parser = argparse.ArgumentParser(
        description="Prepare NISAR products for MintPy",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=EXAMPLE,
    )

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

    parser = arg_utils.add_subset_argument(parser, geo=True)

    return parser


def cmd_line_parse(iargs=None):
    """Create the command line parser."""
    parser = _create_parser()
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
