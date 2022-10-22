#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import sys

from mintpy.utils.arg_utils import create_argument_parser

##################################################################################################
NOTE = """
  For each binary file (unwrapped/wrapped interferogram, spatial coherence file), there are 2 .rsc files:
  1) basic metadata file and 2) baseline parameter file. This script find those two rsc files based on
  input binary file name, and merge those two metadata files into one.

  For example, if input binary file is filt_100901-110117-sim_HDR_4rlks_c10.unw, this script will find
  1) filt_100901-110117-sim_HDR_4rlks_c10.unw.rsc and 2) 100901-110117_baseline.rsc and merge 1) and 2) into
  one file: filt_100901-110117-sim_HDR_4rlks_c10.unw.rsc
"""

EXAMPLE = """example:
  prep_roipac.py  filt_100901-110117-sim_HDR_4rlks_c10.unw
  prep_roipac.py  ./interferograms/*/filt_*.unw
  prep_roipac.py  ./interferograms/*/filt_*rlks.cor
  prep_roipac.py  ./interferograms/*/filt_*rlks.int
  prep_roipac.py  ./interferograms/*/filt_*_snap_connect.byt
"""


def create_parser(subparsers=None):
    synopsis = 'Prepare attributes file for ROI_PAC products.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis+NOTE, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', nargs='+', help='Gamma file(s)')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


##################################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.prep_roipac import prep_roipac

    # run
    prep_roipac(inps)


###################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
