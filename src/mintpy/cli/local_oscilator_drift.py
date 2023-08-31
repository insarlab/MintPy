#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Heresh Fattahi, Aug 2022      #
############################################################


import sys

from mintpy.defaults.template import get_template_content
from mintpy.utils.arg_utils import create_argument_parser

#########################################################################################
TEMPLATE = get_template_content('correct_LOD')

REFERENCE = """reference:
  Marinkovic, P., and Y. Larsen (2013), Consequences of long-term ASAR local oscillator
    frequency decay - An empirical study of 10 years of data, 2013 Living Planet Symposium,
    Edinburgh, U.K.
"""

EXAMPLE = """example:
  local_oscilator_drift.py  timeseries.h5                 inputs/geometryRadar.h5
  local_oscilator_drift.py  filt_101020_110220_4rlks.unw  inputs/geometryRadar.h5
"""

def create_parser(subparsers=None):
    synopsis = 'Local Oscillator Drift (LOD) correction of Envisat'
    epilog = REFERENCE + '\n' + TEMPLATE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument(dest='file', help='timeseries / interferograms file, i.e. timeseries.h5')
    parser.add_argument(dest='range_dist_file',
                        help='Slant range distance file, i.e. inputs/geometryRadar.h5, inputs/geometryGeo.h5\n' +
                        'or use range_distance.py to generate it.')
    parser.add_argument('-o', '--output', dest='out_file',
                        help='Output file name for corrected file.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


#########################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.local_oscilator_drift import correct_local_oscilator_drift

    # run
    inps.outfile = correct_local_oscilator_drift(
        inps.file,
        rg_dist_file=inps.range_dist_file,
        out_file=inps.out_file,
    )


#########################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
