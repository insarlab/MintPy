#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Yunmeng Cao, Aug 2022         #
############################################################


import sys

from mintpy.utils.arg_utils import create_argument_parser

################################################################################
EXAMPLE = '''examples:
    lookup_geo2radar.py geometryGeo.h5
    lookup_geo2radar.py geometryGeo.h5 -w geometryRadar.h5
    lookup_geo2radar.py geometryGeo.h5 -w geometryRadar.h5 --parallel 4
'''

def create_parser(subparsers=None):
    synopsis = 'Convert lookup table from geo-coord (GAMMA, ROI_PAC) into radar-coord (ISCE)'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('geom_geo_file', help='The geometryGeo.h5 file which includes geo-coordinates based lookup-table')
    parser.add_argument('-w','--write', dest='write', metavar='FILE', default='geometryRadar.h5',
                        help='update geometryRadar.h5 file by adding the radar-coordinates based lookup-table.')
    parser.add_argument('--parallel', dest='parallelNumb', type=int, metavar='NUM', default=1,
                        help='Enable parallel processing and specify the the used processor number.[default: 1]')

    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.lookup_geo2radar import run_lookup_geo2radar

    # run
    run_lookup_geo2radar(inps)


##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
