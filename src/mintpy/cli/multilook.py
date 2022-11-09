#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import sys

from mintpy.utils.arg_utils import create_argument_parser

##################################################################################################
EXAMPLE = """example:
  multilook.py  velocity.h5  -r 15 -a 15
  multilook.py  srtm30m.dem  -r 10 -a 10  -o srtm30m_300m.dem

  # Ignore / skip marginal pixels
  multilook.py ../../geom_reference/hgt.rdr.full -r 300 -a 100 --margin 58 58 58 58 -o hgt.rdr
"""


def create_parser(subparsers=None):
    synopsis = 'Multilook the input file'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', nargs='+', help='File(s) to multilook')
    parser.add_argument('-r','--range','-x', dest='lks_x', type=int, default=1,
                        help='number of multilooking in range  /x direction (default: %(default)s).')
    parser.add_argument('-a','--azimuth','-y', dest='lks_y', type=int, default=1,
                        help='number of multilooking in azimuth/y direction (default: %(default)s).')
    parser.add_argument('-o', '--outfile',
                        help='Output file name. Disabled when more than 1 input files')
    parser.add_argument('-m','--method', dest='method', type=str, default='mean', choices=['mean', 'median', 'nearest'],
                        help='downsampling method (default: %(default)s) \n'
                             'e.g. nearest for geometry, average for observations')
    parser.add_argument('--margin', dest='margin', type=int, nargs=4, metavar=('TOP','BOTTOM','LEFT','RIGHT'),
                        default=[0,0,0,0], help='number of pixels on the margin to skip, (default: %(default)s).')
    return parser


def cmd_line_parse(iargs=None):
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # import
    from mintpy.utils import utils1 as ut

    # check: input file list
    inps.file = ut.get_file_list(inps.file)

    # check: -x/y options (num of multilooks)
    if inps.lks_x == 1 and inps.lks_y == 1:
        raise SystemExit('ERROR: no multilooking specified: lks_x/y=1!')

    # check: -o / --outfile (output file name)
    if len(inps.file) > 1 and inps.outfile:
        inps.outfile = None
        print('more than one file is input, disable custom output filename.')

    return inps


##################################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.multilook import multilook_file

    # run
    for infile in inps.file:
        multilook_file(
            infile,
            lks_y=inps.lks_y,
            lks_x=inps.lks_x,
            outfile=inps.outfile,
            method=inps.method,
            margin=inps.margin)
    print('Done.')


###################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
