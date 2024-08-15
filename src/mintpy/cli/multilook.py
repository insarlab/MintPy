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
  multilook.py velocity.h5 -r 15 -a 15
  multilook.py srtm30m.dem -x 10 -y 10 -o srtm300m.dem
  multilook.py filt_fine.int -r 2 -a 2 -o filt_fine_mli.int

  # support GDAL VRT file from ISCE2 as input
  multilook.py lat.rdr.full.vrt lon.rdr.full.vrt -x 9 -y 3

  # --off-file option: use as reference to adjust for the irregular size from isce2 dense offsets
  multilook.py lat.rdr.full.vrt -x 128 -y 64 -o lat.rdr.mli --off-file dense_offsets.bil
  multilook.py ../../geom_reference/lat.rdr.full -x 300 -y 100 -o lat.rdr --off-file offset.bip
"""


def create_parser(subparsers=None):
    synopsis = 'Multilook the input file'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    # basic
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

    # offset
    ampcor = parser.add_argument_group('Ampcor options', 'Ampcor options for dense offsets to account for the extra margin')
    ampcor.add_argument('--search','--search-win', dest='search_win', type=int, nargs=2, metavar=('X','Y'),
                        help='Ampcor (half) search window in (width, height) in pixel, e.g. 20 x 20.')
    ampcor.add_argument('--xcorr','--xcorr-win', dest='xcorr_win', type=int, nargs=2, metavar=('X','Y'),
                        help='Ampcor cross-correlation window in (width, height) in pixel e.g. 32 x 32.')
    ampcor.add_argument('--margin', dest='margin', type=int, default=0,
                        help='Ampcor margin offset (default: %(default)s).')
    ampcor.add_argument('--off-file', dest='off_file', type=str,
                        help='Ampcor offset file as reference for the size.')

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
            search_win=inps.search_win,
            xcorr_win=inps.xcorr_win,
            margin=inps.margin,
            off_file=inps.off_file,
        )
    print('Done.')


###################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
