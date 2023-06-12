#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Heresh Fattahi, Aug 2022      #
############################################################


import sys

from mintpy.utils.arg_utils import create_argument_parser

#############################################################################################
EXAMPLE = """example:
  image_stitch.py  vel_AlosAT422.h5  vel_AlosAT423.h5  vel_AlosAT424.h5  vel_AlosAT425.h5 -o  vel_AlosA.h5
  image_stitch.py geom_AlosAT422.h5 geom_AlosAT423.h5 geom_AlosAT424.h5 geom_AlosAT425.h5 -o geom_AlosA.h5 --no-offset
"""

NOTE = """
  The function automatically:
  1) finds the common area between adjacent input files
  2) calculates the average offset between them
  3) apply this average offset to the later file
"""

def create_parser(subparsers=None):
    synopsis = 'Stitch/mosaic multiple geocoded datasets into one.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis+NOTE, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file1', help='file to stitch')
    parser.add_argument('file2s', nargs='+', metavar='file2', help='file(s) to stitch')
    parser.add_argument('-o', '--output', dest='outfile', required=True, help='output file name')

    # stitch option
    parser.add_argument('--no-offset','--no-off', dest='apply_offset', action='store_false',
                        help='Do not apply offset if 1) data sets are merely to be stitched '
                             'AND 2) no adjustment of values needs to be made\n'
                             '(i.e., for two coherence maps), use this flag')

    # plot options
    parser.add_argument('--scale', dest='disp_scale', type=float, default=1.,
                        help='scale the data when plotting.')
    parser.add_argument('-v','--vlim', dest='disp_vlim', type=float, default=None, nargs=2,
                        help='vmin and vmax when plotting.')
    parser.add_argument('-c','--cmap', dest='disp_cmap', type=str, default=None,
                        help='colormap when plotting.')
    parser.add_argument('--nodisplay', dest='disp_fig', action='store_false',
                        help='do not display the result plotting.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


#############################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.image_stitch import stitch_files

    # run
    stitch_files(
        fnames=[inps.file1] + inps.file2s,
        out_file=inps.outfile,
        apply_offset=inps.apply_offset,
        disp_fig=inps.disp_fig,
        disp_scale=inps.disp_scale,
        disp_vlim=inps.disp_vlim,
        disp_cmap=inps.disp_cmap,
    )


#############################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
