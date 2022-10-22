#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import os
import sys

from mintpy.utils.arg_utils import create_argument_parser

##############################################################################
EXAMPLE = """example:
  load_gbis.py invert_1_2_C.mat
  load_gbis.py invert_1_2_C.mat --nodisplay
"""

def create_parser(subparsers=None):
    synopsis = 'Load GBIS inversion result to HDF5 format.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', help='GBIS inversion mat file.')
    parser.add_argument('-o', '--output', dest='outfile', help='output file name.')
    parser.add_argument('--nodisplay', dest='disp_fig', action='store_false', help='do not display the figure')
    return parser


def cmd_line_parse(iargs=None):
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # import
    import matplotlib.pyplot as plt

    # check: input file path
    inps.file = os.path.abspath(inps.file)

    # check: --nodisplay (matplotlib backend setting)
    if not inps.disp_fig:
        plt.switch_backend('Agg')

    return inps


##############################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.load_gbis import gbis_mat2hdf5

    # run
    gbis_mat2hdf5(inps.file, display=inps.disp_fig)


##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
