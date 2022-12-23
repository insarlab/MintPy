#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import sys

from mintpy.utils.arg_utils import create_argument_parser

#######################################################################################
EXAMPLE = """example:
  image_math.py  velocity.h5            '+'  0.5
  image_math.py  geo_080212_101120.cor  '-'  0.2
  image_math.py  timeseries.h5          '*'  1.5
  image_math.py  velocity.h5            '/'  2.0
  image_math.py  velocity.h5            '^'  2.0

  # update isce2 offset files
  image_math.py  azimuth_*.off '+' 0.0077 --overwrite
"""


def create_parser(subparsers=None):
    synopsis = 'Basic Mathmatic Operation of file'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', nargs='+', help='input file(s).')
    parser.add_argument('-o', '--output', dest='outfile',
                        help='output file name.')
    parser.add_argument('--overwrite', '--overwrite-input', dest='overwrite', action='store_true',
                        help='overwrite input file (disabled if -o/--output is turned on)')
    parser.add_argument('operator', choices=[
                        '+', '-', '*', '/', '^'], help='mathmatical operator')
    parser.add_argument('operand', metavar='VALUE', type=float,
                        help='value to be operated with input file')
    return parser


def cmd_line_parse(iargs=None):
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check: --overwrite vs --output
    if inps.outfile and inps.overwrite:
        print('-o/--output is turned ON, disable --overwrite option and continue.')
        inps.overwrite = False

    return inps


#######################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.image_math import file_operation

    # run
    num_file = len(inps.file)
    for i, fname in enumerate(inps.file):
        if num_file > 1:
            print('-'*20 + f' {i+1}/{num_file} {fname}')
        out_file = fname if inps.overwrite else inps.outfile
        file_operation(fname, inps.operator, inps.operand, out_file)


#######################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
