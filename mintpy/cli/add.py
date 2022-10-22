#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import sys

from mintpy.utils.arg_utils import create_argument_parser

################################################################################
EXAMPLE = """example:
  add.py mask_1.h5 mask_2.h5 mask_3.h5        -o mask_all.h5
  add.py 081008_100220.unw  100220_110417.unw -o 081008_110417.unw
  add.py timeseries_ERA5.h5 inputs/ERA5.h5    -o timeseries.h5
  add.py timeseriesRg.h5    inputs/TECsub.h5  -o timeseriesRg_TECsub.h5 --force
"""


def create_parser(subparsers=None):
    """ Command line parser """
    synopsis = 'Generate the sum of multiple input files.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', nargs='+', help='files (2 or more) to be added')
    parser.add_argument('-o', '--output', dest='outfile', help='output file name')
    parser.add_argument('--force', action='store_true',
                        help='Enforce the adding for the shared dates only for time-series files')
    return parser


def cmd_line_parse(iargs=None):
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # import
    from mintpy.utils import readfile

    # check
    num_file = len(inps.file)
    atr = readfile.read_attribute(inps.file[0])

    # check: number of files a) >= 2 and b) == 2 for time-series
    if num_file < 2:
        parser.print_usage()
        sys.exit('ERROR: At least 2 input files needed!')

    elif num_file > 2 and atr['FILE_TYPE'] == 'timeseries':
        raise ValueError(f'Only TWO files are supported for time-series, input has {num_file}.')

    return inps


################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.add import add_file

    # run
    add_file(inps.file, inps.outfile, force=inps.force)


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
