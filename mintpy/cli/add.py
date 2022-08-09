############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
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
    from ..utils import readfile
    
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    if len(inps.file) < 2:
        parser.print_usage()
        sys.exit('ERROR: At least 2 input files needed!')

    # only two input files are supported for time-series type
    atr = readfile.read_attribute(inps.file[0])
    if atr['FILE_TYPE'] == 'timeseries' and len(inps.file) != 2:
        raise ValueError('Only TWO files are supported for time-series, input has {}'.format(len(inps.file)))

    return inps


################################################################################
def main(iargs=None):
    from ..add import add_file
    inps = cmd_line_parse(iargs)
    print('input files to be added: ({})\n{}'.format(len(inps.file), inps.file))
    add_file(inps.file, inps.outfile, force=inps.force)
    print('Done.')


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
