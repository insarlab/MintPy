############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import sys
from mintpy.utils.arg_utils import create_argument_parser


############################################################################
EXAMPLE = """example:
  temporal_derivative.py  timeseries.h5 
"""

def create_parser(subparsers=None):
    synopsis = 'Calculate the temporal derivative of time-series.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', type=str, help='time-series displacement file.')
    parser.add_argument('-o','--output', dest='outfile', type=str, help='output derivative time-series file.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


############################################################################
def main(iargs=None):
    from ..temporal_derivative import temporal_derivative
    inps = cmd_line_parse(iargs)
    return temporal_derivative(inps.file, inps.outfile)

############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
