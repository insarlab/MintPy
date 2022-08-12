############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import sys
from mintpy.utils.arg_utils import create_argument_parser


############################################################
REFERENCE="""reference:
  Wikipedia: https://en.wikipedia.org/wiki/Gaussian_blur
"""

EXAMPLE = """example:
 temporal_filter.py timeseries_ERA5_demErr.h5
 temporal_filter.py timeseries_ERA5_demErr.h5 -t 0.1
"""

def create_parser(subparsers=None):
    synopsis = 'Smoothing timeseries in time domain with a moving filter'
    epilog = REFERENCE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('timeseries_file',
                        help='timeseries file to be smoothed.')
    parser.add_argument('-t', '--time-win', dest='time_win', type=float, default=0.1,
                        help='time window in years, default: 0.1 (Sigma of the assmued Gaussian distribution.)')
    parser.add_argument('-o', '--outfile', help='Output file name.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


############################################################
def main(iargs=None):
    from ..temporal_filter import temporal_filter
    inps = cmd_line_parse(iargs)
    temporal_filter(inps.timeseries_file, inps.outfile, inps.time_win)


############################################################
if __name__ == '__main__':
    main(sys.argv[1:])
