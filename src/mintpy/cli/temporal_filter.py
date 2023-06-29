#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, Antonio Valentino, 2013          #
############################################################


import sys

from mintpy.utils.arg_utils import create_argument_parser

############################################################
REFERENCE="""reference:
  Wikipedia: https://en.wikipedia.org/wiki/Gaussian_blur
  Scipy: https://docs.scipy.org/doc/scipy/reference/ndimage.html
"""

EXAMPLE = """example:
  temporal_filter.py timeseries_ERA5_demErr.h5 -f gaussian -t 1
  temporal_filter.py timeseries_ERA5_demErr.h5 -f median   -t 5
"""

def create_parser(subparsers=None):
    synopsis = 'Smoothing timeseries in time domain with a moving filter'
    epilog = REFERENCE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('timeseries_file', help='timeseries file to be smoothed.')
    parser.add_argument('-f', '--filter', dest='filter_type', type=str, default='median',
                        choices={'median', 'gaussian'},
                        help='temporal filter type (default: %(default)s).')
    parser.add_argument('-t', '--time-win', dest='time_win', type=float, default=5,
                        help='time window size (default: %(default)s).\n'
                             'number of months       for gaussian filter\n'
                             'number of acquisitions for median   filter')
    parser.add_argument('-o', '--outfile', help='Output file name.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


############################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.objects import timeseries

    # run
    ts_obj = timeseries(inps.timeseries_file)
    ts_obj.temporal_filter(
        time_win=inps.time_win,
        filter_type=inps.filter_type,
        out_file=inps.outfile)

    return


############################################################
if __name__ == '__main__':
    main(sys.argv[1:])
