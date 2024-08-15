#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Heresh Fattahi, Aug 2022      #
############################################################


import sys

from mintpy.utils.arg_utils import create_argument_parser

################################################################################################
REFERENCE = """references:
  Bekaert, D. P. S., Handwerger, A. L., Agram, P., & Kirschbaum, D. B. (2020). InSAR-based detection
    method for mapping and monitoring slow-moving landslides in remote regions with steep and
    mountainous terrain: An application to Nepal. Remote Sensing of Environment, 249, 111983.
    doi:10.1016/j.rse.2020.111983
"""

EXAMPLE = """example:
  spatial_filter.py  velocity.h5
  spatial_filter.py  timeseries.h5 -f lowpass_avg       -p 5
  spatial_filter.py  velocity.h5   -f lowpass_avg       -p 5
  spatial_filter.py  velocity.h5   -f highpass_gaussian -p 3
  spatial_filter.py  velocity.h5   -f sobel
  spatial_filter.py  ifgramStack.h5 unwrapPhase
  spatial_filter.py  ifgramStack.h5 unwrapPhase -f lowpass_avg -p 5
  spatial_filter.py  ifgramStack.h5 unwrapPhase -f double_difference -p 1 10
"""


def create_parser(subparsers=None):
    synopsis = 'Spatial filtering of 2D image.'
    epilog = REFERENCE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', help='File to be filtered')
    parser.add_argument('dset', type=str, nargs='*', default=[],
                        help='optional - dataset(s) to filter (default: %(default)s).')
    parser.add_argument('-f', dest='filter_type', default='lowpass_gaussian',
                        choices=['lowpass_gaussian', 'highpass_gaussian',
                                 'lowpass_avg', 'highpass_avg',
                                 'sobel', 'roberts', 'canny', 'double_difference', 'median'],
                        help='Filter type (default: %(default)s).\n' +
                             'Check Bekaert et al. (2020) for double_difference;\n' +
                             'Check scikit-image as below for the other filters:\n' +
                             '    http://scikit-image.org/docs/dev/api/skimage.filters.html')
    parser.add_argument('-p', '--filter_par', dest='filter_par', nargs='*', type=float,
                        help='Filter parameters for filters. Default:\n' +
                             '    Sigma         for low/high pass gaussian filter, default: 3.0\n' +
                             '    Kernel Size   for low/high pass average  filter, default: 5\n' +
                             '    Kernel Radius for double difference local and regional filters, default: 1 10\n'+
                             '    Kernel Radius for median filters, default: 5\n')
    parser.add_argument('-o', '--outfile',default=None, help='Output file name.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


################################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.spatial_filter import filter_file

    # run
    filter_file(
        fname=inps.file,
        ds_names=inps.dset,
        filter_type=inps.filter_type,
        filter_par=inps.filter_par,
        fname_out=inps.outfile,
    )


################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
