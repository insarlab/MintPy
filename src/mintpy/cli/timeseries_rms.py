#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import sys

from mintpy.defaults.template import get_template_content
from mintpy.utils.arg_utils import create_argument_parser

######################################################################################################
TEMPLATE = get_template_content('residual_RMS')

REFERENCE="""reference:
  Yunjun, Z., Fattahi, H. and Amelung, F. (2019), Small baseline InSAR time series analysis:
    Unwrapping error correction and noise reduction, Computers & Geosciences, 133, 104331,
    doi:10.1016/j.cageo.2019.104331.
  Rousseeuw, P. J., and M. Hubert (2011), Robust statistics for outlier detection,
    Wiley Interdisciplinary Reviews: Data Mining and Knowledge Discovery, 1(1),
    73-79, doi:doi:10.1002/widm.2.
"""

EXAMPLE = """example:
  timeseries_rms.py  timeseriesResidual.h5
  timeseries_rms.py  timeseriesResidual.h5  --template smallbaselineApp.cfg
  timeseries_rms.py  timeseriesResidual.h5  -m maskTempCoh.h5  --cutoff 3
"""


def create_parser(subparsers=None):
    synopsis = 'Calculate Root Mean Square (RMS) of deramped residual phase time-series.'
    epilog = TEMPLATE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('timeseries_file', help='Timeseries file')
    parser.add_argument('-t', '--template', dest='template_file',
                        help='template file with options')
    parser.add_argument('-m', '--mask', dest='maskFile', default='maskTempCoh.h5',
                        help='mask file for estimation')
    parser.add_argument('-r','--ramp','--deramp', dest='deramp', default='quadratic',
                        help='ramp type to be remove for RMS calculation.\n' +
                             'Default - quadratic; no - do not remove ramp')
    parser.add_argument('--cutoff', dest='cutoff', default='3', type=float,
                        help='M-score used for outlier detection based on standardised residuals\n'+
                             'Recommend range: [3, 4], default is 3.')
    parser.add_argument('--figsize', dest='fig_size', metavar=('WID', 'LEN'),
                        type=float, nargs=2, default=[5., 3.],
                        help='figure size in inches - width and length')
    parser.add_argument('--tick-year-num', dest='tick_year_num',
                        type=int, default=1, help='Year number per major tick')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # save argv (to check the manually specified arguments)
    # use iargs        for python call
    # use sys.argv[1:] for command line call
    inps.argv = iargs if iargs else sys.argv[1:]

    return inps


######################################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.timeseries_rms import run_timeseries_rms

    # run
    run_timeseries_rms(inps)


######################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
