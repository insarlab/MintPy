#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Mar 2024                           #
############################################################


import os
import sys

from mintpy.defaults.template import get_template_content
from mintpy.utils.arg_utils import create_argument_parser

# key configuration parameter name
key_prefix = 'mintpy.ionosphericDelay.'

###############################################################
TEMPLATE = get_template_content('correct_ionosphere')

REFERENCE = """reference:
  Gomba, G., Parizzi, A., Zan, F. D., Eineder, M., & Bamler, R. (2016). Toward Operational Compensation
    of Ionospheric Effects in SAR Interferograms: The Split-Spectrum Method. IEEE Transactions on
    Geoscience and Remote Sensing, 54(3), 1446-1461. doi:10.1109/TGRS.2015.2481079

  # ISCE-2 topsApp/topsStack
  Liang, C., Agram, P., Simons, M., & Fielding, E. J. (2019). Ionospheric Correction of InSAR Time Series
    Analysis of C-band Sentinel-1 TOPS Data. IEEE Trans. Geosci. Remote Sens., 59(9), 6755-6773,
    doi:10.1109/TGRS.2019.2908494.

  # ISCE-2 stripmapApp/stripmapStack
  Fattahi, H., Simons, M., & Agram, P. (2017). InSAR Time-Series Estimation of the Ionospheric Phase
    Delay: An Extension of the Split Range-Spectrum Technique. IEEE Trans. Geosci. Remote Sens., 55(10),
    5984-5996, doi:10.1109/TGRS.2017.2718566

  # ISCE-2 alos2App/alosStack
  Liang, C., Liu, Z., Fielding, E. J., & Bürgmann, R. (2018). InSAR Time Series Analysis of L-Band
    Wide-Swath SAR Data Acquired by ALOS-2. IEEE Trans. Geosci. Remote Sens., 56(8), 4492-4506,
    doi:10.1109/TGRS.2018.2821150
"""

EXAMPLE = """example:
  # CAUTION: Failed split spectrum estimations are common, check & locate pairs with failed estimations,
  # exclude them using mintpy.ionosphericDelay.excludeDate(12) options.
  # For alosStack, check plots in the "fig_ion" folder.

  # estimate ionospheric delay time-series and correct time-series file
  iono_split_spectrum.py -t JinshaAlos2A148.txt -f timeseries_SET.h5

  # estimate ionospheric delay time-series
  iono_split_spectrum.py -t JinshaAlos2A148.txt
"""


def create_parser(subparsers=None):
    synopsis = 'Ionospheric correction using split spectrum (from ISCE-2 stack processing)'
    epilog = REFERENCE + '\n' + TEMPLATE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    # inputs
    parser.add_argument('-t', '--template', dest='template_file', required=True,
                        help='template file with ionospheric delay options.')
    parser.add_argument('--iono-stack-file', dest='iono_stack_file', default='./inputs/ionStack.h5',
                        help='input ionospheric delay stack file (default: %(default)s).')
    parser.add_argument('-f', '--file', dest='dis_file',
                        help='time-series HDF5 file to be corrected, e.g. timeseries.h5')

    # modify network
    parser.add_argument('--ex-date12','--exclude-date12', dest='excludeDate12', nargs='*',
                        help='pair(s) to remove/drop in YYYYMMDD_YYYYMMDD format.')
    parser.add_argument('--ex-date','--exclude-date', dest='excludeDate', nargs='*',
                        help='date(s) to remove/drop, all pairs included date(s) will be removed')

    # outputs
    parser.add_argument('--iono-file', dest='iono_file', default='ion.h5',
                        help='output ionospheric delay time series file (default: %(default)s).')
    parser.add_argument('-o', '--output', dest='cor_dis_file',
                        help='output corrected time-series file, e.g. timeseries_ion.h5')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    # check: existence of the input file
    if inps.dis_file and not os.path.isfile(inps.dis_file):
        raise FileNotFoundError(f'input file not exist: {inps.dis_file}')

    # default: -o / --output option
    if inps.dis_file and not inps.cor_dis_file:
        fbase, fext = os.path.splitext(inps.dis_file)
        inps.cor_dis_file = f'{fbase}_ion{fext}'

    # default: use absolute path for all files
    inps.template_file    = os.path.abspath(inps.template_file)
    inps.iono_stack_file  = os.path.abspath(inps.iono_stack_file)
    inps.iono_file        = os.path.abspath(inps.iono_file)
    if inps.dis_file:
        inps.dis_file     = os.path.abspath(inps.dis_file)
        inps.cor_dis_file = os.path.abspath(inps.cor_dis_file)

    return inps


def read_template2inps(template_file, inps):
    """Read input template options into Namespace inps"""
    print('read input option from template file:', template_file)

    from mintpy.utils import ptime, readfile, utils1 as ut

    iDict = vars(inps)
    template = readfile.read_template(template_file, skip_chars=['[', ']'])
    template = ut.check_template_auto_value(template)

    key_list = [i for i in list(iDict.keys()) if key_prefix + i in template.keys()]
    for key in key_list:
        value = template[key_prefix + key]
        if value:
            if key == 'excludeDate':
                iDict[key] = ptime.yyyymmdd(value.split(','))
            elif key == 'excludeDate12':
                iDict[key] = ptime.yyyymmdd_date12(value.split(','))

    return inps


###############################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.iono_split_spectrum import run_iono_split_spectrum

    # run
    run_iono_split_spectrum(inps)

###############################################################
if __name__ == '__main__':
    main(sys.argv[1:])
