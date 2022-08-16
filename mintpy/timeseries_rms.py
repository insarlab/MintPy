#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2017                               #
############################################################


import os
import sys
import numpy as np

from mintpy.defaults.template import get_template_content
from mintpy.utils import readfile, utils as ut, plot as pp
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
    return inps


def read_template2inps(templateFile, inps):
    """Update inps with mintpy.residualRMS.* option from templateFile"""
    if not inps:
        inps = cmd_line_parse()
    inpsDict = vars(inps)
    print('read options from template file: '+os.path.basename(templateFile))
    template = readfile.read_template(templateFile)
    template = ut.check_template_auto_value(template)

    prefix = 'mintpy.residualRMS.'
    keyList = [i for i in list(inpsDict.keys()) if prefix+i in template.keys()]
    for key in keyList:
        value = template[prefix+key]
        if value:
            if key in ['maskFile', 'deramp']:
                inpsDict[key] = value
            elif key in ['cutoff']:
                inpsDict[key] = float(value)
    return inps


def analyze_rms(date_list, rms_list, inps):
    # reference date
    ref_idx = np.argmin(rms_list)
    print('-'*50+'\ndate with min RMS: {} - {:.4f}'.format(date_list[ref_idx],
                                                           rms_list[ref_idx]))
    ref_date_file = 'reference_date.txt'
    if ut.run_or_skip(out_file=ref_date_file,
                      in_file=[inps.timeseries_file, inps.maskFile, inps.template_file],
                      check_readable=False) == 'run':
        with open(ref_date_file, 'w') as f:
            f.write(date_list[ref_idx]+'\n')
        print('save date to file: '+ref_date_file)

    # exclude date(s) - outliers
    # equivalent calculation using numpy assuming Gaussian distribution as:
    # rms_threshold = np.median(rms_list) / .6745 * inps.cutoff
    rms_threshold = ut.median_abs_deviation_threshold(rms_list, center=0., cutoff=inps.cutoff)

    ex_idx = [rms_list.index(i) for i in rms_list if i > rms_threshold]
    print('-'*50)
    print(f'date(s) with RMS > {inps.cutoff} * median RMS ({rms_threshold:.4f})')

    ex_date_file = 'exclude_date.txt'
    if ex_idx:
        # print
        for i in ex_idx:
            print('{} - {:.4f}'.format(date_list[i], rms_list[i]))
        # save to text file
        with open(ex_date_file, 'w') as f:
            for i in ex_idx:
                f.write(date_list[i]+'\n')
        print('save date(s) to file: '+ex_date_file)
    else:
        print('None.')
        if os.path.isfile(ex_date_file):
            os.remove(ex_date_file)

    return inps


######################################################################################################
def main(iargs=None):

    # read inputs
    inps = cmd_line_parse(iargs)
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    # calculate timeseries of residual Root Mean Square
    inps.rms_list, inps.date_list, inps.rms_file = ut.get_residual_rms(
        inps.timeseries_file,
        mask_file=inps.maskFile,
        ramp_type=inps.deramp,
    )

    # analyze RMS: generate reference/exclude_date.txt files
    analyze_rms(inps.date_list, inps.rms_list, inps)

    # plot RMS
    pp.plot_timeseries_rms(
        rms_file=inps.rms_file,
        cutoff=inps.cutoff,
        out_fig=os.path.splitext(inps.rms_file)[0]+'.pdf',
        disp_fig=False,
        fig_size=inps.fig_size,
        tick_year_num=inps.tick_year_num,
    )

    return


######################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
