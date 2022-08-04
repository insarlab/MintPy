#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2017                               #
############################################################


import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mintpy.defaults.template import get_template_content
from mintpy.utils import readfile, ptime, utils as ut, plot as pp
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
    try:
        rms_threshold = ut.median_abs_deviation_threshold(rms_list, center=0., cutoff=inps.cutoff)
    except:
        # equivalent calculation using numpy assuming Gaussian distribution
        rms_threshold = np.median(rms_list) / .6745 * inps.cutoff

    ex_idx = [rms_list.index(i) for i in rms_list if i > rms_threshold]
    print(('-'*50+'\ndate(s) with RMS > {} * median RMS'
           ' ({:.4f})'.format(inps.cutoff, rms_threshold)))
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

    # plot bar figure and save
    fig_file = os.path.splitext(inps.rms_file)[0]+'.pdf'
    fig, ax = plt.subplots(figsize=inps.fig_size)
    print('create figure in size:', inps.fig_size)
    ax = plot_rms_bar(ax, date_list, np.array(rms_list)*1000., cutoff=inps.cutoff)
    fig.savefig(fig_file, bbox_inches='tight', transparent=True)
    print('save figure to file: '+fig_file)
    return inps


def plot_rms_bar(ax, date_list, rms, cutoff=3., font_size=12,
                 tick_year_num=1, legend_loc='best',
                 disp_legend=True, disp_side_plot=True, disp_thres_text=False,
                 ylabel='Residual phase RMS [mm]'):
    """ Bar plot Phase Residual RMS
    Parameters: ax : Axes object
                date_list : list of string in YYYYMMDD format
                rms    : 1D np.array of float for RMS value in mm
                cutoff : cutoff value of MAD outlier detection
                tick_year_num : int, number of years per major tick
                legend_loc : 'upper right' or (0.5, 0.5)
    Returns:    ax : Axes object
    """
    dates, datevector = ptime.date_list2vector(date_list)
    dates = np.array(dates)
    try:
        bar_width = min(ut.most_common(np.diff(dates).tolist(), k=2))*3/4
    except:
        bar_width = np.min(np.diff(dates).tolist())*3/4
    rms = np.array(rms)

    # Plot all dates
    ax.bar(dates, rms, bar_width.days, color=pp.mplColors[0])

    # Plot reference date
    ref_idx = np.argmin(rms)
    ax.bar(dates[ref_idx], rms[ref_idx], bar_width.days, color=pp.mplColors[1], label='Reference date')

    # Plot exclude dates
    rms_threshold = ut.median_abs_deviation_threshold(rms, center=0., cutoff=cutoff)
    ex_idx = rms > rms_threshold
    if not np.all(ex_idx==False):
        ax.bar(dates[ex_idx], rms[ex_idx], bar_width.days, color='darkgray', label='Exclude date')

    # Plot rms_threshold line
    (ax, xmin, xmax) = pp.auto_adjust_xaxis_date(ax, datevector, font_size, every_year=tick_year_num)
    ax.plot(np.array([xmin, xmax]), np.array([rms_threshold, rms_threshold]), '--k',
            label='Median Abs Dev * {}'.format(cutoff))

    # axis format
    ax = pp.auto_adjust_yaxis(ax, np.append(rms, rms_threshold), font_size, ymin=0.0)
    #ax.set_xlabel('Time [years]', fontsize=font_size)
    ax.set_ylabel(ylabel, fontsize=font_size)
    ax.tick_params(which='both', direction='in', labelsize=font_size,
                   bottom=True, top=True, left=True, right=True)

    # 2nd axes for circles
    if disp_side_plot:
        divider = make_axes_locatable(ax)
        ax2 = divider.append_axes("right", "10%", pad="2%")
        ax2.plot(np.ones(rms.shape, np.float32) * 0.5, rms, 'o', mfc='none', color=pp.mplColors[0])
        ax2.plot(np.ones(rms.shape, np.float32)[ref_idx] * 0.5, rms[ref_idx], 'o', mfc='none', color=pp.mplColors[1])
        if not np.all(ex_idx==False):
            ax2.plot(np.ones(rms.shape, np.float32)[ex_idx] * 0.5, rms[ex_idx], 'o', mfc='none', color='darkgray')
        ax2.plot(np.array([0, 1]), np.array([rms_threshold, rms_threshold]), '--k')

        ax2.set_ylim(ax.get_ylim())
        ax2.set_xlim([0, 1])
        ax2.tick_params(which='both', direction='in', labelsize=font_size,
                        bottom=True, top=True, left=True, right=True)
        ax2.get_xaxis().set_ticks([])
        ax2.get_yaxis().set_ticklabels([])

    if disp_legend:
        ax.legend(loc=legend_loc, frameon=False, fontsize=font_size)

    # rms_threshold text
    if disp_thres_text:
        ymin, ymax = ax.get_ylim()
        yoff = (ymax - ymin) * 0.1
        if (rms_threshold - ymin) > 0.5 * (ymax - ymin):
            yoff *= -1.
        ax.annotate('Median Abs Dev * {}'.format(cutoff),
                    xy=(xmin + (xmax-xmin)*0.05, rms_threshold + yoff ),
                    color='k', xycoords='data', fontsize=font_size)
    return ax


######################################################################################################
def main(iargs=None):
    plt.switch_backend('Agg')  # Backend setting

    inps = cmd_line_parse(iargs)
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    # calculate timeseries of residual Root Mean Square
    (inps.rms_list,
     inps.date_list,
     inps.rms_file) = ut.get_residual_rms(inps.timeseries_file,
                                          mask_file=inps.maskFile,
                                          ramp_type=inps.deramp)

    analyze_rms(inps.date_list, inps.rms_list, inps)
    return


######################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
