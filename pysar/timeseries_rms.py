#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2017-2018, Zhang Yunjun                     #
# Author:  Zhang Yunjun                                    #
############################################################


import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pysar.utils import readfile, ptime, utils as ut, plot as pp


######################################################################################################
TEMPLATE = """
## calculate the deramped Root Mean Square (RMS) for each epoch of timeseries residual
## To get rid of long wavelength component in space, a ramp is removed for each epoch.
## Set optimal reference date to date with min RMS
## Set exclude dates (outliers) to dates with RMS > cutoff * median RMS (Median Absolute Deviation)
pysar.residualRms.maskFile = auto  #[file name / no], auto for maskTempCoh.h5, mask for ramp estimation
pysar.residualRms.ramp     = auto  #[quadratic / plane / no], auto for quadratic
pysar.residualRms.cutoff   = auto  #[0.0-inf], auto for 3
"""

EXAMPLE = """example:
  timeseries_rms.py  timeseriesResidual.h5 
  timeseries_rms.py  timeseriesResidual.h5  --template pysarApp_template.txt
  timeseries_rms.py  timeseriesResidual.h5  -m maskTempCoh.h5  --cutoff 3
"""

REFERENCE="""reference:
Rousseeuw, P. J., and M. Hubert (2011), Robust statistics for outlier detection,
    Wiley Interdisciplinary Reviews: Data Mining and Knowledge Discovery, 1(1),
    73-79, doi:doi:10.1002/widm.2.
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Calculate Root Mean Square (RMS) of deramped time series.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('timeseries_file', help='Timeseries file')
    parser.add_argument('-t', '--template', dest='template_file',
                        help='template file with options below:\n'+TEMPLATE+'\n')
    parser.add_argument('-m', '--mask', dest='mask_file', default='maskTempCoh.h5',
                        help='mask file for estimation')
    parser.add_argument('-s', dest='ramp_type', default='quadratic',
                        help='ramp type to be remove for RMS calculation.\n' +
                             'default - quadratic; no - do not remove ramp')
    parser.add_argument('--cutoff', dest='cutoff', default='3', type=float,
                        help='M-score used for outlier detection based on standardised residuals\n'+
                             'Recommend range: [3, 4], default is 3.')
    parser.add_argument('--figsize', dest='fig_size', metavar=('WID', 'LEN'), type=float, nargs=2,
                        help='figure size in inches - width and length')
    parser.add_argument('--tick-year-num', dest='tick_year_num',
                        type=int, default=1, help='Year number per major tick')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


def read_template2inps(templateFile, inps=None):
    """Update inps with pysar.residualRms.* option from templateFile"""
    if not inps:
        inps = cmd_line_parse()
    inpsDict = vars(inps)
    print('read options from template file: '+os.path.basename(templateFile))
    template = readfile.read_template(templateFile)
    template = ut.check_template_auto_value(template)

    prefix = 'pysar.residualRms.'
    keyList = [i for i in list(inpsDict.keys()) if prefix+i in template.keys()]
    for key in keyList:
        value = template[prefix+key]
        if value:
            if key in ['maskFile', 'ramp']:
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
    if ut.update_file(ref_date_file, [inps.timeseries_file,
                                      inps.mask_file,
                                      inps.template_file], check_readable=False):
        with open(ref_date_file, 'w') as f:
            f.write(date_list[ref_idx]+'\n')
        print('save date to file: '+ref_date_file)

    # exclude date(s) - outliers
    try:
        rms_threshold = median_abs_deviation_threshold(rms_list,
                                                       center=0.,
                                                       kappa=inps.cutoff)
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
            rmCmd = 'rm {}'.format(ex_date_file)
            print(rmCmd)
            os.system(rmCmd)

    # plot bar figure and save
    fig_file = os.path.splitext(inps.rms_file)[0]+'.pdf'
    if ut.update_file(fig_file, [ex_date_file,
                                 ref_date_file,
                                 inps.template_file], check_readable=False):
        fig, ax = plt.subplots(figsize=inps.fig_size)
        ax = plot_rms_bar(ax, date_list, rms_list, rms_threshold, cutoff=inps.cutoff)
        fig.savefig(fig_file, bbox_inches='tight', transparent=True)
        print('save figure to file: '+fig_file)
    return inps


def median_abs_deviation_threshold(data, center=0., kappa=3.):
    """calculate rms_threshold based on the standardised residual
    outlier detection with median absolute deviation.
    
    With the default input arguments, it's equivalent to:
        np.median(data) / .6745 * 3.0
    """
    from statsmodels.robust import mad
    rms_mad = mad(data, c=0.67448975019608171, center=center)
    rms_threshold = center + kappa * rms_mad
    return rms_threshold


def plot_rms_bar(ax, date_list, rms_list, rms_threshold,
                 unit_scale=1000., font_size=12,
                 tick_year_num=1, legend_loc='best', cutoff=3.):
    """
        legend_loc - 'upper right' or (0.5, 0.5)
    """
    dates, datevector = ptime.date_list2vector(date_list)
    try:
        bar_width = min(ut.most_common(np.diff(dates).tolist(), k=2))*3/4
    except:
        bar_width = np.min(np.diff(dates).tolist())*3/4
    datex = np.array(dates) - bar_width / 2
    rms = np.array(rms_list)

    # Plot all dates
    ax.bar(datex, rms * unit_scale, bar_width.days, color=pp.mplColors[0])

    # Plot reference date
    ref_idx = np.argmin(rms)
    ax.bar(datex[ref_idx],
           rms[ref_idx] * unit_scale,
           bar_width.days,
           color=pp.mplColors[1],
           label='Reference date')

    # Plot exclude dates
    ex_idx = rms > rms_threshold
    if not np.all(ex_idx==False):
        ax.bar(datex[ex_idx],
               rms[ex_idx] * unit_scale,
               bar_width.days,
               color='darkgray',
               label='Exclude date')

    # Plot rms_threshold line
    (ax, xmin, xmax) = pp.auto_adjust_xaxis_date(ax, datevector, font_size,
                                                 every_year=tick_year_num)
    ax.plot(np.array([xmin, xmax]),
            np.array([rms_threshold, rms_threshold]) * unit_scale,
            '--k',
            label='RMS Threshold')

    # axis format
    ax = pp.auto_adjust_yaxis(ax, np.append(rms, rms_threshold) * unit_scale,
                              font_size, ymin=0.0)
    ax.set_xlabel('Time [years]', fontsize=font_size)
    ax.set_ylabel('Phase residual RMS [mm]', fontsize=font_size)
    #ax.yaxis.set_ticks_position('both')
    ax.tick_params(which='both', direction='in', labelsize=font_size,
                   bottom=True, top=True, left=True, right=True)

    # 2nd axes for circles
    divider = make_axes_locatable(ax)
    ax2 = divider.append_axes("right", "10%", pad="2%")
    ax2.plot(np.ones(rms.shape, np.float32) * 0.5,
             rms * unit_scale,
             'o', mfc='none',
             color=pp.mplColors[0])
    ax2.plot(np.ones(rms.shape, np.float32)[ref_idx] * 0.5,
             rms[ref_idx] * unit_scale,
             'o', mfc='none',
             color=pp.mplColors[1])
    if not np.all(ex_idx==False):
        ax2.plot(np.ones(rms.shape, np.float32)[ex_idx] * 0.5,
                 rms[ex_idx] * unit_scale,
                 'o', mfc='none',
                 color='darkgray')
    ax2.plot(np.array([0, 1]),
             np.array([rms_threshold, rms_threshold]) * unit_scale,
             '--k')

    ax2.set_ylim(ax.get_ylim())
    ax2.set_xlim([0, 1])
    ax2.tick_params(which='both', direction='in', labelsize=font_size,
                    bottom=True, top=True, left=True, right=True)
    ax2.get_xaxis().set_ticks([])
    ax2.get_yaxis().set_ticklabels([])

    ax.legend(loc=legend_loc, fontsize=font_size)
    ymin, ymax = ax.get_ylim()
    ax.annotate('median RMS * {}'.format(cutoff),
                xy=(xmin + (xmax-xmin)*0.05, rms_threshold * unit_scale - (ymax-ymin)*0.1),
                color='k', xycoords='data', fontsize=font_size)
    return ax


######################################################################################################
def main(iargs=None):
    plt.switch_backend('Agg')  # Backend setting

    inps = cmd_line_parse(iargs)
    if inps.template_file:
        inps = read_template2inps(inps.template_file)

    # calculate timeseries of residual Root Mean Square
    (inps.rms_list,
     inps.date_list,
     inps.rms_file) = ut.get_residual_rms(inps.timeseries_file,
                                          inps.mask_file,
                                          inps.ramp_type)

    analyze_rms(inps.date_list, inps.rms_list, inps)
    return


######################################################################################################
if __name__ == '__main__':
    main()
