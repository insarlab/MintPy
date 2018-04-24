#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2017, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################


import os
import sys
import argparse
import numpy as np
import matplotlib as mpl; mpl.use('Agg')
import matplotlib.pyplot as plt
from pysar.utils import readfile, datetime as ptime, utils as ut, plot as pp


######################################################################################################
def read_template2inps(templateFile, inps=None):
    """Update inps with pysar.residualRms.* option from templateFile"""
    if not inps:
        inps = cmd_line_parse()

    template = readfile.read_template(templateFile)
    prefix = 'pysar.residualRms.'

    key = prefix+'maskFile'
    if key in template.keys():
        value = template[key]
        if value == 'auto':
            inps.mask_file = 'maskTempCoh.h5'
        elif value == 'no':
            inps.mask_file = None
        else:
            inps.mask_file = value

    key = prefix+'ramp'
    if key in template.keys():
        value = template[key]
        if value == 'auto':
            inps.ramp_type = 'quadratic'
        else:
            inps.ramp_type = value

    key = prefix+'threshold'
    if key in template.keys():
        value = template[key]
        if value == 'auto':
            inps.min_rms = 0.02
        else:
            inps.min_rms = float(value)

    return inps


######################################################################################################
TEMPLATE = """
## calculate the deramped Root Mean Square (RMS) for each epoch of timeseries residual from DEM error inversion
## To get rid of long wavelength component in space, a ramp is removed for each epoch.
pysar.residualRms.maskFile        = auto  #[file name / no], auto for maskTempCoh.h5, mask for ramp estimation
pysar.residualRms.ramp            = auto  #[quadratic / plane / no], auto for quadratic
pysar.residualRms.threshold       = auto  #[0.0-inf], auto for 0.02, minimum RMS in meter for exclude date(s)
"""

EXAMPLE = """example:
  timeseries_rms.py  timeseriesResidual.h5 
  timeseries_rms.py  timeseriesResidual.h5  --template pysarApp_template.txt
  timeseries_rms.py  timeseriesResidual.h5  -m maskTempCoh.h5  --min-rms 0.03
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Calculate Root Mean Square (RMS) of deramped time series.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('timeseries_file', help='Timeseries file')
    parser.add_argument('-t','--template', dest='template_file',\
                        help='template file with options below:\n'+TEMPLATE+'\n')
    parser.add_argument('-m','--mask', dest='mask_file', default='maskTempCoh.h5',\
                        help='mask file for estimation')
    parser.add_argument('-s', dest='ramp_type', default='quadratic',\
                        help='ramp type to be remove for RMS calculation.\n'+\
                             'default - quadratic; no - do not remove ramp')
    parser.add_argument('--min-rms', dest='min_rms', default='0.02', type=float,\
                        help='minimum RMS in m, threshold used to exclude dates, default: 0.02 m')
    parser.add_argument('--figsize', dest='fig_size', metavar=('WID','LEN'), type=float, nargs=2,\
                        help='figure size in inches - width and length')
    parser.add_argument('--tick-year-num', dest='tick_year_num', type=int, default=1, help='Year number per major tick')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


def save_date2txt_file(inps):
    inps.refDateFile = 'reference_date.txt'
    if ut.update_file(inps.refDateFile, [inps.timeseries_file, inps.mask_file, inps.template_file], check_readable=False):
        f = open(inps.refDateFile, 'w')
        f.write(inps.refDate+'\n')
        f.close()
        print('save date to file: '+inps.refDateFile)

    inps.exDateFile = 'exclude_date.txt'
    if inps.exIdxList and ut.update_file(inps.exDateFile, [inps.timeseries_file, inps.mask_file, inps.template_file],\
                                         check_readable=False):
        f = open(inps.exDateFile, 'w')
        for i in inps.exIdxList:
            f.write(inps.dateList[i]+'\n')
        f.close()
        print('save date(s) to file: '+inps.exDateFile)
    else:
        print('None.')
    return inps


def plot_bar4date_rms(inps):
    inps.figName = os.path.splitext(inps.rmsFile)[0]+'.pdf'
    if ut.update_file(inps.figName, [inps.exDateFile, inps.refDateFile, inps.template_file], check_readable=False):
        if inps.fig_size:
            fig = plt.figure(figsize=inps.fig_size)
        else:
            fig = plt.figure()
        ax = fig.add_subplot(111)
        font_size = 12

        dates, datevector = ptime.date_list2vector(inps.dateList)
        try:    bar_width = ut.most_common(np.diff(dates).tolist())*3/4
        except: bar_width = np.min(np.diff(dates).tolist())*3/4
        x_list = [i-bar_width/2 for i in dates]

        inps.rmsList = [i*1000. for i in inps.rmsList]
        min_rms = inps.min_rms * 1000.
        # Plot all dates
        ax.bar(x_list, inps.rmsList, bar_width.days)

        # Plot reference date
        ax.bar(x_list[inps.refDateIndex], inps.rmsList[inps.refDateIndex], bar_width.days, label='Reference date')

        # Plot exclude dates
        if inps.exIdxList:
            ex_x_list = [x_list[i] for i in inps.exIdxList]
            inps.exRmsList = [inps.rmsList[i] for i in inps.exIdxList]
            ax.bar(ex_x_list, inps.exRmsList, bar_width.days, color='darkgray', label='Exclude date(s)')

        # Plot min_rms line
        ax, xmin, xmax = pp.auto_adjust_xaxis_date(ax, datevector, font_size, every_year=inps.tick_year_num)
        ax.plot(np.array([xmin, xmax]), np.array([min_rms, min_rms]), '--k')

        # axis format
        ax = pp.auto_adjust_yaxis(ax, inps.rmsList+[min_rms], font_size, ymin=0.0)
        ax.set_xlabel('Time [years]',fontsize=font_size)
        ax.set_ylabel('Root Mean Square [mm]',fontsize=font_size)
        ax.yaxis.set_ticks_position('both')
        ax.tick_params(labelsize=font_size)
        plt.legend(fontsize=font_size)

        # save figure
        fig.savefig(inps.figName, bbox_inches='tight', transparent=True)
        print('save figure to file: '+inps.figName)
    return inps


######################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    if inps.template_file:
        inps = read_template2inps(inps.template_file)

    ##### calculate timeseries of residual Root Mean Square
    inps.rmsList, inps.dateList, inps.rmsFile = ut.get_residual_rms(inps.timeseries_file, inps.mask_file, inps.ramp_type)

    ##### reference date
    inps.refDateIndex = np.argmin(inps.rmsList)
    inps.refDate = inps.dateList[inps.refDateIndex]
    print('-'*50)
    print('date with minimum residual RMS: %s - %.4f' % (inps.refDate, inps.rmsList[inps.refDateIndex]))

    ##### exclude date(s)
    inps.exIdxList = [inps.rmsList.index(i) for i in inps.rmsList if i > inps.min_rms]
    print('-'*50)
    print('date(s) with residual RMS > {}'.format(inps.min_rms))
    if inps.exIdxList:
        for i in inps.exIdxList:
            print('%s - %.4f' % (inps.dateList[i], inps.rmsList[i]))

    ##### Save to text file
    inps = save_date2txt_file(inps)

    ##### Plot
    inps = plot_bar4date_rms(inps)

    return


######################################################################################################
if __name__ == '__main__':
    main()


