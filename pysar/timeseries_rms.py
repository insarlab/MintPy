#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2017, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################
#


import os
import sys
import argparse

import numpy as np
import matplotlib as mpl; mpl.use('Agg')
import matplotlib.pyplot as plt

import pysar._datetime as ptime
import pysar._network as pnet
import pysar._readfile as readfile
import pysar._pysar_utilities as ut


######################################################################################################
def read_template2inps(templateFile, inps=None):
    '''Update inps with pysar.residualRms.* option from templateFile'''
    if not inps:
        inps = cmdLineParse()

    template = readfile.read_template(templateFile)
    key_list = template.keys()

    prefix = 'pysar.residualRms.'

    key = prefix+'maskFile'
    if key in key_list:
        value = template[key]
        if value == 'auto':
            inps.mask_file = 'maskTempCoh.h5'
        elif value == 'no':
            inps.mask_file = None
        else:
            inps.mask_file = value

    key = prefix+'ramp'
    if key in key_list:
        value = template[key]
        if value == 'auto':
            inps.ramp_type = 'quadratic'
        else:
            inps.ramp_type = value

    key = prefix+'threshold'
    if key in key_list:
        value = template[key]
        if value == 'auto':
            inps.min_rms = 0.02
        else:
            inps.min_rms = float(value)

    key = prefix+'saveRefDate'
    if key in key_list:
        value = template[key]
        if value in ['auto','yes']:
            inps.save_reference_date = True
        else:
            inps.save_reference_date = False

    key = prefix+'saveExcludeDate'
    if key in key_list:
        value = template[key]
        if value in ['auto','yes']:
            inps.save_exclude_date = True
        else:
            inps.save_exclude_date = False

    return inps


######################################################################################################
TEMPLATE='''
## calculate the deramped Root Mean Square (RMS) for each epoch of timeseries residual from DEM error inversion
## To get rid of long wavelength component in space, a ramp is removed for each epoch.
pysar.residualRms.maskFile        = auto  #[file name / no], auto for maskTempCoh.h5, mask for ramp estimation
pysar.residualRms.ramp            = auto  #[quadratic / plane / no], auto for quadratic
pysar.residualRms.threshold       = auto  #[0.0-inf], auto for 0.02, minimum RMS in meter for exclude date(s)
pysar.residualRms.saveRefDate     = auto  #[yes / no], auto for yes, save date with min RMS to txt file.
pysar.residualRms.saveExcludeDate = auto  #[yes / no], auto for yes, save date(s) with RMS > threshold to txt file.
'''

EXAMPLE='''example:
  timeseries_rms.py  timeseriesResidual.h5 
  timeseries_rms.py  timeseriesResidual.h5  --template pysarApp_template.txt
  timeseries_rms.py  timeseriesResidual.h5  -m maskTempCoh.h5  --min-rms 0.03
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Calculate Root Mean Square (RMS) of deramped time series.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('timeseries_file', help='Timeseries file')
    parser.add_argument('--template','-t', dest='template_file',\
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
    inps = parser.parse_args()
    inps.save_reference_date = True
    inps.save_exclude_date = True
    return inps


######################################################################################################
def main(argv):
    inps = cmdLineParse()
    if inps.template_file:
        inps = read_template2inps(inps.template_file)

    ##### calculate timeseries of residual Root Mean Square
    #std_list, date_list = ut.get_residual_std(inps.timeseries_file, inps.mask_file, inps.ramp_type)
    rms_list, date_list = ut.get_residual_rms(inps.timeseries_file, inps.mask_file, inps.ramp_type)

    ##### reference_date.txt
    print '------------------------------------------------------------'
    ref_idx = np.argmin(rms_list)
    ref_date = date_list[ref_idx]
    print 'date with minimum residual RMS: %s - %.4f' % (ref_date, rms_list[ref_idx])

    refTxtFile = 'reference_date.txt'
    if (inps.save_reference_date and \
        ut.update_file(refTxtFile, [inps.timeseries_file, inps.mask_file, inps.template_file],\
                       check_readable=False)):
        f = open(refTxtFile, 'w')
        f.write(ref_date+'\n')
        f.close()
        print 'save date to file: '+refTxtFile

    ##### exclude_date.txt
    print '------------------------------------------------------------'
    ex_idx_list = [rms_list.index(i) for i in rms_list if i > inps.min_rms]
    print 'date(s) with residual RMS > '+str(inps.min_rms)
    exTxtFile = 'exclude_date.txt'
    if ex_idx_list:
        if (inps.save_exclude_date and \
            ut.update_file(exTxtFile, [inps.timeseries_file, inps.mask_file, inps.template_file],\
                           check_readable=False)):
            f = open(exTxtFile, 'w')
            for i in ex_idx_list:
                print '%s - %.4f' % (date_list[i], rms_list[i])
                f.write(date_list[i]+'\n')
            f.close()
            print 'save date(s) to file: '+exTxtFile
    else:
        print 'None.'

    ##### Plot
    fig_name = os.path.dirname(os.path.abspath(inps.timeseries_file))+\
               '/rms_'+os.path.splitext(inps.timeseries_file)[0]
    if inps.ramp_type != 'no':
        fig_name += '_'+inps.ramp_type
    fig_name += '.pdf'

    if ut.update_file(fig_name, [exTxtFile, refTxtFile, inps.template_file], check_readable=False):
        if inps.fig_size:
            fig = plt.figure(figsize=inps.fig_size)
        else:
            fig = plt.figure()
        ax = fig.add_subplot(111)
        font_size = 12

        dates, datevector = ptime.date_list2vector(date_list)
        try:    bar_width = ut.mode(np.diff(dates).tolist())*3/4
        except: bar_width = np.min(np.diff(dates).tolist())*3/4
        x_list = [i-bar_width/2 for i in dates]

        rms_list = [i*1000. for i in rms_list]
        min_rms = inps.min_rms * 1000.
        # Plot all dates
        ax.bar(x_list, rms_list, bar_width.days)
        #ax.bar(x_list, rms_list, bar_width.days)

        # Plot reference date
        #if inps.save_reference_date:
        ax.bar(x_list[ref_idx], rms_list[ref_idx], bar_width.days, label='Reference date')

        # Plot exclude dates
        #if ex_idx_list and inps.save_exclude_date:
        if ex_idx_list:
            ex_x_list = [x_list[i] for i in ex_idx_list]
            ex_rms_list = [rms_list[i] for i in ex_idx_list]
            ax.bar(ex_x_list, ex_rms_list, bar_width.days, color='darkgray', label='Exclude date(s)')

        # Plot min_rms line
        ax, xmin, xmax = ptime.auto_adjust_xaxis_date(ax, datevector, font_size, every_year=inps.tick_year_num)
        ax.plot(np.array([xmin, xmax]), np.array([min_rms, min_rms]), '--k')

        # axis format
        ax = pnet.auto_adjust_yaxis(ax, rms_list+[min_rms], font_size, ymin=0.0)
        ax.set_xlabel('Time [years]',fontsize=font_size)
        ax.set_ylabel('Root Mean Square [mm]',fontsize=font_size)
        ax.yaxis.set_ticks_position('both')
        ax.tick_params(labelsize=font_size)

        if inps.save_reference_date or inps.save_exclude_date:
            plt.legend(fontsize=font_size)

        # save figure
        fig.savefig(fig_name, bbox_inches='tight', transparent=True)
        print 'save figure to file: '+fig_name

    return


######################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])


