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
    '''Update inps with pysar.residualStd.* option from templateFile'''
    if not inps:
        inps = cmdLineParse()

    template = readfile.read_template(templateFile)
    key_list = template.keys()

    prefix = 'pysar.residualStd.'

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
            inps.min_std = 0.02
        else:
            inps.min_std = float(value)

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
## calculate the deramped standard deviation (STD) for each epoch of timeseries residual from DEM error inversion
## To get rid of long wavelength component in space, a ramp is removed for each epoch.
pysar.residualStd.maskFile        = auto  #[file name / no], auto for maskTempCoh.h5, mask for ramp estimation
pysar.residualStd.ramp            = auto  #[quadratic / plane / no], auto for quadratic
pysar.residualStd.threshold       = auto  #[0.0-inf], auto for 0.02, minimum STD in meter for exclude date(s)
pysar.residualStd.saveRefDate     = auto  #[yes / no], auto for yes, save date with min RSD to txt/pdf file.
pysar.residualStd.saveExcludeDate = auto  #[yes / no], auto for yes, save date(s) with RSD > minStd to txt/pdf file.
'''

EXAMPLE='''example:
  timeseries_std.py  timeseries_ECMWF_demErrInvResid.h5 
  timeseries_std.py  timeseries_ECMWF_demErrInvResid.h5  --template pysarApp_template.txt
  timeseries_std.py  timeseries_ECMWF_demErrInvResid.h5  -m maskTempCoh.h5  --min-std 0.03
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Calculate standard deviation of deramped phase residual for time series.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('timeseries_file', help='Timeseries file')
    parser.add_argument('--template', dest='template_file', help='template file with options below:\n'+TEMPLATE+'\n')
    parser.add_argument('-m', dest='mask_file', default='maskTempCoh.h5', help='mask file for estimation')
    parser.add_argument('-s', dest='ramp_type', default='quadratic',\
                        help='ramp type to be remove for STD calculation.\n'+\
                             'default - quadratic; no - do not remove ramp')
    parser.add_argument('--min-std', dest='min_std', default='0.02', type=float,\
                        help='minimum standard deviation in m, threshold used to exclude dates, default: 0.02 m')
    inps = parser.parse_args()
    inps.save_reference_date = True
    inps.save_exclude_date = True
    return inps


######################################################################################################
def main(argv):
    inps = cmdLineParse()
    if inps.template_file:
        inps = read_template2inps(inps.template_file)

    ##### calculate timeseries of residual standard deviation
    std_list, date_list = ut.get_residual_std(inps.timeseries_file, inps.mask_file, inps.ramp_type)

    ##### reference_date.txt
    print '------------------------------------------------------------'
    ref_idx = np.argmin(std_list)
    ref_date = date_list[ref_idx]
    print 'date with minimum residual std: %s - %.4f' % (ref_date, std_list[ref_idx])

    if inps.save_reference_date:
        txtFile = 'reference_date.txt'
        f = open(txtFile, 'w')
        f.write(ref_date+'\n')
        f.close()
        print 'save date to file: '+txtFile

    ##### exclude_date.txt
    print '------------------------------------------------------------'
    ex_idx_list = [std_list.index(i) for i in std_list if i > inps.min_std]
    print 'date(s) with residual std > '+str(inps.min_std)

    if ex_idx_list:
        if inps.save_exclude_date:
            txtFile = 'exclude_date.txt'
            f = open(txtFile, 'w')
            for i in ex_idx_list:
                print '%s - %.4f' % (date_list[i], std_list[i])
                f.write(date_list[i]+'\n')
            f.close()
            print 'save date(s) to file: '+txtFile
    else:
        print 'None.'

    ##### Plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    font_size = 12

    dates, datevector = ptime.date_list2vector(date_list)
    bar_width = np.min(np.diff(dates))*3/4
    x_list = [i-bar_width/2 for i in dates]

    # Plot all dates
    ax.bar(x_list, std_list, bar_width.days)

    # Plot reference date
    if inps.save_reference_date:
        ax.bar(x_list[ref_idx], std_list[ref_idx], bar_width.days, label='Reference date')

    # Plot exclude dates
    if ex_idx_list and inps.save_exclude_date:
        ex_x_list = [x_list[i] for i in ex_idx_list]
        ex_std_list = [std_list[i] for i in ex_idx_list]
        ax.bar(ex_x_list, ex_std_list, bar_width.days, color='darkgray', label='Exclude date(s)')

    # Plot min_std line
    ax, xmin, xmax = ptime.auto_adjust_xaxis_date(ax, datevector, font_size)
    ax.plot(np.array([xmin, xmax]), np.array([inps.min_std, inps.min_std]), '-')

    # axis format
    ax = pnet.auto_adjust_yaxis(ax, std_list+[inps.min_std], font_size, ymin=0.0)
    ax.set_xlabel('Time [years]',fontsize=font_size)
    ax.set_ylabel('Standard Deviation (m)',fontsize=font_size)
    ax.yaxis.set_ticks_position('both')

    if inps.save_reference_date or inps.save_exclude_date:
        plt.legend()

    # save figure
    if inps.ramp_type != 'no':
        fig_name = os.path.splitext(inps.timeseries_file)[0]+'_'+inps.ramp_type+'_std.pdf'
    else:
        fig_name = os.path.splitext(inps.timeseries_file)[0]+'_std.pdf'
    fig.savefig(fig_name, bbox_inches='tight')
    print 'save figure to file: '+fig_name

    return


######################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])


