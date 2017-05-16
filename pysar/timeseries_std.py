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

    template_dict = readfile.read_template(templateFile)
    key_list = template_dict.keys()

    prefix = 'pysar.residualStd.'
    if prefix+'maskFile' in key_list:  inps.mask_file = template_dict[prefix+'maskFile']
    if prefix+'ramp'     in key_list:  inps.ramp_type = template_dict[prefix+'ramp']
    if prefix+'minStd'   in key_list:  inps.min_std   = template_dict[prefix+'minStd']

    return inps


######################################################################################################
TEMPLATE='''
pysar.residualStd.maskFile = maskTempCoh_aoi.h5
pysar.residualStd.ramp     = quadratic
pysar.residualStd.minStd   = 0.02
'''

EXAMPLE='''example:
  timeseries_std.py  timeseries_ECMWF_demErrInvResid.h5 
  timeseries_std.py  timeseries_ECMWF_demErrInvResid.h5  --template KujuAlosAT422F650.template
  timeseries_std.py  timeseries_ECMWF_demErrInvResid.h5  -m maskTempCoh.h5  --min-std 0.03
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Calculate standard deviation of deramped phase residual for time series.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('timeseries_file', help='Timeseries file')
    parser.add_argument('--template', dest='template_file', help='template file with options below:\n'+TEMPLATE+'\n')
    parser.add_argument('-m', dest='mask_file', default='maskTempCoh.h5', help='mask file for estimation')
    parser.add_argument('-s', dest='ramp_type', default='quadratic', help='ramp type to be remove for STD calculation.')
    parser.add_argument('--min-std', dest='min_std', default='0.02', type=float,\
                        help='minimum standard deviation in m, threshold used to exclude dates, default: 0.02 m')
    inps = parser.parse_args()
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
    ax.bar(x_list[ref_idx], std_list[ref_idx], bar_width.days, label='Reference date')

    # Plot exclude dates
    if ex_idx_list:
        ex_x_list = [x_list[i] for i in ex_idx_list]
        ex_std_list = [std_list[i] for i in ex_idx_list]
        ax.bar(ex_x_list, ex_std_list, bar_width.days, color='darkgray', label='Exclude date')

    # Plot min_std line
    ax, xmin, xmax = ptime.auto_adjust_xaxis_date(ax, datevector, font_size)
    ax.plot(np.array([xmin, xmax]), np.array([inps.min_std, inps.min_std]), '-')

    # axis format
    ax = pnet.auto_adjust_yaxis(ax, std_list+[inps.min_std], font_size, ymin=0.0)
    ax.set_xlabel('Time [years]',fontsize=font_size)
    ax.set_ylabel('Standard Deviation (m)',fontsize=font_size)

    plt.legend()

    # save figure
    fig_name = os.path.splitext(inps.timeseries_file)[0]+'_'+inps.ramp_type+'_std.pdf'
    fig.savefig(fig_name, bbox_inches='tight')
    print 'save figure to file: '+fig_name

    return


######################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])


