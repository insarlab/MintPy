#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Aug 2015: Add -m/M/d option
# Yunjun, Jun 2016: Add -t option
# Yunjun, Aug 2015: Support drop_date txt file input


import os
import sys
import time
import datetime
import argparse

import numpy as np
import h5py

import pysar._datetime as ptime
import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar._pysar_utilities as ut


############################################################################
def get_exclude_date(inps, date_list_all):
    '''Get inps.ex_date full list
    Inputs:
        inps          - Namespace, 
        date_list_all - list of string for all available date in YYYYMMDD format
    Output:
        inps.ex_date  - list of string for exclude date in YYYYMMDD format
    '''
    yy_list_all = ptime.yyyymmdd2years(date_list_all)

    # 1. template_file
    if inps.template_file:
        print 'read option from template file: '+inps.template_file
        inps = read_template2inps(inps.template_file, inps)

    # 2. ex_date
    input_ex_date = list(inps.ex_date)
    inps.ex_date = []
    if input_ex_date:
        for ex_date in input_ex_date:
            if os.path.isfile(ex_date):
                ex_date = ptime.read_date_list(ex_date)
            else:
                ex_date = [ptime.yyyymmdd(ex_date)]
            inps.ex_date += list(set(ex_date) - set(inps.ex_date))
        # delete dates not existed in input file
        inps.ex_date = list(set(inps.ex_date).intersection(date_list_all))
        print 'exclude date:'+str(inps.ex_date)

    # 3. min_date
    if inps.min_date:
        print 'start date: '+inps.min_date
        yy_min = ptime.yyyymmdd2years(ptime.yyyymmdd(inps.min_date))
        for i in range(len(date_list_all)):
            date = date_list_all[i]
            if yy_list_all[i] < yy_min and date not in inps.ex_date:
                print '  remove date: '+date
                inps.ex_date.append(date)

    # 4. max_date
    if inps.max_date:
        print 'end date: '+inps.max_date
        yy_max = ptime.yyyymmdd2years(ptime.yyyymmdd(inps.max_date))
        for i in range(len(date_list_all)):
            date = date_list_all[i]
            if yy_list_all[i] > yy_max and date not in inps.ex_date:
                print '  remove date: '+date
                inps.ex_date.append(date)

    return inps.ex_date


def get_velocity_filename(timeseries_file, template_file=None, vel_file='velocity.h5', inps=None):
    '''Get output velocity filename
    Example: velocity_file = get_output_filename('timeseries_ECMWF_demErr_refDate.h5', 'KujuAlosAT422F650.template')
    '''
    if not inps:
        inps = cmdLineParse()
    inps.template_file = template_file

    h5 = h5py.File(timeseries_file, 'r')
    date_list_all = sorted(h5['timeseries'].keys())
    h5.close()

    ex_date = get_exclude_date(inps, date_list_all)
    if ex_date:
        vel_file = os.path.splitext(vel_file)[0]+'Ex.h5'
    return vel_file


def read_template2inps(template_file, inps=None):
    '''Read input template file into inps.ex_date'''
    if not inps:
        inps = cmdLineParse()
    template = readfile.read_template(template_file)
    key_list = template.keys()

    # Read template option
    prefix = 'pysar.velocity.'
    key = prefix+'excludeDate'
    if key in key_list:
        value = template[key]
        if value == 'auto':
            inps.ex_date = ['exclude_date.txt']
        elif value == 'no':
            inps.ex_date = []
        else:
            inps.ex_date = value.replace(',',' ').split()

    key = prefix+'startDate'
    if key in key_list:
        value = template[key]
        if value not in ['auto','no']:
            inps.min_date = ptime.yyyymmdd(value)

    key = prefix+'endDate'
    if key in key_list:
        value = template[key]
        if value not in ['auto','no']:
            inps.max_date = ptime.yyyymmdd(value)

    return inps


############################################################################
EXAMPLE='''example:
  timeseries2velocity.py  timeSeries_ECMWF_demCor.h5
  timeseries2velocity.py  timeseries_ECMWF_demCor_plane.h5  --template KyushuT73F2980_2990AlosD.template
  timeseries2velocity.py  timeseries.h5  --start-date 20080201
  timeseries2velocity.py  timeseries.h5  --start-date 20080201  --end-date 20100508
  timeseries2velocity.py  timeseries.h5  --exclude-date 20040502 20060708 20090103
  timeseries2velocity.py  timeseries.h5  --exclude-date exclude_date.txt
'''

TEMPLATE='''
## estimate linear velocity from timeseries, and from tropospheric delay file if exists.
pysar.velocity.excludeDate = auto   #[exclude_date.txt / 20080520,20090817 / no], auto for exclude_date.txt
pysar.velocity.startDate   = auto   #[20070101 / no], auto for no
pysar.velocity.endDate     = auto   #[20101230 / no], auto for no
'''

DROP_DATE_TXT='''exclude_date.txt:
20040502
20060708
20090103
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Inverse velocity from time series.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=TEMPLATE+'\n'+EXAMPLE)

    parser.add_argument('timeseries_file', help='Time series file for velocity inversion.')
    parser.add_argument('--start-date', dest='min_date', help='start date for velocity estimation')
    parser.add_argument('--end-date', dest='max_date', help='end date for velocity estimation')
    parser.add_argument('--exclude','--ex', dest='ex_date', nargs='+',\
                        help='date(s) not included in velocity estimation, could be list of string or text file, i.e.:\n'+\
                             '--exclude 20040502 20060708 20090103\n'+\
                             '--exclude exclude_date.txt\n'+DROP_DATE_TXT)
    parser.add_argument('--template', dest='template_file',\
                        help='template file with the following items:'+TEMPLATE)
    parser.add_argument('-o','--output', dest='outfile', help='output file name')

    inps = parser.parse_args()
    if not inps.ex_date:
        inps.ex_date = []
    return inps


############################################################################
def main(argv):
    inps = cmdLineParse()

    #print '\n********** Inversion: Time Series to Velocity ***********'
    atr = readfile.read_attribute(inps.timeseries_file)
    k = atr['FILE_TYPE']
    print 'input '+k+' file: '+inps.timeseries_file
    if not k == 'timeseries':
        sys.exit('ERROR: input file is not timeseries!') 
    h5file = h5py.File(inps.timeseries_file)

    #####################################
    ## Date Info
    dateListAll = sorted(h5file[k].keys())
    print '--------------------------------------------'
    print 'Dates from input file: '+str(len(dateListAll))
    print dateListAll

    inps.ex_date = get_exclude_date(inps, dateListAll)

    dateList = sorted(list(set(dateListAll) - set(inps.ex_date)))
    print '--------------------------------------------'
    if len(dateList) == len(dateListAll):
        print 'using all dates to calculate the velocity'
    else:
        print 'Dates used to estimate the velocity: '+str(len(dateList))
        print dateList
    print '--------------------------------------------'

    # Date Aux Info
    dates, datevector = ptime.date_list2vector(dateList)

    #####################################
    ## Inversion
    # Design matrix
    B = np.ones([len(datevector),2])
    B[:,0] = datevector
    #B_inv = np.linalg.pinv(B)
    B_inv = np.dot(np.linalg.inv(np.dot(B.T,B)), B.T)
    B_inv = np.array(B_inv, np.float32)

    # Loading timeseries
    print "Loading time series file: "+inps.timeseries_file+' ...'
    width = int(atr['WIDTH'])
    length = int(atr['FILE_LENGTH'])
    dateNum = len(dateList)
    timeseries = np.zeros([dateNum,length*width],np.float32)
    prog_bar = ptime.progress_bar(maxValue=dateNum, prefix='loading: ')
    for i in range(dateNum):
        date = dateList[i]
        timeseries[i,:] = h5file[k].get(date)[:].flatten()
        prog_bar.update(i+1, suffix=date)
    prog_bar.close()
    h5file.close()

    # Velocity Inversion
    print 'Calculating velocity ...'
    X = np.dot(B_inv, timeseries)
    velocity = np.reshape(X[0,:], [length,width])

    print 'Calculating rmse ...'
    timeseries_linear = np.dot(B, X)
    timeseries_residual = timeseries - timeseries_linear
    rmse = np.reshape(np.sqrt((np.sum((timeseries_residual)**2,0))/dateNum), [length,width])
    
    print 'Calculating the standard deviation of the estimated velocity ...'
    s1 = np.sqrt(np.sum(timeseries_residual**2,0) / (dateNum-2))
    s2 = np.sqrt(np.sum((datevector-np.mean(datevector))**2))
    std = np.reshape(s1/s2, [length,width])

    # SSt=np.sum((timeseries-np.mean(timeseries,0))**2,0)
    # SSres=np.sum(residual**2,0)
    # SS_REG=SSt-SSres
    # Rsquared=np.reshape(SS_REG/SSt,[length,width])
    ######################################################  
    # covariance of the velocities

    #####################################
    # Output file name
    if not inps.outfile:
        inps.outfile = 'velocity.h5'

    inps.outfile_rmse = os.path.splitext(inps.outfile)[0]+'Rmse'+os.path.splitext(inps.outfile)[1]
    inps.outfile_std = os.path.splitext(inps.outfile)[0]+'Std'+os.path.splitext(inps.outfile)[1]
    inps.outfile_r2 = os.path.splitext(inps.outfile)[0]+'R2'+os.path.splitext(inps.outfile)[1]

    # Attributes
    atr['date1'] = datevector[0]
    atr['date2'] = datevector[dateNum-1]

    # File Writing
    print '--------------------------------------'
    atr['FILE_TYPE'] = 'velocity'
    print 'writing >>> '+inps.outfile
    writefile.write(velocity, atr, inps.outfile)
    
    #atr['FILE_TYPE'] = 'rmse'
    print 'writing >>> '+inps.outfile_rmse
    writefile.write(rmse, atr, inps.outfile_rmse)
    
    #atr['FILE_TYPE'] = 'rmse'
    print 'writing >>> '+inps.outfile_std
    writefile.write(std, atr, inps.outfile_std)

    print 'Done.\n'
    return inps.outfile


############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])

