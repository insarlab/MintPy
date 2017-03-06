#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
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

import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar._datetime as ptime
import pysar._pysar_utilities as ut


############################################################################
def yyyymmdd2years(date):
    d = datetime.datetime(*time.strptime(date,"%Y%m%d")[0:5])   
    yy = np.float(d.year) + np.float(d.month-1)/12 + np.float(d.day-1)/365
    return yy


def update_inps_from_template(inps, template_file):
    '''Update inps.ex_date with input template file'''    
    tmpl = readfile.read_template(template_file)
    if 'pysar.network.dropDate' in tmpl.keys():
        dropDateList = tmpl['pysar.network.dropDate'].replace(',',' ').split()
        if dropDateList:
            inps.ex_date += list(set(dropDate) - set(inps.ex_date))
    return inps


############################################################################
EXAMPLE='''example:
  timeseries2velocity.py  timeSeries_ECMWF_demCor.h5
  timeseries2velocity.py  timeseries_ECMWF_demCor_plane.h5 -t KyushuT73F2980_2990AlosD.template
  timeseries2velocity.py  timeseries.h5 -m 20080201
  timeseries2velocity.py  timeseries.h5 -m 20080201 -M 20100508
  timeseries2velocity.py  timeseries.h5 -E 20040502,20060708,20090103
  timeseries2velocity.py  timeseries.h5 -E drop_date.txt
'''

TEMPLATE='''
pysar.network.dropDate = 20040502 20060708 20090103
pysar.network.dropDate = drop_date.txt
'''

DROP_DATE_TXT='''drop_date.txt:
20040502
20060708
20090103
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Inverse velocity from time series.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)
    
    parser.add_argument('timeseries_file', help='Time series file for velocity inversion.')
    parser.add_argument('-m', dest='min_date', help='earliest date for velocity estimation')
    parser.add_argument('-M', dest='max_date', help='latest   date for velocity estimation')
    parser.add_argument('--ex', dest='ex_date', nargs='+',\
                        help='date(s) not included in velocity estimation, could be list of string or text file, i.e.:\n'+\
                             '--ex 20040502 20060708 20090103\n'+\
                             '--ex drop_date.txt\n'+DROP_DATE_TXT)
    parser.add_argument('-t','--template', dest='template_file',\
                        help='template file with the following items:'+TEMPLATE)
    parser.add_argument('-o','--output', dest='outfile', help='output file name')
    
    inps = parser.parse_args()
    if not inps.ex_date:  inps.ex_date = []
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
    dateListAll = ptime.yyyymmdd(dateListAll)
    yyListAll = ptime.yyyymmdd2years(dateListAll)
    print '--------------------------------------------'
    print 'Dates from input file: '+str(len(dateListAll))
    print dateListAll
    
    # Extrac exclude dates from input arguments
    inps.datesNot2include = []
    # 1. template_file
    if inps.template_file:
        inps = update_inps_from_template(inps, inps.template_file)
    
    # 2. ex_date 
    if inps.ex_date:
        for ex_date in inps.ex_date:
            if os.path.isfile(ex_date):
                ex_date = ptime.read_date_list(ex_date)
            else:
                ex_date = [ptime.yyyymmdd(ex_date)]
            inps.datesNot2include += list(set(ex_date) - set(inps.datesNot2include))
        # delete dates not existed in input file
        inps.datesNot2include = list(set(inps.datesNot2include).intersection(dateListAll))
        print 'date excluded:'+str(inps.datesNot2include)
    
    # 3. min_date
    if inps.min_date:
        inps.min_date = ptime.yyyymmdd(inps.min_date)
        print 'minimum date: '+inps.min_date
        yy_min = ptime.yyyymmdd2years(inps.min_date)
        for i in range(len(dateListAll)):
            date = dateListAll[i]
            if yyListAll[i] < yy_min and date not in inps.datesNot2include:
                print '  remove date: '+date
                inps.datesNot2include.append(date)
    
    # 4. max_date
    if inps.max_date:
        inps.max_date = ptime.yyyymmdd(inps.max_date)
        print 'minimum date: '+inps.max_date
        yy_max = ptime.yyyymmdd2years(inps.max_date)
        for i in range(len(dateListAll)):
            date = dateListAll[i]
            if yyListAll[i] > yy_max and date not in inps.datesNot2include:
                print '  remove date: '+date
                inps.datesNot2include.append(date)
    
    # Summary
    dateList = sorted(list(set(dateListAll) - set(inps.datesNot2include)))
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
    B=np.ones([len(datevector),2])
    B[:,0]=datevector
    #B1 = np.linalg.pinv(B)
    B1 = np.dot(np.linalg.inv(np.dot(B.T,B)),B.T)
    B1 = np.array(B1,np.float32)
    
    # Loading timeseries
    print "Loading time series file: "+inps.timeseries_file+' ...'
    width = int(atr['WIDTH'])
    length = int(atr['FILE_LENGTH'])
    dateNum = len(dateList)
    timeseries = np.zeros([dateNum,length*width],np.float32)
    for i in range(dateNum):
        date = dateList[i]
        ut.print_progress(i+1, dateNum, prefix='loading:', suffix=date)
        timeseries[i,:] = h5file[k].get(date)[:].flatten()
    h5file.close()

    # Velocity Inversion
    print 'Calculating velocity ...'
    x = np.dot(B1,timeseries)
    velocity = np.reshape(x[0,:],[length,width])
    
    print 'Calculating rmse ...'
    timeseries_linear = np.dot(B,x)
    rmse = np.reshape(np.sqrt((np.sum((timeseries_linear-timeseries)**2,0))/dateNum),[length,width])
    
    print 'Calculating the standard deviation of the estimated velocity ...'
    residual = timeseries_linear - timeseries
    s1 = np.sqrt(np.sum(residual**2,0)/(dateNum-2))
    s2 = np.sqrt(np.sum((datevector-np.mean(datevector))**2))
    std = np.reshape(s1/s2,[length,width])
     
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
        if inps.datesNot2include:
            inps.outfile = os.path.splitext(inps.outfile)[0]+'_ex'+os.path.splitext(inps.outfile)[1]

    inps.outfile_rmse = 'rmse_'+inps.outfile
    inps.outfile_std = 'std_'+inps.outfile
    inps.outfile_r2 = 'R2_'+inps.outfile
    
    # Attributes
    atr['date1'] = datevector[0]
    atr['date2'] = datevector[dateNum-1]
    
    # File Writing
    print '--------------------------------------'
    atr['FILE_TYPE'] = 'velocity'
    print 'writing >>> '+inps.outfile
    writefile.write(velocity, atr, inps.outfile)
    
    atr['FILE_TYPE'] = 'rmse'
    print 'writing >>> '+inps.outfile_rmse
    writefile.write(rmse, atr, inps.outfile_rmse)
    
    atr['FILE_TYPE'] = 'rmse'
    print 'writing >>> '+inps.outfile_std
    writefile.write(std, atr, inps.outfile_std)
    
    print 'Done.'
    return inps.outfile


############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])

