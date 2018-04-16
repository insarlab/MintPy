#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Aug 2015: Add -m/M/d option
# Yunjun, Jun 2016: Add -t option
# Yunjun, Aug 2015: Support drop_date txt file input


import os, sys
import time, datetime
import argparse
import h5py
import numpy as np
from pysar.utils import readfile, writefile, datetime as ptime, utils as ut
from pysar.objects import timeseries

dataType = np.float32


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

def createParser():
    parser = argparse.ArgumentParser(description='Inverse velocity from time series.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=TEMPLATE+'\n'+EXAMPLE)

    parser.add_argument('timeseries_file', help='Time series file for velocity inversion.')
    parser.add_argument('--start-date', dest='min_date', help='start date for velocity estimation')
    parser.add_argument('--end-date', dest='max_date', help='end date for velocity estimation')
    parser.add_argument('--exclude','--ex', dest='ex_date', nargs='+', default=[],\
                        help='date(s) not included in velocity estimation, could be list of string or text file, i.e.:\n'+\
                             '--exclude 20040502 20060708 20090103\n'+\
                             '--exclude exclude_date.txt\n'+DROP_DATE_TXT)
    parser.add_argument('--template','-t', dest='template_file',\
                        help='template file with the following items:'+TEMPLATE)
    parser.add_argument('-o','--output', dest='outfile', help='output file name')
    return parser


def cmdLineParse(iargs=None):
    '''Command line parser.'''
    parser = createParser()
    inps = parser.parse_args(args=iargs)
    if 'timeseries' != readfile.read_attribute(inps.timeseries_file)['FILE_TYPE']:
        print('ERROR: input file is not timeseries!')
        sys.exit(2)
    return inps


############################################################################
def read_exclude_date(inps, dateListAll):
    ##### Merge ex_date/min_date/max_date into ex_date
    yy_list_all = ptime.yyyymmdd2years(dateListAll)
    exDateList = []
    # 1. template_file
    if inps.template_file:
        print('read option from template file: '+inps.template_file)
        inps = read_template2inps(inps.template_file, inps)

    # 2. ex_date
    input_ex_date = list(inps.ex_date)
    if input_ex_date:
        for ex_date in input_ex_date:
            if os.path.isfile(ex_date):
                ex_date = ptime.read_date_list(ex_date)
            else:
                ex_date = [ptime.yyyymmdd(ex_date)]
            exDateList += list(set(ex_date) - set(exDateList))
        # delete dates not existed in input file
        exDateList = list(set(exDateList).intersection(dateListAll))
        print('exclude date:'+str(exDateList))

    # 3. min_date
    if inps.min_date:
        print('start date: '+inps.min_date)
        yy_min = ptime.yyyymmdd2years(ptime.yyyymmdd(inps.min_date))
        for i in range(len(dateListAll)):
            date = dateListAll[i]
            if yy_list_all[i] < yy_min and date not in exDateList:
                print('  remove date: '+date)
                exDateList.append(date)

    # 4. max_date
    if inps.max_date:
        print('end date: '+inps.max_date)
        yy_max = ptime.yyyymmdd2years(ptime.yyyymmdd(inps.max_date))
        for i in range(len(dateListAll)):
            date = dateListAll[i]
            if yy_list_all[i] > yy_max and date not in exDateList:
                print('  remove date: '+date)
                exDateList.append(date)
    exDateList = list(set(exDateList))
    return exDateList


def read_date_info(inps):
    '''Get inps.ex_date full list
    Inputs:
        inps          - Namespace, 
    Output:
        inps.ex_date  - list of string for exclude date in YYYYMMDD format
    '''
    tsobj = timeseries(inps.timeseries_file)
    tsobj.open()
    inps.ex_date = read_exclude_date(inps, tsobj.dateList)

    ##### Date used for estimation inps.dateList
    inps.dateList = [i for i in tsobj.dateList if i not in inps.ex_date]
    inps.numDate = len(inps.dateList)
    print('-'*50)
    print('dates from input file: {}\n{}'.format(tsobj.numDate, tsobj.dateList))
    print('-'*50)
    if len(inps.dateList) == len(tsobj.dateList):
        print('using all dates to calculate the velocity')
    else:
        print('dates used to estimate the velocity: {}\n{}'.format(inps.numDate, inps.dateList))
    print('-'*50)
    # Date Aux Info
    inps.dropDate = np.array([i not in inps.ex_date for i in tsobj.dateList], dtype=np.bool_)
    inps.years = np.array(tsobj.yearList)[inps.dropDate]
    tsobj.close()
    return inps


def read_template2inps(template_file, inps=None):
    '''Read input template file into inps.ex_date'''
    if not inps:
        inps = cmdLineParse()
    template = readfile.read_template(template_file)

    # Read template option
    prefix = 'pysar.velocity.'
    key = prefix+'excludeDate'
    if key in template.keys():
        value = template[key]
        if value == 'auto':
            inps.ex_date = ['exclude_date.txt']
        elif value == 'no':
            inps.ex_date = []
        else:
            inps.ex_date = value.replace(',',' ').split()

    key = prefix+'startDate'
    if key in template.keys():
        value = template[key]
        if value not in ['auto','no']:
            inps.min_date = ptime.yyyymmdd(value)

    key = prefix+'endDate'
    if key in template.keys():
        value = template[key]
        if value not in ['auto','no']:
            inps.max_date = ptime.yyyymmdd(value)

    return inps


def design_matrix(years):
    '''design matrix/function model of linear velocity estimation
    Parameters: years : 1D array of float in size of (numDate,),
                    date in years, e.g. 2015.1830937713896
    Returns:    A : 2D array of int in size of (numDate, 2)
    '''
    A = np.ones([len(years),2], dtype=dataType)
    A[:,0] = years
    #A_inv = np.dot(np.linalg.inv(np.dot(A.T,A)), A.T)
    #A_inv = np.array(A_inv, dataType)
    return A


def estimateVelocity(inps):
    A = design_matrix(inps.years)
    A_inv = np.array(np.linalg.pinv(A), dataType)
    #A_inv = np.array(np.linalg.inv(A.T.dot(A)).dot(A.T), dataType)  #Give wrong and different result, find reason.

    tsobj = timeseries(inps.timeseries_file)
    tsData = tsobj.read()[inps.dropDate,:,:].reshape(inps.numDate, -1)
    dsShape = (tsobj.length, tsobj.width)

    X = np.dot(A_inv, tsData)
    V = np.reshape(X[0,:], dsShape)

    tsResidual = tsData - np.dot(A, X)
    timeStd = np.sqrt(np.sum((inps.years - np.mean(inps.years))**2))
    Vstd = np.sqrt(np.sum(tsResidual**2, axis=0) / (inps.numDate-2)) / timeStd
    Vstd = Vstd.reshape(dsShape)

    ## Write h5 file
    if not inps.outfile:
        if inps.ex_date:
            inps.outfile = 'velocityEx.h5'
        else:
            inps.outfile = 'velocity.h5'
    print('create HDF5 file: {} with w mode'.format(inps.outfile))
    f = h5py.File(inps.outfile, 'w')

    dsName = 'velocity'
    print('create dataset /{:<12} of {:<10} in size of {}'.format(dsName, str(dataType), dsShape))
    ds = f.create_dataset(dsName, data=V, dtype=dataType, chunks=True, compression="gzip")
    ds.attrs['Title'] = dsName
    ds.attrs['MinValue'] = np.nanmin(V)
    ds.attrs['MaxValue'] = np.nanmax(V)

    dsName = 'velocityStd'
    print('create dataset /{:<12} of {:<10} in size of {}'.format(dsName, str(dataType), dsShape))
    ds = f.create_dataset(dsName, data=Vstd, dtype=dataType, chunks=True, compression="gzip")
    ds.attrs['Title'] = dsName
    ds.attrs['MinValue'] = np.nanmin(Vstd)
    ds.attrs['MaxValue'] = np.nanmax(Vstd)

    atr = tsobj.metadata.copy()
    atr['FILE_TYPE'] = 'velocity'
    atr['UNIT'] = 'm/year'
    for key, value in atr.items():
        f.attrs[key] = value
    f.close()
    print('finished writing to {}'.format(inps.outfile))
    return inps.outfile


############################################################################
def main(iargs=None):
    inps = cmdLineParse(iargs)
    inps = read_date_info(inps)
    inps.outfile = estimateVelocity(inps)
    print('Done.')
    return inps.outfile


############################################################################
if __name__ == '__main__':
    main()

