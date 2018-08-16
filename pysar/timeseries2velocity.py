#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013-2018, Heresh Fattahi, Zhang Yunjun     #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import os
import argparse
import numpy as np
from pysar.objects import timeseries, giantTimeseries, HDFEOS
from pysar.utils import readfile, writefile, ptime, utils as ut

dataType = np.float32


############################################################################
EXAMPLE = """example:
  timeseries2velocity.py  timeSeries_ECMWF_demErr.h5
  timeseries2velocity.py  timeseries_ECMWF_demErr_ramp.h5  --template KyushuT73F2980_2990AlosD.template
  timeseries2velocity.py  timeseries.h5  --start-date 20080201
  timeseries2velocity.py  timeseries.h5  --start-date 20080201  --end-date 20100508
  timeseries2velocity.py  timeseries.h5  --exclude-date exclude_date.txt

  timeseries2velocity.py  LS-PARAMS.h5
  timeseries2velocity.py  NSBAS-PARAMS.h5
  timeseries2velocity.py  TS-PARAMS.h5
"""

TEMPLATE = """
## estimate linear velocity from timeseries, and from tropospheric delay file if exists.
pysar.velocity.excludeDate = auto   #[exclude_date.txt / 20080520,20090817 / no], auto for exclude_date.txt
pysar.velocity.startDate   = auto   #[20070101 / no], auto for no
pysar.velocity.endDate     = auto   #[20101230 / no], auto for no
"""

DROP_DATE_TXT = """exclude_date.txt:
20040502
20060708
20090103
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Inverse velocity from time series.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=TEMPLATE+'\n'+EXAMPLE)

    parser.add_argument('timeseries_file',
                        help='Time series file for velocity inversion.')
    parser.add_argument('--start-date', dest='startDate',
                        help='start date for velocity estimation')
    parser.add_argument('--end-date', dest='endDate',
                        help='end date for velocity estimation')
    parser.add_argument('--exclude', '--ex', dest='excludeDate', nargs='+', default=[],
                        help='date(s) not included in velocity estimation, could be list of string or text file, i.e.:\n' +
                             '--exclude 20040502 20060708 20090103\n' +
                             '--exclude exclude_date.txt\n'+DROP_DATE_TXT)
    parser.add_argument('--template', '-t', dest='template_file',
                        help='template file with the following items:'+TEMPLATE)
    parser.add_argument('-o', '--output', dest='outfile',
                        help='output file name')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    inps.key = readfile.read_attribute(inps.timeseries_file)['FILE_TYPE']
    if inps.key not in ['timeseries', 'giantTimeseries', 'HDFEOS']:
        raise Exception('input file is {}, NOT timeseries!'.format(inps.key))
    return inps


def read_template2inps(template_file, inps=None):
    """Read input template file into inps.excludeDate"""
    if not inps:
        inps = cmd_line_parse()
    inpsDict = vars(inps)
    print('read options from template file: '+os.path.basename(template_file))
    template = readfile.read_template(inps.template_file)
    template = ut.check_template_auto_value(template)

    # Read template option
    prefix = 'pysar.velocity.'
    keyList = [i for i in list(inpsDict.keys()) if prefix+i in template.keys()]
    for key in keyList:
        value = template[prefix+key]
        if value:
            if key in ['startDate', 'endDate']:
                inpsDict[key] = ptime.yyyymmdd(value)
            elif key in ['excludeDate']:
                inpsDict[key] = ptime.yyyymmdd(value.replace(',', ' ').split())
    return inps


############################################################################
def read_exclude_date(inps, dateListAll):
    # Merge ex_date/startDate/endDate into ex_date
    yy_list_all = ptime.yyyymmdd2years(dateListAll)
    exDateList = []
    # 1. template_file
    if inps.template_file:
        print('read option from template file: '+inps.template_file)
        inps = read_template2inps(inps.template_file, inps)

    # 2. ex_date
    input_ex_date = list(inps.excludeDate)
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

    # 3. startDate
    if inps.startDate:
        print('start date: '+inps.startDate)
        yy_min = ptime.yyyymmdd2years(ptime.yyyymmdd(inps.startDate))
        for i in range(len(dateListAll)):
            date = dateListAll[i]
            if yy_list_all[i] < yy_min and date not in exDateList:
                print('  remove date: '+date)
                exDateList.append(date)

    # 4. endDate
    if inps.endDate:
        print('end date: '+inps.endDate)
        yy_max = ptime.yyyymmdd2years(ptime.yyyymmdd(inps.endDate))
        for i in range(len(dateListAll)):
            date = dateListAll[i]
            if yy_list_all[i] > yy_max and date not in exDateList:
                print('  remove date: '+date)
                exDateList.append(date)
    exDateList = list(set(exDateList))
    return exDateList


def read_date_info(inps):
    """Get inps.excludeDate full list
    Inputs:
        inps          - Namespace, 
    Output:
        inps.excludeDate  - list of string for exclude date in YYYYMMDD format
    """
    if inps.key == 'timeseries':
        tsobj = timeseries(inps.timeseries_file)
    elif inps.key == 'giantTimeseries':
        tsobj = giantTimeseries(inps.timeseries_file)
    elif inps.key == 'HDFEOS':
        tsobj = HDFEOS(inps.timeseries_file)
    tsobj.open()
    inps.excludeDate = read_exclude_date(inps, tsobj.dateList)

    # Date used for estimation inps.dateList
    inps.dateList = [i for i in tsobj.dateList if i not in inps.excludeDate]
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
    inps.dropDate = np.array([i not in inps.excludeDate for i in tsobj.dateList], dtype=np.bool_)
    inps.years = np.array(tsobj.yearList)[inps.dropDate]

    # output file name
    if not inps.outfile:
        outname = 'velocity'
        #if inps.excludeDate:
        #    outname += 'Ex'
        if inps.key == 'giantTimeseries':
            prefix = os.path.basename(inps.timeseries_file).split('PARAMS')[0]
            outname = prefix + outname
        outname += '.h5'
        inps.outfile = outname
    return inps


def design_matrix(years):
    """design matrix/function model of linear velocity estimation
    Parameters: years : 1D array of float in size of (numDate,),
                    date in years, e.g. 2015.1830937713896
    Returns:    A : 2D array of int in size of (numDate, 2)
    """
    A = np.ones([len(years), 2], dtype=dataType)
    A[:, 0] = years
    #A_inv = np.dot(np.linalg.inv(np.dot(A.T,A)), A.T)
    #A_inv = np.array(A_inv, dataType)
    return A


def estimate_linear_velocity(inps):
    A = design_matrix(inps.years)
    A_inv = np.array(np.linalg.pinv(A), dataType)
    # A_inv = np.array(np.linalg.inv(A.T.dot(A)).dot(A.T), dataType)  #Give wrong and different result, find reason.

    tsData, atr = readfile.read(inps.timeseries_file)
    tsData = tsData[inps.dropDate, :, :].reshape(inps.numDate, -1)
    if atr['UNIT'] == 'mm':
        tsData *= 1./1000.
    dsShape = (int(atr['LENGTH']), int(atr['WIDTH']))

    X = np.dot(A_inv, tsData)
    V = np.reshape(X[0, :], dsShape)

    tsResidual = tsData - np.dot(A, X)
    timeStd = np.sqrt(np.sum((inps.years - np.mean(inps.years))**2))
    Vstd = np.sqrt(np.sum(tsResidual**2, axis=0) / (inps.numDate-2)) / timeStd
    Vstd = Vstd.reshape(dsShape)

    # Write h5 file
    # outfileStd = '{}Std{}'.format(os.path.splitext(inps.outfile)[0],
    #                               os.path.splitext(inps.outfile)[1])

    atr['FILE_TYPE'] = 'velocity'
    atr['UNIT'] = 'm/year'
    atr['START_DATE'] = inps.dateList[0]
    atr['END_DATE'] = inps.dateList[-1]
    atr['DATE12'] = '{}_{}'.format(inps.dateList[0], inps.dateList[-1])

    dsDict = dict()
    dsDict['velocity'] = V
    dsDict['velocityStd'] = Vstd
    writefile.write(dsDict, out_file=inps.outfile, metadata=atr)
    return inps.outfile


############################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    inps = read_date_info(inps)
    inps.outfile = estimate_linear_velocity(inps)
    print('Done.')
    return inps.outfile


############################################################################
if __name__ == '__main__':
    main()
