#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2016, Yunjun Zhang                          #
# Author:  Yunjun Zhang                                    #
############################################################
# Based on scripts writen by Heresh Fattahi
# Yunjun, Aug 2016: add read_date_list()
# Yunjun, Oct 2016: update yymmdd() for string and list input
#
# Recommended Usage:
#   import pysar._datetime as ptime
#   date_list = ptime.ifgram_date_list('unwrapIfgram.h5')


import sys
import re
import time
import datetime
from datetime import datetime as dt

import h5py
import numpy as np
import matplotlib.dates as mdates


################################################################
def yyyymmdd2years(dates):
    if isinstance(dates, basestring):
        d = dt(*time.strptime(dates,"%Y%m%d")[0:5])
        day_of_year = d.timetuple().tm_yday
        yy = float(d.year)+float(day_of_year-1)/365.25
        #yy = float(d.year) + float(d.month-1)/12 + float(d.day-1)/365.25
    elif isinstance(dates, list):
        yy = []
        for date in dates:
            d = dt(*time.strptime(date,"%Y%m%d")[0:5])
            day_of_year = d.timetuple().tm_yday
            yy.append(float(d.year)+float(day_of_year-1)/365.25)
    else:
        print 'Unrecognized date format. Only string and list supported.'
        sys.exit(1)
    return yy


def yymmdd2yyyymmdd(date):
    if date[0] == '9':      date = '19'+date
    else:                   date = '20'+date
    return date


def yyyymmdd(dates):
    if isinstance(dates,basestring):
        if len(dates) == 6:  datesOut = yymmdd2yyyymmdd(dates)
        else:                datesOut = dates
    elif isinstance(dates,list):
        datesOut = []
        for date in dates:
            if len(date) == 6:   date = yymmdd2yyyymmdd(date)
            datesOut.append(date)
    else:
        #print 'Un-recognized date input!'
        return None
    return datesOut


def yymmdd(dates):
    if isinstance(dates,basestring):
        if len(dates) == 8:  datesOut = dates[2:8]
        else:                datesOut = dates
    elif isinstance(dates,list):
        datesOut = []
        for date in dates:
            if len(date) == 8:   date = date[2:8]
            datesOut.append(date)
    else:
        #print 'Un-recognized date input!'
        return None
    return datesOut


#################################################################
def ifgram_date_list(ifgramFile, fmt='YYYYMMDD'):
    '''Read Date List from Interferogram file
        for timeseries file, use h5file['timeseries'].keys() directly
    Inputs:
        ifgramFile - string, name/path of interferograms file
        fmt        - string, output date format, choices=['YYYYMMDD','YYMMDD']
    Output:
        date_list  - list of string, date included in ifgramFile in YYYYMMDD or YYMMDD format
    '''
    if not ifgramFile:
        return []

    h5 = h5py.File(ifgramFile, 'r')
    k = h5.keys()
    if 'interferograms' in k: k = 'interferograms'
    elif 'coherence' in k: k = 'coherence'
    elif 'wrapped' in k: k = 'wrapped'
    if k not in  ['interferograms','coherence','wrapped']:
        raise ValueError('Only interferograms / coherence / wrapped are supported. Input is '+str(k))

    # Get date_list in YYMMDD format
    date_list = []
    ifgram_list = sorted(h5[k].keys())
    for ifgram in  ifgram_list:
        date12 = h5[k][ifgram].attrs['DATE12'].split('-')
        date_list.append(date12[0])
        date_list.append(date12[1])
    h5.close()
    date_list = sorted(list(set(date_list)))

    if fmt == 'YYYYMMDD':
        date_list = yyyymmdd(date_list)
    return date_list


#################################################################
def read_date_list(date_list_file):
    '''Read Date List from txt file'''
    fl = open(date_list_file,'r')
    dateList = fl.read().splitlines()
    fl.close()

    dateList = yyyymmdd(dateList)
    dateList.sort()

    return dateList


################################################################
def date_index(dateList):
    dateIndex={}
    for ni in range(len(dateList)):
        dateIndex[dateList[ni]] = ni
    return dateIndex

################################################################
def date_list2tbase(dateList):
    '''Get temporal Baseline in days with respect to the 1st date
    Input: dateList - list of string, date in YYYYMMDD or YYMMDD format
    Output:
        tbase    - list of int, temporal baseline in days
        dateDict - dict with key   - string, date in YYYYMMDD format
                             value - int, temporal baseline in days
    '''
    dateList = yyyymmdd(dateList)
    tbase=[]
    d1 = dt(*time.strptime(dateList[0],"%Y%m%d")[0:5])
    for ni in range(len(dateList)):
        d2 = dt(*time.strptime(dateList[ni],"%Y%m%d")[0:5])
        diff = d2-d1
        tbase.append(diff.days)

    ## Dictionary: key - date, value - temporal baseline
    dateDict = {}
    for i in range(len(dateList)):
        dateDict[dateList[i]] = tbase[i]
  
    return tbase, dateDict


################################################################
def date_list2vector(dateList):
    '''Get time in datetime format: datetime.datetime(2006, 5, 26, 0, 0)
    Input: dateList - list of string, date in YYYYMMDD or YYMMDD format
    Outputs:
        dates      - list of datetime.datetime objects, i.e. datetime.datetime(2010, 10, 20, 0, 0)
        datevector - list of float, years, i.e. 2010.8020547945205
    '''
    dateList = yyyymmdd(dateList)
    dates=[]
    for ni in range(len(dateList)):
        d = dt(*time.strptime(dateList[ni],"%Y%m%d")[0:5])
        dates.append(d)

    ## date in year - float format
    datevector=[]
    for i in range(len(dates)):
        dvector = dates[i].year + (dates[i].month-1)/12.0 + (dates[i].day-1)/365.0
        datevector.append(dvector)
    datevector2=[round(i,2) for i in datevector]
  
    return dates, datevector


################################################################
def auto_adjust_xaxis_date(ax, datevector, fontSize=12):
    '''Adjust X axis
    Input:
        ax : matplotlib figure axes object
        datevector : list of float, date in years
                     i.e. [2007.013698630137, 2007.521917808219, 2007.6463470319634]
    Output:
        ax  - matplotlib figure axes object
        dss - datetime.date object, xmin
        dee - datetime.date object, xmax
    '''

    # Min/Max
    ts=datevector[0] -0.2;  ys=int(ts);  ms=int((ts-ys)*12.0)
    te=datevector[-1]+0.3;  ye=int(te);  me=int((te-ye)*12.0)
    if ms>12:   ys = ys+1;   ms=1
    if me>12:   ye = ye+1;   me=1
    if ms<1:    ys = ys-1;   ms=12
    if me<1:    ye = ye-1;   me=12
    dss=datetime.date(ys,ms,1)
    dee=datetime.date(ye,me,1)
    ax.set_xlim(dss,dee)

    # Label/Tick format
    ax.fmt_xdata = mdates.DateFormatter('%Y-%m-%d %H:%M:%S')
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())

    # Label font size
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontSize)
    #fig2.autofmt_xdate()     #adjust x overlap by rorating, may enble again
    return ax, dss, dee


def list_ifgram2date12(ifgram_list):
    '''Convert ifgram list into date12 list
    Input:
        ifgram_list  - list of string in *YYMMDD-YYMMDD* or *YYMMDD_YYMMDD* format
    Output:
        date12_list  - list of string in YYMMDD-YYMMDD format
    Example:
        h5 = h5py.File('unwrapIfgram.h5','r')
        ifgram_list = sorted(h5['interferograms'].keys())
        date12_list = ptime.list_ifgram2date12(ifgram_list)
    '''
    date12_list = [str(re.findall('\d{6}[-_]\d{6}', i)[0]).replace('_','-') for i in ifgram_list]
    return date12_list


###########################Simple progress bar######################
class progress_bar:
    '''Creates a text-based progress bar. Call the object with 
    the simple `print'command to see the progress bar, which looks 
    something like this:
    [=======> 22% ]
    You may specify the progress bar's width, min and max values on init.
    
    note:
        modified from PyAPS release 1.0 (http://earthdef.caltech.edu/projects/pyaps/wiki/Main)
        Code originally from http://code.activestate.com/recipes/168639/
    
    example:
    import pysar._datetime as ptime
    date12_list = ptime.list_ifgram2date12(ifgram_list)
    prog_bar = ptime.progress_bar(maxValue=1000, prefix='calculating:')
    for i in range(1000):
        prog_bar.update(i+1, suffix=date)
        prog_bar.update(i+1, suffix=date12_list[i])
    prog_bar.close()
    '''

    def __init__(self, maxValue=100, prefix='', minValue=0, totalWidth=60):
        self.progBar = "[]" # This holds the progress bar string
        self.min = minValue
        self.max = maxValue
        self.span = maxValue - minValue
        self.width = totalWidth
        self.suffix = ''
        self.prefix = prefix
        self.reset()

    def reset(self):
        self.start_time = time.time()
        self.amount = 0 # When amount == max, we are 100% done
        self.update_amount(0) # Build progress bar string

    def update_amount(self, newAmount=0, suffix=''):
        """ Update the progress bar with the new amount (with min and max
        values set at initialization; if it is over or under, it takes the
        min or max value as a default. """
        if newAmount < self.min:
            newAmount = self.min
        if newAmount > self.max:
            newAmount = self.max
        self.amount = newAmount

        # Figure out the new percent done, round to an integer
        diffFromMin = np.float(self.amount - self.min)
        percentDone = (diffFromMin / np.float(self.span)) * 100.0
        percentDone = np.int(np.round(percentDone))

        # Figure out how many hash bars the percentage should be
        allFull = self.width - 2 - 18
        numHashes = (percentDone / 100.0) * allFull
        numHashes = np.int(np.round(numHashes))

        # Build a progress bar with an arrow of equal signs; special cases for
        # empty and full
        if numHashes == 0:
            self.progBar = '%s[>%s]' % (self.prefix, ' '*(allFull-1))
        elif numHashes == allFull:
            self.progBar = '%s[%s]' % (self.prefix, '='*allFull)
            if suffix:
                self.progBar += ' %s' % (suffix)
        else:
            self.progBar = '[%s>%s]' % ('='*(numHashes-1), ' '*(allFull-numHashes))
            # figure out where to put the percentage, roughly centered
            percentPlace = (len(self.progBar) / 2) - len(str(percentDone))
            percentString = ' ' + str(percentDone) + '% '
            # slice the percentage into the bar
            self.progBar = ''.join([self.progBar[0:percentPlace], percentString,
                    self.progBar[percentPlace+len(percentString):], ])
            # prefix and suffix
            self.progBar = self.prefix + self.progBar
            if suffix:
                self.progBar += ' %s' % (suffix)
            # time info - elapsed time and estimated remaining time
            if percentDone > 0:
                elapsed_time = time.time() - self.start_time
                self.progBar += '%5ds / %5ds' % (int(elapsed_time),
                        int(elapsed_time*(100./percentDone-1)))

    def update(self, value, every=1, suffix=''):
        """ Updates the amount, and writes to stdout. Prints a
         carriage return first, so it will overwrite the current
          line in stdout."""
        if value % every == 0 or value >= self.max:
            self.update_amount(newAmount=value, suffix=suffix)
            sys.stdout.write('\r' + self.progBar)
            sys.stdout.flush()

    def close(self):
        """Prints a blank space at the end to ensure proper printing
        of future statements."""
        print ' '
################################End of progress bar class####################################

