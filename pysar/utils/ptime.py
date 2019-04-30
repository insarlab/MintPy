############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2016-2018, Zhang Yunjun, Heresh Fattahi     #
# Author:  Zhang Yunjun, Heresh Fattahi, 2016              #
############################################################
# Recommend import:
#   from pysar.utils import ptime

import os
import sys
import time
from datetime import datetime as dt, timedelta

import h5py
import numpy as np


################################################################
def datenum2datetime(datenum):
    """Convert Matlab datenum into Python datetime.
    Parameters: datenum : Date in datenum format, i.e. 731763.5
    Returns:    datetime: Date in datetime.datetime format, datetime.datetime(2003, 7, 1, 12, 0)
    """
    return dt.fromordinal(int(datenum)) \
           + timedelta(days=datenum % 1) \
           - timedelta(days=366)


def decimal_year2datetime(x):
    """read date in 2002.40657084 to datetime format"""
    x = float(x)
    year = np.floor(x).astype(int)
    yday = np.floor((x - year) * 365.25).astype(int) + 1
    x2 = '{:d}-{:d}'.format(year, yday)
    try:
        xt = dt(*time.strptime(x2, "%Y-%j")[0:5])
    except:
        raise ValueError('wrong format: ',x)
    return xt


def yyyymmdd2years(dates):
    if isinstance(dates, str):
        d = dt(*time.strptime(dates, "%Y%m%d")[0:5])
        yy = float(d.year)+float(d.timetuple().tm_yday-1)/365.25
    elif isinstance(dates, list):
        yy = []
        for date in dates:
            d = dt(*time.strptime(date, "%Y%m%d")[0:5])
            yy.append(float(d.year)+float(d.timetuple().tm_yday-1)/365.25)
    else:
        raise ValueError('Unrecognized date format. Only string and list supported.')
    return yy


def yymmdd2yyyymmdd(date):
    if date[0] == '9':
        date = '19'+date
    else:
        date = '20'+date
    return date


def yyyymmdd(dates):
    if isinstance(dates, str):
        if len(dates) == 6:
            datesOut = yymmdd2yyyymmdd(dates)
        else:
            datesOut = dates
    elif isinstance(dates, list):
        datesOut = []
        for date in dates:
            if len(date) == 6:
                date = yymmdd2yyyymmdd(date)
            datesOut.append(date)
    else:
        # print 'Un-recognized date input!'
        return None
    return datesOut


def yymmdd(dates):
    if isinstance(dates, str):
        if len(dates) == 8:
            datesOut = dates[2:8]
        else:
            datesOut = dates
    elif isinstance(dates, list):
        datesOut = []
        for date in dates:
            if len(date) == 8:
                date = date[2:8]
            datesOut.append(date)
    else:
        # print 'Un-recognized date input!'
        return None
    return datesOut


def yyyymmdd_date12(date12_list):
    """Convert date12 into YYYYMMDD_YYYYMMDD format"""
    m_dates = yyyymmdd([i.replace('-', '_').split('_')[0] for i in date12_list])
    s_dates = yyyymmdd([i.replace('-', '_').split('_')[1] for i in date12_list])
    date12_list = ['{}_{}'.format(m, s) for m, s in zip(m_dates, s_dates)]
    return date12_list

def yymmdd_date12(date12_list):
    """Convert date12 into YYMMDD-YYMMDD format"""
    m_dates = yymmdd([i.replace('-', '_').split('_')[0] for i in date12_list])
    s_dates = yymmdd([i.replace('-', '_').split('_')[1] for i in date12_list])
    date12_list = ['{}-{}'.format(m, s) for m, s in zip(m_dates, s_dates)]
    return date12_list


#################################################################
def ifgram_date_list(ifgramFile, fmt='YYYYMMDD'):
    """Read Date List from Interferogram file
        for timeseries file, use h5file['timeseries'].keys() directly
    Inputs:
        ifgramFile - string, name/path of interferograms file
        fmt        - string, output date format, choices=['YYYYMMDD','YYMMDD']
    Output:
        date_list  - list of string, date included in ifgramFile in YYYYMMDD or YYMMDD format
    """
    if not ifgramFile:
        return []

    h5 = h5py.File(ifgramFile, 'r')
    k = list(h5.keys())
    if 'interferograms' in k: k = 'interferograms'
    elif 'coherence' in k: k = 'coherence'
    elif 'wrapped' in k: k = 'wrapped'
    else: k = k[0]
    if k not in  ['interferograms','coherence','wrapped']:
        raise ValueError('Only interferograms / coherence / wrapped are supported. Input is '+str(k))

    # Get date_list in YYMMDD format
    date_list = []
    ifgram_list = sorted(h5[k].keys())
    for ifgram in ifgram_list:
        date12 = h5[k][ifgram].attrs['DATE12']
        try:
            date12 = date12.decode('utf-8')
        except:
            pass
        date_list.append(date12.split('_')[0])
        date_list.append(date12.split('_')[1])
    h5.close()
    date_list = sorted(list(set(date_list)))

    if fmt == 'YYYYMMDD':
        date_list = yyyymmdd(date_list)
    return date_list


#################################################################
def read_date_txt(date_file):
    """Read Date List from txt file"""
    # read text file
    #dateList = np.loadtxt(date_file, dtype=bytes).astype(str)
    with open(date_file, 'r') as f:
        dateList = f.read().splitlines()
    # format
    dateList = sorted(yyyymmdd(dateList))
    return dateList


def read_date_list(date_list_in, date_list_all=None):
    """Read Date List
    Parameters: date_list_in  : list of str / text file
                date_list_all : list of str in YYYYMMDD format
    Returns:    date_list_out : list of str in YYYYMMDD format
    """
    if not date_list_in:
        return []
    elif isinstance(date_list_in, str):
        date_list_in = [date_list_in]

    # read date_list_in
    date_list_out = []
    for d in date_list_in:
        if os.path.isfile(d):
            ds = read_date_txt(d)
        else:
            ds = [d]
        date_list_out += ds
    date_list_out = sorted(yyyymmdd(list(set(date_list_out))))

    # exclude date not in date_list_ref
    if date_list_all:
        date_list_out = list(set(date_list_out).intersection(date_list_all))

    return date_list_out


################################################################
def date_index(dateList):
    dateIndex = {}
    for ni in range(len(dateList)):
        dateIndex[dateList[ni]] = ni
    return dateIndex

################################################################
def date_list2tbase(dateList):
    """Get temporal Baseline in days with respect to the 1st date
    Input: dateList - list of string, date in YYYYMMDD or YYMMDD format
    Output:
        tbase    - list of int, temporal baseline in days
        dateDict - dict with key   - string, date in YYYYMMDD format
                             value - int, temporal baseline in days
    """
    dateList = yyyymmdd(dateList)
    dates = [dt(*time.strptime(i, "%Y%m%d")[0:5]) for i in dateList]
    tbase = [(i-dates[0]).days for i in dates]

    # Dictionary: key - date, value - temporal baseline
    dateDict = {}
    for i in range(len(dateList)):
        dateDict[dateList[i]] = tbase[i]
    return tbase, dateDict


################################################################
def date_list2vector(dateList):
    """Get time in datetime format: datetime.datetime(2006, 5, 26, 0, 0)
    Input: dateList - list of string, date in YYYYMMDD or YYMMDD format
    Outputs:
        dates      - list of datetime.datetime objects, i.e. datetime.datetime(2010, 10, 20, 0, 0)
        datevector - list of float, years, i.e. 2010.8020547945205
    """
    dateList = yyyymmdd(dateList)
    dates = [dt(*time.strptime(i, "%Y%m%d")[0:5]) for i in dateList]
    # date in year - float format
    datevector = [i.year + (i.timetuple().tm_yday - 1)/365.25 for i in dates]
    datevector2 = [round(i, 2) for i in datevector]
    return dates, datevector


################################################################
def closest_weather_product_time(sar_acquisition_time, grib_source='ECMWF'):
    """Find closest available time of weather product from SAR acquisition time
    Inputs:
        sar_acquisition_time - string, SAR data acquisition time in seconds
        grib_source - string, Grib Source of weather reanalysis product
    Output:
        grib_hr - string, time of closest available weather product 
    Example:
        '06' = closest_weather_product_time(atr['CENTER_LINE_UTC'])
        '12' = closest_weather_product_time(atr['CENTER_LINE_UTC'], 'NARR')
    """
    # Get hour/min of SAR acquisition time
    sar_time = float(sar_acquisition_time)

    # Find closest time in available weather products
    grib_hr_list = [0, 6, 12, 18]
    grib_hr = int(min(grib_hr_list, key=lambda x: abs(x-sar_time/3600.)))

    # Adjust time output format
    grib_hr = "%02d" % grib_hr
    return grib_hr


###########################Simple progress bar######################
class progressBar:
    """Creates a text-based progress bar. Call the object with 
    the simple print command to see the progress bar, which looks 
    something like this:
    [=======> 22%       ]
    You may specify the progress bar's min and max values on init.

    note:
        modified from PyAPS release 1.0 (http://earthdef.caltech.edu/projects/pyaps/wiki/Main)
        Code originally from http://code.activestate.com/recipes/168639/

    example:
        from pysar.utils import ptime
        date12_list = ptime.list_ifgram2date12(ifgram_list)
        prog_bar = ptime.progressBar(maxValue=1000, prefix='calculating:')
        for i in range(1000):
            prog_bar.update(i+1, suffix=date)
            prog_bar.update(i+1, suffix=date12_list[i])
        prog_bar.close()
    """

    def __init__(self, maxValue=100, prefix='', minValue=0, totalWidth=70, print_msg=True):
        self.prog_bar = "[]"  # This holds the progress bar string
        self.min = minValue
        self.max = maxValue
        self.span = maxValue - minValue
        self.suffix = ''
        self.prefix = prefix

        self.print_msg = print_msg
        ## calculate total width based on console width
        #rows, columns = os.popen('stty size', 'r').read().split()
        #self.width = round(int(columns) * 0.7 / 10) * 10
        self.width = totalWidth
        self.reset()

    def reset(self):
        self.start_time = time.time()
        self.amount = 0  # When amount == max, we are 100% done
        self.update_amount(0)  # Build progress bar string

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
            self.prog_bar = '%s[>%s]' % (self.prefix, ' '*(allFull-1))
        elif numHashes == allFull:
            self.prog_bar = '%s[%s]' % (self.prefix, '='*allFull)
            if suffix:
                self.prog_bar += ' %s' % (suffix)
        else:
            self.prog_bar = '[%s>%s]' % ('='*(numHashes-1), ' '*(allFull-numHashes))
            # figure out where to put the percentage, roughly centered
            percentPlace = int(len(self.prog_bar)/2 - len(str(percentDone)))
            percentString = ' ' + str(percentDone) + '% '
            # slice the percentage into the bar
            self.prog_bar = ''.join([self.prog_bar[0:percentPlace],
                                     percentString,
                                     self.prog_bar[percentPlace+len(percentString):]])
            # prefix and suffix
            self.prog_bar = self.prefix + self.prog_bar
            if suffix:
                self.prog_bar += ' %s' % (suffix)
            # time info - elapsed time and estimated remaining time
            if percentDone > 0:
                elapsed_time = time.time() - self.start_time
                self.prog_bar += '%5ds / %5ds' % (int(elapsed_time),
                                                  int(elapsed_time * (100./percentDone-1)))

    def update(self, value, every=1, suffix=''):
        """ Updates the amount, and writes to stdout. Prints a
         carriage return first, so it will overwrite the current
          line in stdout."""
        if value % every == 0 or value >= self.max:
            self.update_amount(newAmount=value, suffix=suffix)
            if self.print_msg:
                sys.stdout.write('\r' + self.prog_bar)
                sys.stdout.flush()

    def close(self):
        """Prints a blank space at the end to ensure proper printing
        of future statements."""
        if self.print_msg:
            print(' ')
################################End of progress bar class####################################
