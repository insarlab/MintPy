############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2016               #
############################################################
# Recommend import:
#   from mintpy.utils import ptime

import os
import sys
import re
import time
import datetime as dt
import numpy as np
from mintpy.objects.progress import progressBar


################################################################
def get_date_str_format(date_str):
    """Get the datetime string format as defined in:
    https://docs.python.org/3.7/library/datetime.html#strftime-and-strptime-behavior

    Parameters: date_str - str, date in one of the following formats:
                            YYYYMMDDTHHMM
                            YYYYMMDD
                            YYMMDD
    Returns:    date_str_format - str, datetime string format
    """
    if isinstance(date_str, list):
        date_str = date_str[0]

    try:
        date_str = date_str.decode('utf8')
    except:
        pass

    date_str_format = None
    if len(re.findall('\d{8}T\d{6}', date_str)) > 0:
        date_str_format = '%Y%m%dT%H%M%S'

    elif len(re.findall('\d{8}T\d{4}', date_str)) > 0:
        date_str_format = '%Y%m%dT%H%M'

    elif len(re.findall('\d{6}T\d{4}', date_str)) > 0:
        date_str_format = '%y%m%dT%H%M'

    elif len(re.findall('\d{8}', date_str)) > 0:
        date_str_format = '%Y%m%d'

    elif len(re.findall('\d{6}', date_str)) > 0:
        date_str_format = '%y%m%d'

    else:
        raise ValueError('un-recognized date string format for "{}"!'.format(date_str))

    return date_str_format


def round_seconds(datetime_obj):
    """Round datetime object to the nearest second.
    Link: https://stackoverflow.com/questions/47792242/rounding-time-off-to-the-nearest-second-python
    """
    datetime_obj_out = datetime_obj
    if datetime_obj_out.microsecond >= 5e5:
        datetime_obj_out += dt.timedelta(seconds=1)
    return datetime_obj_out.replace(microsecond=0)


def yyyymmdd2season(date_str):
    """Determine the season of input date in YYYYMMDD format

    Parameters: date_str - str, date in YYYYMMDD format
    Returns:    season   - str, season in ['WINTER', 'SPRING', 'SUMMER', 'FALL']
    """
    # get day of the year
    date_str = yyyymmdd(date_str)
    yday = dt.datetime.strptime(date_str, "%Y%m%d").timetuple().tm_yday

    # determine the season
    season = None
    if yday < 60 or yday > 330:
        season = 'WINTER'
    elif yday < 152:
        season = 'SPRING'
    elif yday < 244:
        season = 'SUMMER'
    else:
        season = 'FALL'
    return season


def datenum2datetime(datenum):
    """Convert Matlab datenum into Python datetime.
    Parameters: datenum : Date in datenum format, i.e. 731763.5
    Returns:    datetime: Date in datetime.datetime format, datetime.datetime(2003, 7, 1, 12, 0)
    """
    return dt.datetime.fromordinal(int(datenum)) \
           + dt.timedelta(days=datenum % 1) \
           - dt.timedelta(days=366)


def decimal_year2datetime(years):
    """read date in 2002.40657084 to datetime format
    Parameters: years    : (list of) float or str for years
    Returns:    years_dt : (list of) datetime.datetime objects
    """
    def decimal_year2datetime1(x):
        x = float(x)
        year = np.floor(x).astype(int)
        yday = np.floor((x - year) * 365.25).astype(int) + 1
        x2 = '{:d}-{:d}'.format(year, yday)
        try:
            xt = dt.datetime.strptime(x2, "%Y-%j")
        except:
            raise ValueError('wrong format: ',x)
        return xt

    if isinstance(years, (float, str)):
        years_dt = decimal_year2datetime1(years)

    elif isinstance(years, list):
        years_dt = []
        for year in years:
            years_dt.append(decimal_year2datetime1(year))

    else:
        raise ValueError('unrecognized input format: {}. Only float/str/list are supported.'.format(type(years)))
    return years_dt


def yyyymmdd2years(dates):
    """Convert date(s) string into float number in the unit of year

    Parameters: dates - (list of) str, date in YYYYMMDD format
    Returns:    years - (list of) float, years including the date and time info
    """

    # make a copy in list of input arg
    if isinstance(dates, str):
        date_list = [dates]
    else:
        date_list = list(dates)

    date_format = get_date_str_format(date_list[0])

    years = []
    for date_str in date_list:
        d = dt.datetime.strptime(date_str, date_format)
        y = (d.year + (d.timetuple().tm_yday - 1) / 365.25 +
             d.hour / (365.25 * 24) +
             d.minute / (365.25 * 24 * 60) +
             d.second / (365.25 * 24 * 60 * 60))
        years.append(y)

    if isinstance(dates, str):
        years = years[0]

    return years


def yymmdd2yyyymmdd(date):
    """Convert date str from YYMMDD to YYYYMMDD format"""
    if date[0] == '9':
        date = '19'+date
    else:
        date = '20'+date
    return date


def yy2yyyy(year):
    """Convert year str from YY to YYYY format"""
    if year[0] == '9':
        year = '19'+year
    else:
        year = '20'+year
    return datyeare


def yyyymmdd(dates):
    """Convert date str from (YY)YYMMDD(THHMM) to YYYYMMDD(THHMM) format"""
    if isinstance(dates, str):
        if len(dates.split('T')[0]) == 6:
            datesOut = yymmdd2yyyymmdd(dates)
        else:
            datesOut = dates

    elif isinstance(dates, list):
        datesOut = []
        for date in dates:
            if len(date.split('T')[0]) == 6:
                date = yymmdd2yyyymmdd(date)
            datesOut.append(date)

    else:
        return None
    return datesOut


def yymmdd(dates):
    """Convert date str in (YY)YYMMDD to YYMMDD format"""
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
def read_date_txt(date_file):
    """Read Date List from txt file"""
    #date_list = np.loadtxt(date_file, dtype=bytes).astype(str)
    date_list = []

    if os.path.isfile(date_file):
        # read text file
        with open(date_file, 'r') as f:
            date_list = f.read().splitlines()

        # format
        date_list = sorted(yyyymmdd(date_list))
    return date_list


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
        if d.endswith(('.txt','.cfg','.dat')):
            if os.path.isfile(d):
                ds = read_date_txt(d)
                date_list_out += ds

        else:
            ds = [d]
            date_list_out += ds
    date_list_out = sorted(yyyymmdd(list(set(date_list_out))))

    # exclude date not in date_list_ref
    if date_list_all:
        date_list_out = sorted(list(set(date_list_out).intersection(date_list_all)))

    return date_list_out


################################################################
def date_list2tbase(date_list):
    """Get temporal Baseline in days with respect to the 1st date
    Parameters: date_list - list of string, date in YYYYMMDD or YYMMDD format
    Returns:    tbase     - list of int, temporal baseline in days
                dateDict  - dict with key   - string, date in YYYYMMDD format
                                      value - int, temporal baseline in days
    """
    # date str to dt object
    date_list = yyyymmdd(date_list)
    date_format = get_date_str_format(str(date_list))
    dates = [dt.datetime.strptime(i, date_format) for i in date_list]

    # dt object to time difference in days
    tbase = []
    for date in dates:
        date_delta = date - dates[0]
        tbase_i = date_delta.days + date_delta.seconds / (24 * 60 * 60)
        tbase.append(tbase_i)

    # Dictionary: key - date, value - temporal baseline
    dateDict = {}
    for i in range(len(date_list)):
        dateDict[date_list[i]] = tbase[i]
    return tbase, dateDict


################################################################
def date_list2vector(date_list):
    """Get time in datetime format: datetime.datetime(2006, 5, 26, 0, 0)
    Parameters: date_list  - list of string, date in YYYYMMDD or YYMMDD format
    Returns:    dates      - list of datetime.datetime objects, i.e. datetime.datetime(2010, 10, 20, 0, 0)
                datevector - list of float, years, i.e. 2010.8020547945205
    """
    date_list = yyyymmdd(date_list)
    date_format = get_date_str_format(str(date_list))
    dates = [dt.datetime.strptime(i, date_format) for i in date_list]

    # date in year - float format
    datevector = []
    for d in dates:
        date_vec = (d.year + (d.timetuple().tm_yday - 1) / 365.25 +
                    d.hour / (365.25 * 24) +
                    d.minute / (365.25 * 24 * 60) +
                    d.second / (365.25 * 24 * 60 * 60))
        datevector.append(date_vec)

    return dates, datevector


################################################################
def get_date_range(dmin, dmax):
    """Make a list of dates with one-day interval [dmin, dmax]
    Parameters: dmin      : str in YYYYMMDD format
                dmax      : str in YYYYMMDD format
    Returns:    date_list : list of str in YYYYMMDD format
    """
    t1 = '{}-{}-{}'.format(dmin[:4], dmin[4:6], dmin[6:])
    t2 = '{}-{}-{}'.format(dmax[:4], dmax[4:6], dmax[6:])
    date_objs = np.arange(t1, t2, dtype='datetime64[D]').astype('M8[D]').astype('O')
    date_list = [obj.strftime('%Y%m%d') for obj in date_objs] + [dmax]
    return date_list


def utc2solar_time(utc_time, longitude):
    """Convert UTC time to solar local time.
    Solar time: https://en.wikipedia.org/wiki/Solar_time
    Link: https://stackoverflow.com/questions/13314626

    Parameters: utc_time   - datetime.datetime object for the UTC time
                longitude  - float, longitude of the observer in degrees
    Returns:    solar_time - datetime.datetime object for the local solar time
    Example:    utc_time = dt.datetime(2015, 2, 9, 3, 18, 48)
                solar_time = ptime.utc2solar_time(utc_time, 130.7)
    """
    from math import pi, cos, sin
    # use 366 for leap years
    if utc_time.year % 4 == 0 and utc_time.year % 100 != 0 and utc_time.year % 400 != 0:
        year_len = 366
    else:
        year_len = 365

    gamma = 2 * pi / year_len * (utc_time.timetuple().tm_yday - 1 + float(utc_time.hour - 12) / 24)
    eqtime = 229.18 * (0.000075 + 0.001868 * cos(gamma) - 0.032077 * sin(gamma) \
             - 0.014615 * cos(2 * gamma) - 0.040849 * sin(2 * gamma))
    #decl = 0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) \
    #       - 0.006758 * cos(2 * gamma) + 0.000907 * sin(2 * gamma) \
    #       - 0.002697 * cos(3 * gamma) + 0.00148 * sin(3 * gamma)
    time_offset = eqtime + 4 * longitude
    tst = utc_time.hour * 60 + utc_time.minute + utc_time.second / 60 + time_offset
    solar_time = dt.datetime.combine(utc_time.date(), dt.time(0)) + dt.timedelta(minutes=tst)

    return solar_time

