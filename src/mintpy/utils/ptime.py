"""Utilities for date/time operations."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2016               #
############################################################
# Recommend import:
#   from mintpy.utils import ptime

import datetime as dt
import os
import re

import numpy as np

from mintpy.objects.progress import progressBar


################################################################
def get_compact_isoformat(date_str):
    """Get the "compact-looking" isoformat of the input datetime string.
    Parameters: date_str   - str, an example date string
    Returns:    iso_format - str, date string in "compact" iso format
    """
    date_str = date_str[0] if isinstance(date_str, list) else date_str
    iso_format = get_date_str_format(date_str)
    iso_format = iso_format.replace('y', 'Y')
    iso_format = iso_format.replace('%Y%m%d', '%Y-%m-%d')
    iso_format = iso_format.replace('%H%M%S', '%H:%M:%S')
    iso_format = iso_format.replace('%H%M',   '%H:%M')
    return iso_format


def get_date_str_format(date_str):
    """Get the datetime string format as defined in:
    https://docs.python.org/3.7/library/datetime.html#strftime-and-strptime-behavior

    Parameters: date_str - str, date in one of the following formats:
                            YYYYMMDDTHHMM
                            YYYYMMDD
                            YYMMDD
    Returns:    date_str_format - str, datetime string format
    """
    date_str = date_str[0] if isinstance(date_str, list) else date_str

    try:
        date_str = date_str.decode('utf8')
    except:
        pass

    date_str_format = None
    if len(re.findall(r'\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}', date_str)) > 0:
        date_str_format = '%Y-%m-%dT%H:%M:%S'

    elif len(re.findall(r'\d{4}-\d{2}-\d{2}T\d{2}:\d{2}', date_str)) > 0:
        date_str_format = '%Y-%m-%dT%H:%M'

    elif len(re.findall(r'\d{4}-\d{2}-\d{2}T\d{2}', date_str)) > 0:
        date_str_format = '%Y-%m-%dT%H'

    elif len(re.findall(r'\d{4}-\d{2}-\d{2}', date_str)) > 0:
        date_str_format = '%Y-%m-%d'

    elif len(re.findall(r'\d{8}T\d{6}', date_str)) > 0:
        date_str_format = '%Y%m%dT%H%M%S'

    elif len(re.findall(r'\d{8}T\d{4}', date_str)) > 0:
        date_str_format = '%Y%m%dT%H%M'

    elif len(re.findall(r'\d{6}T\d{4}', date_str)) > 0:
        date_str_format = '%y%m%dT%H%M'

    elif len(re.findall(r'\d{8}', date_str)) > 0:
        date_str_format = '%Y%m%d'

    elif len(re.findall(r'\d{6}', date_str)) > 0:
        date_str_format = '%y%m%d'

    else:
        raise ValueError(f'un-recognized date string format for "{date_str}"!')

    return date_str_format


def get_date12_from_path(file_path):
    """Get date12 str from a given file path.

    Parameters: file_path  - str, path to a file that contains date1/2 info
    Returns:    date12_str - str, date12 in (YY)YYMMDD(THHMM)[-_](YY)YYMMDD(THHMM) format
    """

    # support date string format
    date12_fmts = [
        r'\d{8}T\d{4}[-_]\d{8}T\d{4}',   # %Y%m%dT%H%M
        r'\d{8}[-_]\d{8}',               # %Y%m%d
        r'\d{6}[-_]\d{6}',               # %y%m%d
    ]

    # search date12 pattern part by part in the file path
    date12_str = None
    parts = file_path.split(os.sep)[::-1]
    for part in parts:
        for date12_fmt in date12_fmts:
            if len(re.findall(date12_fmt, part)) > 0:
                date12_str = re.findall(date12_fmt, part)[0]
                break

        # exit remaining parts searching
        if date12_str:
            break

    if not date12_str:
        raise ValueError(f'NO date12 str found in path: {file_path}!')

    return date12_str


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
        x2 = f'{year:d}-{yday:d}'
        try:
            xt = dt.datetime.strptime(x2, "%Y-%j")
        except:
            raise ValueError('wrong format: ',x)
        return xt

    if isinstance(years, (float, np.float32, np.float64, str)):
        years_dt = decimal_year2datetime1(years)

    elif isinstance(years, list):
        years_dt = []
        for year in years:
            years_dt.append(decimal_year2datetime1(year))

    else:
        raise ValueError(f'unrecognized input format: {type(years)}. Only float/str/list are supported.')
    return years_dt


def yyyymmdd2years(dates, seconds=0):
    """Convert date(s) string into float number in the unit of year
    Parameters: dates   - (list of) str, date in YYYYMMDD format
                seconds - float or str, time of the day info in seconds
    Returns:    years   - (list of) float, years including the date and time info
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

        # add time of the day info if:
        # 1) seconds arg is valid AND
        # 2) no time info from dates arg
        if seconds and 'T' not in date_format:
            y += float(seconds) / (365.25 * 24 * 60 * 60)

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
    return year


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


def yyyymmdd_date12(date12_list_in):
    """Convert date12 into YYYYMMDD_YYYYMMDD format
    Parameters: date12_list_in  - (list of) str
    Returns:    date12_list_out - (list of) str in YYYYMMDD_YYYYMMDD format
    """
    # endure list type input
    if isinstance(date12_list_in, str):
        date12_list = [date12_list_in]
    else:
        date12_list = list(date12_list_in)

    # convert
    m_dates = yyyymmdd([i.replace('-', '_').split('_')[0] for i in date12_list])
    s_dates = yyyymmdd([i.replace('-', '_').split('_')[1] for i in date12_list])
    date12_list_out = [f'{m}_{s}' for m, s in zip(m_dates, s_dates)]

    # ensure same type output
    if isinstance(date12_list_in, str):
        date12_list_out = date12_list_out[0]

    return date12_list_out


def yymmdd_date12(date12_list_in):
    """Convert date12 into YYMMDD-YYMMDD format
    Parameters: date12_list_in  - (list of) str
    Returns:    date12_list_out - (list of) str in YYMMDD-YYMMDD format
    """
    # ensure list type input
    if isinstance(date12_list_in, str):
        date12_list = [date12_list_in]
    else:
        date12_list = list(date12_list_in)

    # convert
    m_dates = yymmdd([i.replace('-', '_').split('_')[0] for i in date12_list])
    s_dates = yymmdd([i.replace('-', '_').split('_')[1] for i in date12_list])
    date12_list_out = [f'{m}-{s}' for m, s in zip(m_dates, s_dates)]

    # ensure same type output
    if isinstance(date12_list_in, str):
        date12_list_out = date12_list_out[0]

    return date12_list_out


#################################################################
def read_date_txt(date_file):
    """Read Date List from txt file"""
    #date_list = np.loadtxt(date_file, dtype=bytes).astype(str)
    date_list = []

    if os.path.isfile(date_file):
        # read text file
        with open(date_file) as f:
            date_list = f.read().splitlines()
        # ignore invalid lines, e.g. empty
        date_list = [x for x in date_list if x]

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


def get_exclude_date_list(date_list, start_date=None, end_date=None, exclude_date=None):
    """Get exclude date list from input options (start/end/ex_date).

    Parameters: date_list    - list of str, all dates in YYYYMMDD(THHMM) format
                start_date   - str, starting date
                end_date     - str, ending date
                exclude_date - list of str, exclude date in YYYYMMDD or text file
    Returns:    ex_date_list - list of str, exclude date
    """

    year_list = yyyymmdd2years(date_list)
    ex_date_list = []

    # exclude_date
    if exclude_date:
        ex_date_list += read_date_list(list(exclude_date), date_list_all=date_list)
        print(f'exclude date: {ex_date_list}')

    # start_date
    if start_date:
        print(f'start   date: {start_date}')
        year_min = yyyymmdd2years(yyyymmdd(start_date))
        for year, date_str in zip(year_list, date_list):
            if year < year_min and date_str not in ex_date_list:
                print(f'  remove date: {date_str}')
                ex_date_list.append(date_str)

    # end_date
    if end_date:
        print(f'end     date: {end_date}')
        year_max = yyyymmdd2years(yyyymmdd(end_date))
        for year, date_str in zip(year_list, date_list):
            if year > year_max and date_str not in ex_date_list:
                print(f'  remove date: {date_str}')
                ex_date_list.append(date_str)

    ex_date_list = sorted(list(set(ex_date_list)))

    return ex_date_list



################################################################
def date_list2tbase(date_list, ref_date=None):
    """Get temporal Baseline in days with respect to the 1st date
    Parameters: date_list - list(str), date in YYYYMMDD or YYMMDD format
                ref_date  - str, reference date in YYYYMMDD format
    Returns:    tbase     - list(int), temporal baseline in days
                dateDict  - dict with key   - string, date in YYYYMMDD format
                                      value - int, temporal baseline in days
    """
    date_list = yyyymmdd(date_list)
    ref_date = ref_date if ref_date else date_list[0]
    ref_date = yyyymmdd(ref_date)

    # date str to dt object
    dt_fmt = get_date_str_format(str(date_list))
    dt_list = [dt.datetime.strptime(i, dt_fmt) for i in date_list]
    ref_dt = dt.datetime.strptime(ref_date, get_date_str_format(ref_date))

    # list: dt object to time difference in days
    tbase = []
    for dt_i in dt_list:
        delta_dt = dt_i - ref_dt
        tbase.append(delta_dt.days + delta_dt.seconds / (24 * 60 * 60))

    # dict: key - date, value - temporal baseline
    dateDict = {}
    for i, date_str in enumerate(date_list):
        dateDict[date_str] = tbase[i]
    return tbase, dateDict


def date_list2vector(date_list, seconds=0):
    """Get time in datetime format: datetime.datetime(2006, 5, 26, 0, 0)

    Parameters: date_list  - list of string, date in YYYYMMDD or YYMMDD format
                seconds    - float, float or str, acquisition time of the day info in seconds.
    Returns:    dates      - list of datetime.datetime objects, i.e. datetime.datetime(2010, 10, 20, 0, 0)
                datevector - list of float, years, i.e. 2010.8020547945205
    """
    date_list = yyyymmdd(date_list)
    date_format = get_date_str_format(str(date_list))
    dates = [dt.datetime.strptime(i, date_format) for i in date_list]

    # add time of the day info if:
    # 1) seconds arg is valid AND
    # 2) no time info from dates arg
    if seconds and 'T' not in date_format:
        dates = [x + dt.timedelta(seconds=float(seconds)) for x in dates]

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
def get_date_range(dmin, dmax, dstep=1, dunit='D', out_fmt='%Y%m%d'):
    """Make a list of dates with one-day (or given days) interval for [dmin, dmax]
    Parameters: dmin    - str in format supported by get_date_str_format()
                dmax    - str in format supported by get_date_str_format()
                dstep   - int, interval in number of dunit
                dunit   - str, unit of interval, e.g. Y, M, W, D, h, m, s
                out_fmt - str, output datetime string format
    Returns:    dt_list - list of str in YYYYMMDD format
    """
    # read inputs
    date_str_format = get_date_str_format(dmin)
    t1 = np.datetime64(dt.datetime.strptime(dmin, date_str_format).isoformat())
    t2 = np.datetime64(dt.datetime.strptime(dmax, date_str_format).isoformat())
    tstep = np.timedelta64(dstep, dunit)

    # prepare date range
    dt_objs = np.arange(t1, t2+tstep, tstep, dtype='datetime64').astype('O')
    dt_list = [obj.strftime(out_fmt) for obj in dt_objs]

    return dt_list


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
    from math import cos, pi, sin

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
