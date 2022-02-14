#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################
# Add bootstrap method for std. dev. estimation, Emre Havazli, May 2020.
# Add poly, periodic and step func., Yuan-Kai Liu, Aug 2020.
# Add exp and log func., Yuan-Kai Liu, Jun 2021.


import os
import sys
import time
import argparse
import numpy as np
from scipy import linalg

from mintpy.defaults.template import get_template_content
from mintpy.objects import timeseries, giantTimeseries, HDFEOS, cluster
from mintpy.utils import arg_group, ptime, time_func, readfile, writefile, utils as ut


dataType = np.float32
# key configuration parameter name
key_prefix = 'mintpy.velocity.'
configKeys = [
    'startDate',
    'endDate',
    'excludeDate',
    'bootstrap',
    'bootstrapCount',
]


############################################################################
TEMPLATE = get_template_content('velocity')

REFERENCE = """references:
  Fattahi, H., and F. Amelung (2015), InSAR bias and uncertainty due to the systematic and stochastic
  tropospheric delay, Journal of Geophysical Research: Solid Earth, 120(12), 8758-8773, doi:10.1002/2015JB012419.

  Efron, B., and R. Tibshirani (1986), Bootstrap methods for standard errors, confidence intervals,
  and other measures of statistical accuracy, Statistical science, 54-75, doi:10.1214/ss/1177013815.
"""

EXAMPLE = """example:
  timeseries2velocity.py  timeseries_ERA5_demErr.h5
  timeseries2velocity.py  timeseries_ERA5_demErr_ramp.h5  -t KyushuT73F2980_2990AlosD.template
  timeseries2velocity.py  timeseries.h5  --start-date 20080201  --end-date 20100508
  timeseries2velocity.py  timeseries.h5  --exclude exclude_date.txt

  timeseries2velocity.py  LS-PARAMS.h5
  timeseries2velocity.py  NSBAS-PARAMS.h5
  timeseries2velocity.py  TS-PARAMS.h5

  # bootstrapping for STD calculation
  timeseries2velocity.py timeseries_ERA5_demErr.h5 --bootstrap

  # complex time functions
  timeseries2velocity.py timeseries_ERA5_ramp_demErr.h5 --poly 3 --period 1 0.5 --step 20170910
  timeseries2velocity.py timeseries_ERA5_demErr.h5      --poly 1 --exp 20170910 90
  timeseries2velocity.py timeseries_ERA5_demErr.h5      --poly 1 --log 20170910 60.4
  timeseries2velocity.py timeseries_ERA5_demErr.h5      --poly 1 --log 20170910 60.4 200 --log 20171026 200.7
"""

DROP_DATE_TXT = """exclude_date.txt:
20040502
20060708
20090103
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Estimate velocity / time functions from time-series.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=TEMPLATE+'\n'+REFERENCE+'\n'+EXAMPLE)

    # inputs
    parser.add_argument('timeseries_file',
                        help='Time series file for velocity inversion.')
    parser.add_argument('--template', '-t', dest='template_file', help='template file with options')
    parser.add_argument('--ts-cov-file', dest='ts_cov_file',
                        help='Time-series (co)variance file for velocity STD calculation')

    # outputs
    parser.add_argument('-o', '--output', dest='outfile', help='output file name')
    parser.add_argument('--update', dest='update_mode', action='store_true',
                        help='Enable update mode, and skip estimation if:\n'+
                             '1) output velocity file already exists, readable '+
                             'and newer than input file\n' +
                             '2) all configuration parameters are the same.')

    # reference in time and space
    # for input file without reference info, e.g. ERA5.h5
    parser.add_argument('--ref-lalo', dest='ref_lalo', metavar=('LAT', 'LON'), type=float, nargs=2,
                        help='Change referene point LAT LON for estimation.')
    parser.add_argument('--ref-yx', dest='ref_yx', metavar=('Y', 'X'), type=int, nargs=2,
                        help='Change referene point Y X for estimation.')
    parser.add_argument('--ref-date', dest='ref_date', metavar='DATE',
                        help='Change reference date for estimation.')

    # dates of interest
    date = parser.add_argument_group('dates of interest')
    date.add_argument('--start-date','-s', dest='startDate',
                      help='start date for velocity estimation')
    date.add_argument('--end-date','-e', dest='endDate',
                      help='end date for velocity estimation')
    date.add_argument('--exclude', '--ex', dest='excludeDate', nargs='+', default=[],
                      help='date(s) not included in velocity estimation, i.e.:\n' +
                           '--exclude 20040502 20060708 20090103\n' +
                           '--exclude exclude_date.txt\n'+DROP_DATE_TXT)

    # bootstrap
    bootstrap = parser.add_argument_group('bootstrapping', 'estimating the mean / STD of the velocity estimator')
    bootstrap.add_argument('--bootstrap', '--bootstrapping', dest='bootstrap', action='store_true',
                           help='Enable bootstrapping to estimate the mean and STD of the velocity estimator.')
    bootstrap.add_argument('--bc', '--bootstrap-count', dest='bootstrapCount', type=int, default=400,
                           help='number of iterations for bootstrapping (default: %(default)s).')

    # time functions
    parser = arg_group.add_timefunc_argument(parser)

    # residual file
    resid = parser.add_argument_group('Residual file', 'Save residual displacement time-series to HDF5 file.')
    resid.add_argument('--save-res', '--save_residual', dest='save_res', action='store_true',
                       help='Save the residual displacement time-series to HDF5 file.')
    resid.add_argument('--res-file', '--residual-file', dest='res_file', default='timeseriesResidual.h5',
                       help='Output file name for the residual time-series file (default: %(default)s).')

    # computing
    parser = arg_group.add_memory_argument(parser)

    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    inps.key = readfile.read_attribute(inps.timeseries_file)['FILE_TYPE']
    if inps.key not in ['timeseries', 'giantTimeseries', 'HDFEOS']:
        raise Exception('input file is {}, NOT timeseries!'.format(inps.key))

    # check bootstrap count number
    if inps.bootstrap and inps.bootstrapCount <= 1:
        inps.bootstrap = False
        print('bootstrap-count should be larger than 1, otherwise it does not make sense')
        print('turn OFF bootstrapping and continue without it.')

    if inps.bootstrap:
        print('bootstrapping is turned ON.')
        if (inps.polynomial != 1 or inps.periodic or inps.step or inps.exp or inps.log):
            raise ValueError('bootstrapping currently support polynomial ONLY and ONLY with the order of 1!')

    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    # Initialize the dictionaries of exp and log funcs
    inps = init_exp_log_dicts(inps)

    # --ref-lalo option
    if inps.ref_lalo:
        atr = readfile.read_attribute(inps.timeseries_file)
        coord = ut.coordinate(atr)
        ref_y, ref_x = coord.geo2radar(inps.ref_lalo[0], inps.ref_lalo[1])[:2]
        if ref_y is not None and ref_x is not None:
            inps.ref_yx = [ref_y, ref_x]
            print('input reference point in (lat, lon): ({}, {})'.format(inps.ref_lalo[0], inps.ref_lalo[1]))
            print('corresponding   point in (y, x): ({}, {})'.format(inps.ref_yx[0], inps.ref_yx[1]))

    return inps


def init_exp_log_dicts(inps):
    """Initialize the dictionaries of exp and log funcs
    By trarnslating inps.exp/log into inps.expDict/logDict.
    """
    # --exp option: convert cmd inputs into dict format
    inps.expDict = dict()
    if inps.exp:
        for exp_list in inps.exp:
            onset_time, char_times = exp_list[0], exp_list[1:]
            if len(onset_time) == 8:
                if len(char_times) > 0:
                    inps.expDict[onset_time] = np.array(char_times).astype(float).tolist()

                else:
                    msg = 'NO characteristic time found: {}\n'.format(char_times)
                    msg += 'one or more characteristic time(s) are required for each onset date'
                    msg += ' for the exp function, e.g.:\n'
                    msg += '--exp 20181026 60 OR\n'
                    msg += '--exp 20161231 80.5 200  # append as many char_times as you like!'
                    raise ValueError(msg)
            else:
                raise ValueError('input onset time is NOT in YYYYMMDD format: {}'.format(onset_time))

    # --log option: convert cmd inputs into dict format
    inps.logDict = dict()
    if inps.log:
        for log_list in inps.log:
            onset_time, char_times = log_list[0], log_list[1:]
            if len(onset_time) == 8:
                if len(char_times) > 0:
                    inps.logDict[onset_time] = np.array(char_times).astype(float).tolist()

                else:
                    msg = 'NO characteristic time found: {}\n'.format(char_times)
                    msg += 'one or more characteristic time(s) are required for each onset date'
                    msg += ' for the log function, e.g.:\n'
                    msg += '--exp 20181026 60 OR\n'
                    msg += '--exp 20161231 80.5 200  # append as many char_times as you like!'
                    raise ValueError(msg)
            else:
                raise ValueError('input onset time is NOT in YYYYMMDD format: {}'.format(onset_time))

    return inps


def read_template2inps(template_file, inps=None):
    """Read input template file into inps.excludeDate"""
    if not inps:
        inps = cmd_line_parse()
    iDict = vars(inps)
    print('read options from template file: '+os.path.basename(template_file))
    template = readfile.read_template(inps.template_file)
    template = ut.check_template_auto_value(template)

    # Read template option
    prefix = 'mintpy.velocity.'
    keyList = [i for i in list(iDict.keys()) if prefix+i in template.keys()]
    for key in keyList:
        value = template[prefix+key]
        if key in ['bootstrap']:
            iDict[key] = value
        if value:
            if key in ['startDate', 'endDate']:
                iDict[key] = ptime.yyyymmdd(value)
            elif key in ['excludeDate']:
                value = value.replace('[','').replace(']','').replace(',', ' ')
                iDict[key] = ptime.yyyymmdd(value.split())
            elif key in ['bootstrapCount']:
                iDict[key] = int(value)

    key = 'mintpy.compute.maxMemory'
    if key in template.keys() and template[key]:
        inps.maxMemory = float(template[key])

    return inps


def run_or_skip(inps):
    print('update mode: ON')
    flag = 'skip'

    # check output file
    if not os.path.isfile(inps.outfile):
        flag = 'run'
        print('1) output file {} NOT found.'.format(inps.outfile))
    else:
        print('1) output file {} already exists.'.format(inps.outfile))
        ti = os.path.getmtime(inps.timeseries_file)
        to = os.path.getmtime(inps.outfile)
        if ti > to:
            flag = 'run'
            print('2) output file is NOT newer than input file: {}.'.format(inps.timeseries_file))
        else:
            print('2) output file is newer than input file: {}.'.format(inps.timeseries_file))

    # check configuration
    if flag == 'skip':
        atr = readfile.read_attribute(inps.outfile)
        if any(str(vars(inps)[key]) != atr.get(key_prefix+key, 'None') for key in configKeys):
            flag = 'run'
            print('3) NOT all key configuration parameters are the same: {}.'.format(configKeys))
        else:
            print('3) all key configuration parameters are the same: {}.'.format(configKeys))

    # result
    print('run or skip: {}.'.format(flag))
    return flag


############################################################################
def read_exclude_date(inps, dateListAll):
    # Merge ex_date/startDate/endDate into ex_date
    yy_list_all = ptime.yyyymmdd2years(dateListAll)
    exDateList = []

    # ex_date
    exDateList += ptime.read_date_list(list(inps.excludeDate), date_list_all=dateListAll)
    if exDateList:
        print('exclude date:'+str(exDateList))

    # startDate
    if inps.startDate:
        print('start date: '+inps.startDate)
        yy_min = ptime.yyyymmdd2years(ptime.yyyymmdd(inps.startDate))
        for i in range(len(dateListAll)):
            date = dateListAll[i]
            if yy_list_all[i] < yy_min and date not in exDateList:
                print('  remove date: '+date)
                exDateList.append(date)

    # endDate
    if inps.endDate:
        print('end date: '+inps.endDate)
        yy_max = ptime.yyyymmdd2years(ptime.yyyymmdd(inps.endDate))
        for i in range(len(dateListAll)):
            date = dateListAll[i]
            if yy_list_all[i] > yy_max and date not in exDateList:
                print('  remove date: '+date)
                exDateList.append(date)
    exDateList = sorted(list(set(exDateList)))
    return exDateList


def read_date_info(inps):
    """Read dates used in the estimation and its related info.
    Parameters: inps - Namespace
    Returns:    inps - Namespace
    """
    if inps.key == 'timeseries':
        tsobj = timeseries(inps.timeseries_file)
    elif inps.key == 'giantTimeseries':
        tsobj = giantTimeseries(inps.timeseries_file)
    elif inps.key == 'HDFEOS':
        tsobj = HDFEOS(inps.timeseries_file)
    tsobj.open()
    inps.excludeDate = read_exclude_date(inps, tsobj.dateList)

    # exclude dates without obs data [for offset time-series only for now]
    if os.path.basename(inps.timeseries_file).startswith('timeseriesRg'):
        date_list = timeseries(inps.timeseries_file).get_date_list()
        data, atr = readfile.read(inps.timeseries_file)
        flag = np.nansum(data, axis=(1,2)) == 0
        flag[date_list.index(atr['REF_DATE'])] = 0
        if np.sum(flag) > 0:
            print('number of empty dates to exclude: {}'.format(np.sum(flag)))
            inps.excludeDate += np.array(date_list)[flag].tolist()
            inps.excludeDate = sorted(list(set(inps.excludeDate)))

    # Date used for estimation inps.dateList
    inps.dateList = [i for i in tsobj.dateList if i not in inps.excludeDate]
    inps.numDate = len(inps.dateList)
    inps.startDate = inps.dateList[0]
    inps.endDate = inps.dateList[-1]
    print('-'*50)
    print('dates from input file: {}\n{}'.format(tsobj.numDate, tsobj.dateList))
    print('-'*50)
    if len(inps.dateList) == len(tsobj.dateList):
        print('using all dates to calculate the velocity')
    else:
        print('dates used to estimate the velocity: {}\n{}'.format(inps.numDate, inps.dateList))
    print('-'*50)

    # flag array for ts data reading
    inps.dropDate = np.array([i not in inps.excludeDate for i in tsobj.dateList], dtype=np.bool_)

    # output file name
    if not inps.outfile:
        fbase = os.path.splitext(os.path.basename(inps.timeseries_file))[0]
        outname = 'velocity'
        if inps.key == 'giantTimeseries':
            prefix = os.path.basename(inps.timeseries_file).split('PARAMS')[0]
            outname = prefix + outname
        elif fbase in ['timeseriesRg', 'timeseriesAz']:
            suffix = fbase.split('timeseries')[-1]
            outname = outname + suffix
        outname += '.h5'
        inps.outfile = outname

    return inps


def read_inps2model(inps, date_list=None):
    """get model info from inps"""
    # check model date limits
    if not date_list:
        date_list = inps.dateList
    dmin, dmax = date_list[0], date_list[-1]
    ymin = ptime.yyyymmdd2years(dmin)
    ymax = ptime.yyyymmdd2years(dmax)

    if inps.step:
        for d_step in inps.step:
            y_step = ptime.yyyymmdd2years(d_step)
            if not (ymin < y_step < ymax):
                raise ValueError(f'input step date "{d_step}" exceed date list min/max: {dmin}, {dmax}')

    if inps.expDict:
        for d_onset in inps.expDict.keys():
            y_onset = ptime.yyyymmdd2years(d_onset)
            if y_onset >= ymax:
                raise ValueError(f'input exp onset date "{d_onset}" >= the last date: {dmax}')

    if inps.logDict:
        for d_onset in inps.logDict.keys():
            y_onset = ptime.yyyymmdd2years(d_onset)
            if y_onset >= ymax:
                raise ValueError(f'input log onset date "{d_onset}" >= the last date: {dmax}')

    model = dict()
    model['polynomial'] = inps.polynomial
    model['periodic']   = inps.periodic
    model['step']       = inps.step
    model['exp']        = inps.expDict
    model['log']        = inps.logDict

    # msg
    print('estimate deformation model with the following assumed time functions:')
    for key, value in model.items():
        print('    {:<10} : {}'.format(key, value))

    if 'polynomial' not in model.keys():
        raise ValueError('linear/polynomial model is NOT included! Are you sure?!')

    # number of parameters
    num_param = (
        model['polynomial'] + 1
        + len(model['periodic']) * 2
        + len(model['step'])
        + sum([len(val) for key, val in model['exp'].items()])
        + sum([len(val) for key, val in model['log'].items()])
    )

    return model, num_param


############################################################################
def run_timeseries2time_func(inps):

    # basic info
    atr = readfile.read_attribute(inps.timeseries_file)
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    num_date = inps.numDate
    dates = np.array(inps.dateList)
    seconds = atr.get('CENTER_LINE_UTC', 0)

    # use the 1st date as reference if not found, e.g. timeseriesResidual.h5 file
    if "REF_DATE" not in atr.keys() and not inps.ref_date:
        inps.ref_date = inps.dateList[0]
        print('WARNING: No REF_DATE found in time-series file or input in command line.')
        print('  Set "--ref-date {}" and continue.'.format(inps.dateList[0]))

    # get deformation model from parsers
    model, num_param = read_inps2model(inps)


    ## output preparation

    # time_func_param: attributes
    atrV = dict(atr)
    atrV['FILE_TYPE'] = 'velocity'
    atrV['UNIT'] = 'm/year'
    atrV['START_DATE'] = inps.dateList[0]
    atrV['END_DATE'] = inps.dateList[-1]
    atrV['DATE12'] = '{}_{}'.format(inps.dateList[0], inps.dateList[-1])
    if inps.ref_yx:
        atrV['REF_Y'] = inps.ref_yx[0]
        atrV['REF_X'] = inps.ref_yx[1]
    if inps.ref_date:
        atrV['REF_DATE'] = inps.ref_date

    # time_func_param: config parameter
    print('add/update the following configuration metadata:\n{}'.format(configKeys))
    for key in configKeys:
        atrV[key_prefix+key] = str(vars(inps)[key])

    # time_func_param: instantiate output file
    ds_name_dict, ds_unit_dict = model2hdf5_dataset(model, ds_shape=(length, width))[1:]
    writefile.layout_hdf5(inps.outfile,
                          metadata=atrV,
                          ds_name_dict=ds_name_dict,
                          ds_unit_dict=ds_unit_dict)

    # timeseries_res: attributes + instantiate output file
    if inps.save_res:
        atrR = dict(atr)
        # remove REF_DATE attribute
        for key in ['REF_DATE']:
            if key in atrR.keys():
                atrR.pop(key)
        # prepare ds_name_dict manually, instead of using ref_file, to support --ex option
        date_len = len(inps.dateList[0])
        ds_name_dict = {
            "date"       : [np.dtype(f'S{date_len}'), (num_date,), np.array(inps.dateList, dtype=np.string_)],
            "timeseries" : [np.float32,               (num_date, length, width), None]
        }
        writefile.layout_hdf5(inps.res_file, ds_name_dict=ds_name_dict, metadata=atrR)
    

    ## estimation

    # calc number of box based on memory limit
    memoryAll = (num_date + num_param * 2 + 2) * length * width * 4
    if inps.bootstrap:
        memoryAll += inps.bootstrapCount * num_param * length * width * 4
    num_box = int(np.ceil(memoryAll * 3 / (inps.maxMemory * 1024**3)))
    box_list = cluster.split_box2sub_boxes(box=(0, 0, width, length),
                                           num_split=num_box,
                                           dimension='y',
                                           print_msg=True)

    # loop for block-by-block IO
    for i, box in enumerate(box_list):
        box_wid  = box[2] - box[0]
        box_len = box[3] - box[1]
        num_pixel = box_len * box_wid
        if num_box > 1:
            print('\n------- processing patch {} out of {} --------------'.format(i+1, num_box))
            print('box width:  {}'.format(box_wid))
            print('box length: {}'.format(box_len))

        # initiate output
        m = np.zeros((num_param, num_pixel), dtype=dataType)
        m_std = np.zeros((num_param, num_pixel), dtype=dataType)

        # read input
        print('reading data from file {} ...'.format(inps.timeseries_file))
        ts_data = readfile.read(inps.timeseries_file, box=box)[0]

        # referencing in time and space
        # for file w/o reference info. e.g. ERA5.h5
        if inps.ref_date:
            print('referecing to date: {}'.format(inps.ref_date))
            ref_ind = inps.dateList.index(inps.ref_date)
            ts_data -= np.tile(ts_data[ref_ind, :, :], (ts_data.shape[0], 1, 1))

        if inps.ref_yx:
            print('referencing to point (y, x): ({}, {})'.format(inps.ref_yx[0], inps.ref_yx[1]))
            ref_box = (inps.ref_yx[1], inps.ref_yx[0], inps.ref_yx[1]+1, inps.ref_yx[0]+1)
            ref_val = readfile.read(inps.timeseries_file, box=ref_box)[0]
            ts_data -= np.tile(ref_val.reshape(ts_data.shape[0], 1, 1),
                               (1, ts_data.shape[1], ts_data.shape[2]))

        ts_data = ts_data[inps.dropDate, :, :].reshape(inps.numDate, -1)
        if atrV['UNIT'] == 'mm':
            ts_data *= 1./1000.

        ts_cov = None
        if inps.ts_cov_file:
            print(f'reading time-series covariance matrix from file {inps.ts_cov_file} ...')
            ts_cov = readfile.read(inps.ts_cov_file, box=box)[0]
            if len(ts_cov.shape) == 4:
                # full covariance matrix in 4D --> 3D
                if inps.numDate < ts_cov.shape[0]:
                    ts_cov = ts_cov[inps.dropDate, :, :, :]
                    ts_cov = ts_cov[:, inps.dropDate, :, :]
                ts_cov = ts_cov.reshape(inps.numDate, inps.numDate, -1)

            elif len(ts_cov.shape) == 3:
                # diaginal variance matrix in 3D --> 2D
                if inps.numDate < ts_cov.shape[0]:
                    ts_cov = ts_cov[inps.dropDate, :, :]
                ts_cov = ts_cov.reshape(inps.numDate, -1)

            ## set zero value to a fixed small value to avoid divide by zero
            #epsilon = 1e-5
            #ts_cov[ts_cov<epsilon] = epsilon

        # mask invalid pixels
        print('skip pixels with zero/nan value in all acquisitions')
        ts_stack = np.nanmean(ts_data, axis=0)
        mask = np.multiply(~np.isnan(ts_stack), ts_stack!=0.)
        del ts_stack

        #if ts_cov is not None:
        #    print('skip pxiels with nan STD value in any acquisition')
        #    num_std_nan = np.sum(np.isnan(ts_cov), axis=0)
        #    mask *= num_std_nan == 0
        #    del num_std_nan

        ts_data = ts_data[:, mask]
        num_pixel2inv = int(np.sum(mask))
        idx_pixel2inv = np.where(mask)[0]
        print('number of pixels to invert: {} out of {} ({:.1f}%)'.format(
            num_pixel2inv, num_pixel, num_pixel2inv/num_pixel*100))

        # go to next if no valid pixel found
        if num_pixel2inv == 0:
            continue


        ### estimation / solve Gm = d
        print('estimating time functions via linalg.lstsq ...')

        if inps.bootstrap:
            ## option 1 - least squares with bootstrapping
            # Bootstrapping is a resampling method which can be used to estimate properties
            # of an estimator. The method relies on independently sampling the data set with
            # replacement.
            print('estimating time function STD with bootstrap resampling ({} times) ...'.format(
                inps.bootstrapCount))

            # calc model of all bootstrap sampling
            rng = np.random.default_rng()
            m_boot = np.zeros((inps.bootstrapCount, num_param, num_pixel2inv), dtype=dataType)
            prog_bar = ptime.progressBar(maxValue=inps.bootstrapCount)
            for i in range(inps.bootstrapCount):
                # bootstrap resampling
                boot_ind = rng.choice(inps.numDate, size=inps.numDate, replace=True)
                boot_ind.sort()

                # estimation
                m_boot[i] = time_func.estimate_time_func(
                    model=model,
                    date_list=dates[boot_ind].tolist(),
                    dis_ts=ts_data[boot_ind],
                    seconds=seconds)[1]

                prog_bar.update(i+1, suffix='iteration {} / {}'.format(i+1, inps.bootstrapCount))
            prog_bar.close()
            #del ts_data

            # get mean/std among all bootstrap sampling
            m[:, mask] = m_boot.mean(axis=0).reshape(num_param, -1)
            m_std[:, mask] = m_boot.std(axis=0).reshape(num_param, -1)
            del m_boot

            # get design matrix to calculate the residual time series
            G = time_func.get_design_matrix4time_func(inps.dateList, model=model, ref_date=inps.ref_date, seconds=seconds)


        else:
            ## option 2 - least squares with uncertainty propagation
            G, m[:, mask], e2 = time_func.estimate_time_func(
                model=model,
                date_list=inps.dateList,
                dis_ts=ts_data,
                seconds=seconds)
            #del ts_data

            ## Compute the covariance matrix for model parameters:
            #       G * m = d
            #     C_m_hat = G+ * C_d * G+.T
            #
            # For ordinary least squares estimation:
            #     G+ = (G.T * G)^-1 * G.T                       (option 2.1)
            #
            # For weighted least squares estimation:
            #          G+ = (G.T * C_d^-1 * G)^-1 * G.T * C_d^-1
            # =>  C_m_hat = (G.T * C_d^-1 * G)^-1               (option 2.2)
            #
            # Assuming normality of the observation errors (in the time domain) with a variance of sigma^2
            # we have C_d = sigma^2 * I, then the above equation is simplfied into:
            #     C_m_hat = sigma^2 * (G.T * G)^-1              (option 2.3)
            #
            # Based on the law of integrated expectation, we estimate the obs sigma^2 using
            # the OLS estimation residual as:
            #           e_hat = d - d_hat
            # =>  sigma_hat^2 = (e_hat.T * e_hat) / N
            # =>      sigma^2 = sigma_hat^2 * N / (N - P)       (option 2.4)
            #                 = (e_hat.T * e_hat) / (N - P)
            # which is the equation (10) from Fattahi and Amelung (2015, JGR)

            if ts_cov is not None:
                # option 2.1 - linear propagation from time-series (co)variance matrix
                # TO DO: save the full covariance matrix of the time function parameters
                # only the STD is saved right now
                covar_flag = True if len(ts_cov.shape) == 3 else False
                msg = 'estimating time function STD from time-serries '
                msg += 'covariance pixel-by-pixel ...' if covar_flag else 'variance pixel-by-pixel ...'
                print(msg)

                # calc the common pseudo-inverse matrix
                Gplus = linalg.pinv(G)

                # loop over each pixel
                # or use multidimension matrix multiplication
                # m_cov = Gplus @ ts_cov @ Gplus.T
                prog_bar = ptime.progressBar(maxValue=num_pixel2inv)
                for i in range(num_pixel2inv):
                    idx = idx_pixel2inv[i]

                    # cov: time-series -> time func
                    ts_covi = ts_cov[:, :, idx] if covar_flag else np.diag(ts_cov[:, idx])
                    m_cov = np.linalg.multi_dot([Gplus, ts_covi, Gplus.T])
                    m_std[:, idx] = np.sqrt(np.diag(m_cov))

                    prog_bar.update(i+1, every=200, suffix='{}/{} pixels'.format(i+1, num_pixel2inv))
                prog_bar.close()

            else:
                # option 2.3 - assume obs errors following normal dist. in time
                print('estimating time function STD from time-series fitting residual ...')
                G_inv = linalg.inv(np.dot(G.T, G))
                m_var = e2.reshape(1, -1) / (num_date - num_param)
                m_std[:, mask] = np.sqrt(np.dot(np.diag(G_inv).reshape(-1, 1), m_var))

                # option 2.4 - simplified form for linear velocity (without matrix linear algebra)
                # The STD can also be calculated using Eq. (10) from Fattahi and Amelung (2015, JGR)
                # ts_diff = ts_data - np.dot(G, m)
                # t_diff = G[:, 1] - np.mean(G[:, 1])
                # vel_std = np.sqrt(np.sum(ts_diff ** 2, axis=0) / np.sum(t_diff ** 2)  / (num_date - 2))

        # write - time func params
        block = [box[1], box[3], box[0], box[2]]
        ds_dict = model2hdf5_dataset(model, m, m_std, mask=mask)[0]
        for ds_name, data in ds_dict.items():
            writefile.write_hdf5_block(inps.outfile,
                                       data=data.reshape(box_len, box_wid),
                                       datasetName=ds_name,
                                       block=block)

        # write - residual file
        if inps.save_res:
            block = [0, num_date, box[1], box[3], box[0], box[2]]
            ts_res = np.ones((num_date, box_len*box_wid), dtype=np.float32) * np.nan
            ts_res[:, mask] = ts_data - np.dot(G, m)[:, mask]
            writefile.write_hdf5_block(inps.res_file,
                                       data=ts_res.reshape(num_date, box_len, box_wid),
                                       datasetName='timeseries',
                                       block=block)

    return inps.outfile


def model2hdf5_dataset(model, m=None, m_std=None, mask=None, ds_shape=None):
    """Prepare the estimated model parameters into a list of dicts for HDF5 dataset writing.
    Parameters: model        - dict,
                m            - 2D np.ndarray in (num_param, num_pixel) where num_pixel = 1 or length * width
                m_std        - 2D np.ndarray in (num_param, num_pixel) where num_pixel = 1 or length * width
                mask         - 1D np.ndarray in (num_pixel), mask of valid pixels
                ds_shape     - tuple of 2 int in (length, width)
    Returns:    ds_dict      - dict, dictionary of dataset values,     input for writefile.write_hdf5_block()
                ds_name_dict - dict, dictionary of dataset initiation, input for writefile.layout_hdf5()
                ds_unit_dict - dict, dictionary of dataset unit,       input for writefile.layout_hdf5()
    Examples:   # read input model parameters into dict
                model = read_inps2model(inps, date_list=inps.date_list)
                # for time series cube
                ds_name_dict, ds_name_dict = model2hdf5_dataset(model, ds_shape=(200,300))[1:]
                ds_dict = model2hdf5_dataset(model, m, m_std, mask=mask)[0]
                # for time series point
                ds_unit_dict = model2hdf5_dataset(model)[2]
                ds_dict = model2hdf5_dataset(model, m, m_std)[0]
    """
    # deformation model info
    poly_deg   = model['polynomial']
    num_period = len(model['periodic'])
    num_step   = len(model['step'])
    num_exp    = sum([len(val) for key, val in model['exp'].items()])

    # init output
    ds_dict = {}
    ds_name_dict = {}
    ds_unit_dict = {}

    # assign ds_dict ONLY IF m is not None
    if m is not None:
        num_pixel = m.shape[1] if m.ndim > 1 else 1
        m = m.reshape(-1, num_pixel)
        m_std = m_std.reshape(-1, num_pixel)

        # default mask
        if mask is None:
            mask = np.ones(num_pixel, dtype=np.bool_)
        else:
            mask = mask.flatten()

    # time func 1 - polynomial
    for i in range(1, poly_deg+1):
        # dataset name
        if i == 1:
            dsName = 'velocity'
            unit = 'm/year'
        elif i == 2:
            dsName = 'acceleration'
            unit = 'm/year^2'
        else:
            dsName = 'poly{}'.format(i)
            unit = 'm/year^{}'.format(i)

        # assign ds_dict
        if m is not None:
            ds_dict[dsName] = m[i, :]
            ds_dict[dsName+'Std'] = m_std[i, :]

        # assign ds_name/unit_dict
        ds_name_dict[dsName] = [dataType, ds_shape, None]
        ds_unit_dict[dsName] = unit
        ds_name_dict[dsName+'Std'] = [dataType, ds_shape, None]
        ds_unit_dict[dsName+'Std'] = unit

    # time func 2 - periodic
    p0 = poly_deg + 1
    for i in range(num_period):
        # dataset name
        period = model['periodic'][i]
        dsNameSuffixes = ['Amplitude', 'Phase']
        if period == 1:
            dsNames = [f'annual{x}' for x in dsNameSuffixes]
        elif period == 0.5:
            dsNames = [f'semiAnnual{x}' for x in dsNameSuffixes]
        else:
            dsNames = [f'periodY{period}{x}' for x in dsNameSuffixes]

        # calculate the amplitude and phase of the periodic signal
        # following equation (9-10) in Minchew et al. (2017, JGR)
        if m is not None:
            coef_cos = m[p0 + 2*i, :]
            coef_sin = m[p0 + 2*i + 1, :]
            period_amp = np.sqrt(coef_cos**2 + coef_sin**2)
            period_pha = np.zeros(num_pixel, dtype=dataType)
            # avoid divided by zero warning
            if not np.all(coef_sin[mask] == 0):
                period_pha[mask] = np.arctan(coef_cos[mask] / coef_sin[mask])

            # assign ds_dict
            for dsName, data in zip(dsNames, [period_amp, period_pha]):
                ds_dict[dsName] = data

        # update ds_name/unit_dict
        ds_name_dict[dsNames[0]] = [dataType, ds_shape, None]
        ds_unit_dict[dsNames[0]] = 'm'
        ds_name_dict[dsNames[1]] = [dataType, ds_shape, None]
        ds_unit_dict[dsNames[1]] = 'radian'

    # time func 3 - step
    p0 = (poly_deg + 1) + (2 * num_period)
    for i in range(num_step):
        # dataset name
        dsName = 'step{}'.format(model['step'][i])

        # assign ds_dict
        if m is not None:
            ds_dict[dsName] = m[p0+i, :]
            ds_dict[dsName+'Std'] = m_std[p0+i, :]

        # assign ds_name/unit_dict
        ds_name_dict[dsName] = [dataType, ds_shape, None]
        ds_unit_dict[dsName] = 'm'
        ds_name_dict[dsName+'Std'] = [dataType, ds_shape, None]
        ds_unit_dict[dsName+'Std'] = 'm'

    # time func 4 - exponential
    p0 = (poly_deg + 1) + (2 * num_period) + (num_step)
    i = 0
    for exp_onset in model['exp'].keys():
        for exp_tau in model['exp'][exp_onset]:
            # dataset name
            dsName = 'exp{}Tau{}'.format(exp_onset, exp_tau)

            # assign ds_dict
            if m is not None:
                ds_dict[dsName] = m[p0+i, :]
                ds_dict[dsName+'Std'] = m_std[p0+i, :]

            # assign ds_name/unit_dict
            ds_name_dict[dsName] = [dataType, ds_shape, None]
            ds_unit_dict[dsName] = 'm'
            ds_name_dict[dsName+'Std'] = [dataType, ds_shape, None]
            ds_unit_dict[dsName+'Std'] = 'm'

            # loop because each onset_time could have multiple char_time
            i += 1

    # time func 5 - logarithmic
    p0 = (poly_deg + 1) + (2 * num_period) + (num_step) + (num_exp)
    i = 0
    for log_onset in model['log'].keys():
        for log_tau in model['log'][log_onset]:
            # dataset name
            dsName = 'log{}Tau{}'.format(log_onset, log_tau)

            # assign ds_dict
            if m is not None:
                ds_dict[dsName] = m[p0+i, :]
                ds_dict[dsName+'Std'] = m_std[p0+i, :]

            # assign ds_name/unit_dict
            ds_name_dict[dsName] = [dataType, ds_shape, None]
            ds_unit_dict[dsName] = 'm'
            ds_name_dict[dsName+'Std'] = [dataType, ds_shape, None]
            ds_unit_dict[dsName+'Std'] = 'm'

            # loop because each onset_time could have multiple char_time
            i += 1

    return ds_dict, ds_name_dict, ds_unit_dict


############################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    start_time = time.time()

    inps = read_date_info(inps)

    # --update option
    if inps.update_mode and run_or_skip(inps) == 'skip':
        return inps.outfile

    run_timeseries2time_func(inps)

    # time info
    m, s = divmod(time.time()-start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs.'.format(m, s))

    return inps.outfile


############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
