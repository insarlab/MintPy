#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################
# Add bootstrap method for std. dev. estimation, Emre Havazli, May 2020.
# Add poly / periodic / step func., Yuan-Kai Liu, Aug 2020.


import os
import sys
import time
import argparse
import h5py
import numpy as np
from scipy import linalg

from mintpy.objects import timeseries, giantTimeseries, HDFEOS, cluster
from mintpy.defaults.template import get_template_content
from mintpy.utils import arg_group, readfile, writefile, ptime, utils as ut


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
  timeseries2velocity.py  timeseries.h5  --exclude-date exclude_date.txt

  timeseries2velocity.py  LS-PARAMS.h5
  timeseries2velocity.py  NSBAS-PARAMS.h5
  timeseries2velocity.py  TS-PARAMS.h5

  # bootstrapping for STD calculation
  timeseries2velocity.py timeseries_ERA5_demErr.h5 --bootstrap

  # complex time functions
  timeseries2velocity.py timeseries_ERA5_ramp_demErr.h5 --poly 3 --period 1 0.5 --step 20170910
"""

DROP_DATE_TXT = """exclude_date.txt:
20040502
20060708
20090103
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Inverse velocity from time-series.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=TEMPLATE+'\n'+REFERENCE+'\n'+EXAMPLE)

    parser.add_argument('timeseries_file',
                        help='Time series file for velocity inversion.')
    parser.add_argument('--template', '-t', dest='template_file', help='template file with options')
    parser.add_argument('-o', '--output', dest='outfile', help='output file name')
    parser.add_argument('--update', dest='update_mode', action='store_true',
                        help='Enable update mode, and skip estimation if:\n'+
                             '1) output velocity file already exists, readable '+
                             'and newer than input file\n' +
                             '2) all configuration parameters are the same.')

    # reference in time and space
    # for input file without reference info, e.g. ERA5.h5
    parser.add_argument('--ref-yx', dest='ref_yx', metavar=('Y', 'X'), type=int, nargs=2,
                        help='Change referene point Y X for display')
    parser.add_argument('--ref-date', dest='ref_date', metavar='DATE',
                        help='Change reference date for display')

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

    # time functions
    model = parser.add_argument_group('deformation model', 'a suite of time functions')
    model.add_argument('--polynomial', '--poly', '--poly-order', dest='polynomial', type=int, default=1,
                      help='a polynomial function with the input degree (default: %(default)s). E.g.:\n' +
                           '--polynomial 1            # linear\n' +
                           '--polynomial 2            # quadratic\n' + 
                           '--polynomial 3            # cubic\n')
    model.add_argument('--periodic', '--peri', dest='periodic', type=float, nargs='+', default=[],
                      help='periodic function(s) with period in decimal years (default: %(default)s). E.g.:\n' +
                           '--periodic 1.0            # an annual cycle\n' +
                           '--periodic 1.0 0.5        # an annual cycle plus a semi-annual cycle\n')
    model.add_argument('--step', dest='step', type=str, nargs='+', default=[],
                      help='step function(s) at YYYYMMDD (default: %(default)s). E.g.:\n' +
                           '--step 20061014           # coseismic step  at 2006-10-14\n' +
                           '--step 20110311 20120928  # coseismic steps at 2011-03-11 and 2012-09-28\n')

    # bootstrap
    bootstrap = parser.add_argument_group('bootstrapping', 'estimating the mean / STD of the velocity estimator')
    bootstrap.add_argument('--bootstrap', '--bootstrapping', dest='bootstrap', action='store_true',
                           help='Enable bootstrapping to estimate the mean and STD of the velocity estimator.')
    bootstrap.add_argument('--bc', '--bootstrap-count', dest='bootstrapCount', type=int, default=400,
                           help='number of iterations for bootstrapping (default: %(default)s).')

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
        if (inps.polynomial != 1 or inps.periodic or inps.step):
            raise ValueError('bootstrapping currently support polynomial ONLY and ONLY with the order of 1!')

    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

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
                iDict[key] = ptime.yyyymmdd(value.replace(',', ' ').split())
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
            print('3) NOT all key configration parameters are the same: {}.'.format(configKeys))
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
        outname = 'velocity'
        if inps.key == 'giantTimeseries':
            prefix = os.path.basename(inps.timeseries_file).split('PARAMS')[0]
            outname = prefix + outname
        outname += '.h5'
        inps.outfile = outname
    return inps


def read_inps2model(inps):
    """get model info from inps"""

    model = dict()
    model['polynomial'] = inps.polynomial
    model['periodic'] = inps.periodic
    model['step'] = inps.step

    # msg
    print('estimate deformation model with the following assumed time functions:')
    for key, value in model.items():
        print('{:<10} : {}'.format(key, value))

    if 'polynomial' not in model.keys():
        raise ValueError('linear/polynomial model is NOT included! Are you sure?!')

    # number of parameters
    num_param = (model['polynomial'] + 1
                 + len(model['periodic']) * 2 
                 + len(model['step']))

    return model, num_param


############################################################################
def estimate_time_func(date_list, dis_ts, model):
    """
    Deformation model estimator, using a suite of linear, periodic, step function(s).

    Gm = d

    Parameters: date_list - list of str, dates in YYYYMMDD format
                dis_ts    - 2D np.ndarray, displacement observation in size of (num_date, num_pixel)
                model     - dict of time functions, e.g.:
                    {'polynomial' : 2,            # int, polynomial with 1 (linear), 2 (quadratic), 3 (cubic), etc.
                     'periodic'   : [1.0, 0.5],   # list of float, period(s) in years. 1.0 (annual), 0.5 (semiannual), etc.
                     'step'       : ['20061014'], # list of str, date(s) in YYYYMMDD.
                     ...
                     }
    Returns:    G         - 2D np.ndarray, design matrix           in size of (num_date, num_par)
                m         - 2D np.ndarray, parameter solution      in size of (num_par, num_pixel)
                e2        - 1D np.ndarray, sum of squared residual in size of (num_pixel,)
    """

    G = timeseries.get_design_matrix4time_func(date_list, model)

    # least squares solver
    # Opt. 1: m = np.linalg.pinv(G).dot(dis_ts)
    # Opt. 2: m = scipy.linalg.lstsq(G, dis_ts, cond=1e-15)[0]
    # Numpy is not used because it can not handle NaN value in dis_ts
    m, e2 = linalg.lstsq(G, dis_ts)[:2]

    return G, m, e2


def run_timeseries2time_func(inps):

    # basic info
    atr = readfile.read_attribute(inps.timeseries_file)
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    num_date = inps.numDate
    dates = np.array(inps.dateList)

    # get deformation model from parsers
    model, num_param = read_inps2model(inps)


    ## output preparation

    # attributes
    atr['FILE_TYPE'] = 'velocity'
    atr['UNIT'] = 'm/year'
    atr['START_DATE'] = inps.dateList[0]
    atr['END_DATE'] = inps.dateList[-1]
    atr['DATE12'] = '{}_{}'.format(inps.dateList[0], inps.dateList[-1])
    if inps.ref_yx:
        atr['REF_Y'] = inps.ref_yx[0]
        atr['REF_X'] = inps.ref_yx[1]
    if inps.ref_date:
        atr['REF_DATE'] = inps.ref_date

    # config parameter
    print('add/update the following configuration metadata:\n{}'.format(configKeys))
    for key in configKeys:
        atr[key_prefix+key] = str(vars(inps)[key])

    # instantiate output file
    layout_hdf5(inps.outfile, atr, model)


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
        box_width  = box[2] - box[0]
        box_length = box[3] - box[1]
        num_pixel = box_length * box_width
        if num_box > 1:
            print('\n------- processing patch {} out of {} --------------'.format(i+1, num_box))
            print('box width:  {}'.format(box_width))
            print('box length: {}'.format(box_length))

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
            ts_data -= np.tile(ref_val.reshape(ts_data.shape[0], 1, 1), (1, ts_data.shape[1], ts_data.shape[2]))

        ts_data = ts_data[inps.dropDate, :, :].reshape(inps.numDate, -1)
        if atr['UNIT'] == 'mm':
            ts_data *= 1./1000.

        # mask invalid pixels
        print('skip pixels with zero/nan value in all acquisitions')
        ts_stack = np.nanmean(ts_data, axis=0)
        mask = np.multiply(~np.isnan(ts_stack), ts_stack!=0.)
        del ts_stack

        ts_data = ts_data[:, mask]
        num_pixel2inv = int(np.sum(mask))
        print('number of pixels to invert: {} out of {} ({:.1f}%)'.format(
            num_pixel2inv, num_pixel, num_pixel2inv/num_pixel*100))

        # go to next if no valid pixel found
        if num_pixel2inv == 0:
            block = [box[1], box[3], box[0], box[2]]
            write_hdf5_block(inps.outfile, model, m, m_std,
                             mask=mask,
                             block=block)
            continue


        ### estimation / solve Gm = d

        if inps.bootstrap:
            ## option 1 - least squares with bootstrapping
            # Bootstrapping is a resampling method which can be used to estimate properties
            # of an estimator. The method relies on independently sampling the data set with
            # replacement.

            try:
                from sklearn.utils import resample
            except ImportError:
                raise ImportError('can not import scikit-learn!')
            print('using bootstrap resampling {} times ...'.format(inps.bootstrapCount)) 

            # calc model of all bootstrap sampling
            m_boot = np.zeros((inps.bootstrapCount, num_param, num_pixel2inv), dtype=dataType)
            prog_bar = ptime.progressBar(maxValue=inps.bootstrapCount)
            for i in range(inps.bootstrapCount):
                # bootstrap resampling
                boot_ind = resample(np.arange(inps.numDate),
                                    replace=True,
                                    n_samples=inps.numDate)
                boot_ind.sort()

                # estimation
                m_boot[i] = estimate_time_func(dates[boot_ind].tolist(),
                                               ts_data[boot_ind],
                                               model)[1]

                prog_bar.update(i+1, suffix='iteration {} / {}'.format(i+1, inps.bootstrapCount))
            prog_bar.close()
            del ts_data

            # get mean/std among all bootstrap sampling
            print('calculate mean and standard deviation of bootstrap estimations')
            m[:, mask] = m_boot.mean(axis=0).reshape(num_param, -1)
            m_std[:, mask] = m_boot.std(axis=0).reshape(num_param, -1)
            del m_boot


        else:
            ## option 2 - least squares with uncertainty propagation

            print('estimate time functions via linalg.lstsq ...')
            G, m[:, mask], e2 = estimate_time_func(inps.dateList,
                                                   ts_data,
                                                   model)
            del ts_data

            ## Compute the covariance matrix for model parameters: Gm = d
            # C_m_hat = (G.T * C_d^-1, * G)^-1  # the most generic form
            #         = sigma^2 * (G.T * G)^-1  # assuming the obs error is normally distributed in time.
            # Based on the law of integrated expectation, we estimate the obs sigma^2 using
            # the OLS estimation residual e_hat_i = d_i - d_hat_i
            # sigma^2 = sigma_hat^2 * N / (N - P)
            #         = (e_hat.T * e_hat) / (N - P)  # sigma_hat^2 = (e_hat.T * e_hat) / N

            G_inv = linalg.inv(np.dot(G.T, G))
            m_var = e2.reshape(1, -1) / (num_date - num_param)
            m_std[:, mask] = np.sqrt(np.dot(np.diag(G_inv).reshape(-1, 1), m_var))

            ## for linear velocity, the STD can also be calculated 
            # using Eq. (10) from Fattahi and Amelung (2015, JGR)
            # ts_diff = ts_data - np.dot(G, m)
            # t_diff = G[:, 1] - np.mean(G[:, 1])
            # vel_std = np.sqrt(np.sum(ts_diff ** 2, axis=0) / np.sum(t_diff ** 2)  / (num_date - 2))

        # write
        block = [box[1], box[3], box[0], box[2]]
        write_hdf5_block(inps.outfile, model, m, m_std,
                         mask=mask,
                         block=block)

    return inps.outfile


def layout_hdf5(out_file, atr, model):
    """create HDF5 file for estimated time functions
    with defined metadata and (empty) dataset structure
    """

    # deformation model info
    poly_deg = model['polynomial']
    num_period = len(model['periodic'])
    num_step = len(model['step'])

    # size info
    length = int(atr['LENGTH'])
    width = int(atr['WIDTH'])

    ds_name_dict = {}
    ds_unit_dict = {}

    # time func 1 - polynomial
    for i in range(1, poly_deg + 1):
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

        # update ds_name/unit_dict
        ds_name_dict[dsName] = [dataType, (length, width), None]
        ds_unit_dict[dsName] = unit
        ds_name_dict[dsName+'Std'] = [dataType, (length, width), None]
        ds_unit_dict[dsName+'Std'] = unit

    # time func 2 - periodic
    for i in range(num_period):
        # dataset name
        period = model['periodic'][i]
        if period == 1:
            dsName = 'annualAmp'
        elif period == 0.5:
            dsName = 'semiAnnualAmp'
        else:
            dsName = 'periodY{}Amp'.format(period)

        # update ds_name/unit_dict
        ds_name_dict[dsName] = [dataType, (length, width), None]
        ds_unit_dict[dsName] = 'm'

    # time func 3 - step
    for i in range(num_step):
        # dataset name
        dsName = 'step{}'.format(model['step'][i])

        # update ds_name/unit_dict
        ds_name_dict[dsName] = [dataType, (length, width), None]
        ds_unit_dict[dsName] = 'm'
        ds_name_dict[dsName+'Std'] = [dataType, (length, width), None]
        ds_unit_dict[dsName+'Std'] = 'm'

    # layout hdf5
    writefile.layout_hdf5(out_file, ds_name_dict, metadata=atr)

    # add metadata to HDF5 dataset
    max_digit = max([len(i) for i in ds_unit_dict.keys()])
    with h5py.File(out_file, 'r+') as f:
        for key, value in ds_unit_dict.items():
            f[key].attrs['UNIT'] = value
            print('add /{d:<{w}} attribute: UNIT = {u}'.format(d=key,
                                                               w=max_digit,
                                                               u=value))

    return out_file


def write_hdf5_block(out_file, model, m, m_std, mask=None, block=None):
    """write the estimated time function parameters to file

    Parameters: out_file - str, path of output time func file
                model    - dict, dict of time functions, e.g.:
                    {'polynomial' : 2,            # int, polynomial with
                                                  # e.g.: 1 (linear), 2 (quad), 3 (cubic), etc.
                     'periodic'   : [1.0, 0.5],   # list of float, period(s) in years.
                                                  # e.g.: 1.0 (annual), 0.5 (semiannual), etc.
                     'step'       : ['20061014'], # list of str, date(s) in YYYYMMDD.
                     ...
                     }
                m/m_std  - 2D np.ndarray in float32 in size of (num_param, length*width), time func param. (Std. Dev.)
                mask     - 1D np.ndarray in float32 in size of (length*width), mask of valid pixels
                block    - list of 4 int, for [yStart, yEnd, xStart, xEnd]
    """

    def write_dataset_block(f, dsName, data, block):
        print('write dataset /{:<20} block: {}'.format(dsName, block))
        f[dsName][block[0]:block[1], 
                  block[2]:block[3]] = data.reshape(block[1] - block[0],
                                                    block[3] - block[2])

    # deformation model info
    poly_deg = model['polynomial']
    num_period = len(model['periodic'])
    num_step = len(model['step'])

    length = block[1] - block[0]
    width = block[3] - block[2]
    if mask is None:
        mask = np.ones(length*width, dtype=np.bool_)

    print('open file: {} with "a" mode'.format(out_file))
    with h5py.File(out_file, 'a') as f:

        # time func 1 - polynomial
        for i in range(1, poly_deg+1):
            # dataset name
            if i == 1:
                dsName = 'velocity'
            elif i == 2:
                dsName = 'acceleration'
            else:
                dsName = 'poly{}'.format(i)

            # write
            write_dataset_block(f, 
                                dsName=dsName,
                                data=m[i, :],
                                block=block)
            write_dataset_block(f, 
                                dsName=dsName+'Std',
                                data=m_std[i, :],
                                block=block)

        # time func 2 - periodic
        p0 = poly_deg + 1
        for i in range(num_period):
            # calculate the amplitude and phase of the periodic signal
            # following equation (9-10) in Minchew et al. (2017, JGR)
            coef_cos = m[p0 + 2*i, :]
            coef_sin = m[p0 + 2*i + 1, :]
            period_amp = np.sqrt(coef_cos**2 + coef_sin**2)
            period_pha = np.zeros(length*width, dtype=dataType)
            period_pha[mask] = np.arctan(coef_cos[mask] / coef_sin[mask])

            # dataset name
            period = model['periodic'][i]
            if period == 1:
                dsName = 'annualAmp'
            elif period == 0.5:
                dsName = 'semiAnnualAmp'
            else:
                dsName = 'periodY{}Amp'.format(period)

            # write
            write_dataset_block(f, 
                                dsName=dsName,
                                data=period_amp,
                                block=block)

            # 1. figure out a proper way to save the phase data in radians
            #    and keeping smart display in view.py 
            #    to avoid messy scaling together with dataset in m
            # 2. add code for the error propagation for the periodic amp/pha

        # 3. step
        p0 = (poly_deg + 1) + (2 * num_period)
        for i in range(num_step):
            # dataset name
            dsName = 'step{}'.format(model['step'][i])

            # write
            write_dataset_block(f,
                                dsName=dsName,
                                data=m[p0+i, :],
                                block=block)
            write_dataset_block(f,
                                dsName=dsName+'Std',
                                data=m_std[p0+i, :],
                                block=block)

    print('close HDF5 file {}'.format(out_file))
    return out_file


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
