"""Utilities for time functions."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Yuan-Kai Liu, Aug 2021             #
############################################################
# Recommend import:
#   from mintpy.utils import time_func


import numpy as np
from scipy import linalg

from mintpy.utils import ptime

MODEL_EXAMPLE = """time function configuration:
    model = {
        'polynomial' : 2,                    # int, polynomial degree with 1 (linear), 2 (quadratic), 3 (cubic), etc.
        'periodic'   : [1.0, 0.5],           # list(float), period(s) in years. 1.0 (annual), 0.5 (semiannual), etc.
        'stepDate'   : ['20061014'],         # list(str), date(s) for the onset of step in YYYYMMDD.
        'polyline'   : ['20190101'],         # list(str), date(s) for the onset of extra line segments in YYYYMMDD.
        'exp'        : {'20181026': [60],    # dict, key for onset time in YYYYMMDD(THHMM) and value for char times in integer days.
                        ...
                        },
        'log'        : {'20161231': [80],    # dict, key for onset time in YYYYMMDD(THHMM) and value for char times in integer days.
                        '20190125': [100, 200],
                        ...
                        },
        ...
    }
"""


def estimate_time_func(model, date_list, dis_ts, ref_date=None, seconds=0):
    """Deformation model estimator, using a suite of linear, periodic, step, exponential, and logarithmic function(s).

    Problem setup:
        Gm = d

    Parameters: model     - dict, time functions config, e.g. {cfg}
                date_list - list of str, dates in YYYYMMDD format
                dis_ts    - 2D np.ndarray, displacement observation in size of (num_date, num_pixel)
                ref_date  - reference date from date_list
                seconds   - float or str, acquisition time of the day info in seconds.
    Returns:    G         - 2D np.ndarray, design matrix           in size of (num_date, num_param)
                m         - 2D np.ndarray, parameter solution      in size of (num_param, num_pixel)
                e2        - 1D np.ndarray, sum of squared residual in size of (num_pixel,)
    """.format(cfg=MODEL_EXAMPLE)

    G = get_design_matrix4time_func(date_list, model, ref_date=ref_date, seconds=seconds)

    # least squares solver
    # Opt. 1: m = np.linalg.pinv(G).dot(dis_ts)
    # Opt. 2: m = scipy.linalg.lstsq(G, dis_ts, cond=1e-15)[0]
    # Numpy is not used because it can not handle NaN value in dis_ts
    m, e2 = linalg.lstsq(G, dis_ts, cond=None)[:2]

    # check empty e2 due to the rank-deficient G matrix for sigularities.
    e2 = np.array(e2)
    if e2.size == 0:
        print('\nWarning: empty e2 residues array due to a redundant or rank-deficient G matrix. This can cause sigularities.')
        print('Please check: https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.lstsq.html#scipy.linalg.lstsq')
        print('The issue may be due to:')
        print('\t1) very small char time(s) or tau(s) of the exp/log function(s)')
        print('\t2) the onset time(s) of exp/log are far earlier than the minimum date of the time series.')
        print('Try a different char time, onset time.')
        print('Your G matrix of the temporal model: \n', G)
        raise ValueError('G matrix is redundant/rank-deficient!')

    return G, m, e2


def inps2model(inps, date_list=None, print_msg=True):
    """Convert time function inputs from namespace (inps) into dict object.

    Parameters: inps      - namespace, parsed time function inputs. E.g.:
                            Namespace(
                                date_list=['20141213', '20141225', ...]   # optional
                                polynomial=1,
                                periodic=[1.0, 0.5],
                                stepDate=['20110311', '20120928T1733'],
                                polyline=['20190101'],
                                exp=[['20170910', '60', '200'], ['20171026', '200']],
                                log=[['20170910', '60', '200'], ['20171026', '200']],
                            )
                date_list - list of str, date in YYYYMMDD(THHMM) format
    Returns:    model     - dict, time function configuration, e.g. {cfg}
    """.format(cfg=MODEL_EXAMPLE)

    if not date_list:
        if hasattr(inps, 'dateList'):
            date_list = inps.date_list
        else:
            raise ValueError('date_list NOT given or found in inps!')

    dmin, dmax = date_list[0], date_list[-1]
    ymin = ptime.yyyymmdd2years(dmin)
    ymax = ptime.yyyymmdd2years(dmax)

    ## inps --> model
    model = dict()
    model['polynomial'] = inps.polynomial
    model['periodic']   = inps.periodic
    model['stepDate']   = inps.stepDate
    model['polyline']   = inps.polyline

    if inps.periodic:
        # check 1 - positive value
        if any(x <= 0 for x in inps.periodic):
            raise ValueError(f'Zero or negative input period ({inps.periodic}) found!')

    if inps.stepDate:
        # check 1 - min/max limit
        for d_step in inps.stepDate:
            if not (ymin < ptime.yyyymmdd2years(d_step) < ymax):
                raise ValueError(f'input step date ({d_step}) exceeds date limit: ({dmin} / {dmax})!')

    if inps.polyline:
        # check 1 - min/max limit
        for d_start in inps.polyline:
            if not (ymin < ptime.yyyymmdd2years(d_start) < ymax):
                raise ValueError(f'input polyline date ({d_start}) exceeds date limit: ({dmin} / {dmax})!')

    for func_name, strs_list in zip(['exp', 'log'], [inps.exp, inps.log]):
        func_dict = dict()

        if strs_list:
            for strs in strs_list:
                onset_time, char_times = strs[0], strs[1:]
                char_times = [float(x) for x in char_times]

                # check 1 - onset_time - str format
                date_fmt = ptime.get_date_str_format(onset_time)
                if date_fmt not in ['%Y%m%d', '%Y%m%dT%H%M']:
                    raise ValueError(f'input onset time {onset_time} is NOT in YYYYMMDD(THHMM) format!')

                # check 2 - onset_time - min/max limit
                if ptime.yyyymmdd2years(onset_time) >= ymax:
                    raise ValueError(f'input exp onset date ({onset_time}) >= the last date: {dmax}')

                # check 3 - char_time - input format
                if len(char_times) == 0:
                    msg = 'NO characteristic time inputs found!\n'
                    msg += '1+ characteristic time(s) are required for each onset date'
                    msg += f' for the {func_name} function, e.g.:\n'
                    msg += f'--{func_name} 20181026 60 OR\n'
                    msg += f'--{func_name} 20161231 80 200   # append as many char_times as you like!'
                    raise ValueError(msg)

                elif any(x <= 0 for x in char_times):
                    raise ValueError(f'Zero or negative characteristic time ({char_times}) found!')

                else:
                    int_char_times = [int(x) for x in char_times]
                    for int_char_time, char_time in zip(int_char_times, char_times):
                        if int_char_time != char_time:
                            msg = 'WARNING: float-format characteristic time detected.'
                            msg += '\nIgnore the fraction part of a day and continue: '
                            msg += f'{char_time} --> {int_char_time} days.'
                            print(msg)

                    # save as dict
                    func_dict[onset_time] = int_char_times

        # save to model
        model[func_name] = func_dict

    ## print out model summary
    if print_msg:
        print('estimate deformation model with the following assumed time functions:')
        for key, value in model.items():
            print(f'    {key:<10} : {value}')

    # warning if no polynomial found
    if 'polynomial' not in model.keys():
        raise ValueError('linear/polynomial model is NOT included! Are you sure?!')

    return model


def get_num_param(model):
    """Get the number of unknown parameters from the given time function configuration.

    Parameters: model     - dict, time functions config, e.g. {cfg}
    Returns:    num_param - int, number of unknown parameters
    """.format(cfg=MODEL_EXAMPLE)

    num_param = (
        model['polynomial'] + 1
        + len(model['periodic']) * 2
        + len(model['stepDate'])
        + len(model['polyline'])
        + sum(len(val) for key, val in model['exp'].items())
        + sum(len(val) for key, val in model['log'].items())
    )

    return num_param


#################################### Design Matrices ##########################################

def get_design_matrix4time_func(date_list, model=None, ref_date=None, seconds=0):
    """Design matrix (function model) for time functions parameter estimation.

    Parameters: date_list - list of str in YYYYMMDD format, size=num_date
                model     - dict of time functions, e.g. {cfg}
                ref_date  - reference date from date_list
                seconds   - float or str, acquisition time of the day info in seconds.
    Returns:    A         - 2D np.ndarray of design matrix in size of (num_date, num_param)
                            num_param = (poly_deg + 1) + 2*len(periodic) + len(steps) + len(exp_taus) + len(log_taus)
    """.format(cfg=MODEL_EXAMPLE)

    ## prepare time info
    # convert list of date into array of years in float
    yr_diff = np.array(ptime.yyyymmdd2years(date_list, seconds=seconds))

    # reference date
    if ref_date is None:
        ref_date = date_list[0]
    yr_diff -= yr_diff[date_list.index(ref_date)]

    ## construct design matrix A
    # default model value
    if not model:
        model = {'polynomial' : 1}

    # read the models
    poly_deg   = model.get('polynomial', 0)
    periods    = model.get('periodic', [])
    steps      = model.get('stepDate', [])
    polylines  = model.get('polyline', [])
    exps       = model.get('exp', dict())
    logs       = model.get('log', dict())
    num_period = len(periods)
    num_step   = len(steps)
    num_pline  = len(polylines)
    num_exp    = sum(len(val) for key, val in exps.items())
    num_log    = sum(len(val) for key, val in logs.items())

    num_param = (poly_deg + 1) + (2 * num_period) + num_step + num_pline + num_exp + num_log
    if num_param <= 1:
        raise ValueError('NO time functions specified!')

    # initialize the design matrix
    num_date = len(yr_diff)
    A = np.zeros((num_date, num_param), dtype=np.float32)
    c0 = 0

    # update linear/polynomial term(s)
    # poly_deg of 0 --> offset
    # poly_deg of 1 --> velocity
    # ...
    c1 = c0 + poly_deg + 1
    A[:, c0:c1] = get_design_matrix4polynomial_func(yr_diff, poly_deg)
    c0 = c1

    # update periodic term(s)
    if num_period > 0:
        c1 = c0 + 2 * num_period
        A[:, c0:c1] = get_design_matrix4periodic_func(yr_diff, periods)
        c0 = c1

    # update coseismic/step term(s)
    if num_step > 0:
        c1 = c0 + num_step
        A[:, c0:c1] = get_design_matrix4step_func(date_list, steps, seconds=seconds)
        c0 = c1

    # update polyline term(s)
    if num_pline > 0:
        c1 = c0 + num_pline
        A[:, c0:c1] = get_design_matrix4polyline(date_list, polylines, seconds=seconds)
        c0 = c1

    # update exponential term(s)
    if num_exp > 0:
        c1 = c0 + num_exp
        A[:, c0:c1] = get_design_matrix4exp_func(date_list, exps, seconds=seconds)
        c0 = c1

    # update logarithmic term(s)
    if num_log > 0:
        c1 = c0 + num_log
        A[:, c0:c1] = get_design_matrix4log_func(date_list, logs, seconds=seconds)
        c0 = c1

    return A



def get_design_matrix4polynomial_func(yr_diff, degree):
    """design matrix/function model of linear/polynomial velocity estimation

    d = c0 + c1 * t + 1/2 * c2 * t^2 + 1/6 * c3 * t^3 + ...

    The k! denominator makes the estimated polynomial coefficient (c_k) physically meaningful:
        k=0 makes c0 the offset;
        k=1 makes c1 the velocity;
        k=2 makes c2 the acceleration;
        k=3 makes c3 the acceleration rate;

    Parameters: yr_diff - time difference from ref_date in decimal years
                degree  - polynomial models: 0=offset, 1=linear, 2=quadratic, 3=cubic, etc.
    Returns:    A       - 2D np.ndarray of poly-coeff. in size of (num_date, degree+1)
    """
    A = np.zeros([len(yr_diff), degree + 1], dtype=np.float32)
    for i in range(degree+1):
        A[:,i] = (yr_diff**i) / np.math.factorial(i)

    return A


def get_design_matrix4periodic_func(yr_diff, periods):
    """design matrix/function model of periodic velocity estimation.

    Parameters: yr_diff - 1D array of time difference from ref_date in decimal years
                periods - list of period in years: 1=annual, 0.5=semiannual, etc.
    Returns:    A       - 2D np.ndarray of periodic sine & cosine coeff. in size of (num_date, 2*num_period)
    """
    num_date = len(yr_diff)
    num_period = len(periods)
    A = np.zeros((num_date, 2*num_period), dtype=np.float32)

    for i, period in enumerate(periods):
        c0, c1 = 2*i, 2*i+1
        A[:, c0] = np.cos(2*np.pi/period * yr_diff)
        A[:, c1] = np.sin(2*np.pi/period * yr_diff)

    return A


def get_design_matrix4step_func(date_list, step_date_list, seconds=0):
    """design matrix/function model of coseismic velocity estimation.

    Parameters: date_list      - list of dates in YYYYMMDD format
                step_date_list - Heaviside step function(s) with date in YYYYMMDD
    Returns:    A              - 2D np.ndarray of 1 & 0 in size of (num_date, num_step)
    """
    num_date = len(date_list)
    num_step = len(step_date_list)
    A = np.zeros((num_date, num_step), dtype=np.float32)

    t = np.array(ptime.yyyymmdd2years(date_list, seconds=seconds))
    t_steps = ptime.yyyymmdd2years(step_date_list)
    for i, t_step in enumerate(t_steps):
        A[:, i] = np.array(t > t_step).flatten()

    return A


def get_design_matrix4polyline(date_list, start_date_list, seconds=0):
    """design matrix/function model of polyline (polygonal chain)

    The polyline can be described as:
    d = c + v * t,                                   for t <= t1
      =     ...     + p1 * (t - t1),                 for t1 < t <= t2
      =     ...     +      ...      + p2 * (t - t2), for t2 < t ...
      ...

    Parameters: date_list       - list of dates in YYYYMMDD format
                start_date_list - str, start date(s) for extra line segments in YYYYMMDD format
    Returns:    A               - 2D np.ndarray in size of (num_date, 3)
    """
    # str --> float in years
    t = np.array(ptime.yyyymmdd2years(date_list, seconds=seconds))
    t_start_list = ptime.yyyymmdd2years(start_date_list)

    # construct the design matrix
    num_date = len(t)
    num_start = len(t_start_list)
    A = np.zeros((num_date, num_start), dtype=np.float32)
    for i, t_start in enumerate(t_start_list):
        tbase = t - t_start
        tbase[tbase < 0] = 0
        A[:, i] = tbase

    return A


def get_design_matrix4exp_func(date_list, exp_dict, seconds=0):
    """design matrix/function model of exponential postseismic relaxation estimation

    Reference: Eq. (5) in Hetland et al. (2012, JGR).
    Note that there is a typo in the paper for this equation, based on the MInTS code, it should be:
        Sum_i{ a_i * H(t-Ti) * [1 - e^(-(t-T_i)/tau_i)] }
    instead of the one below shown in the paper:
        Sum_i{ a_i * H(t-Ti) * [1 - e^(-(t)/tau_i)] }
    where:
        a_i         amplitude      of i-th exp term
        T_i         onset time     of i-th exp term
        tau_i       char time      of i-th exp term (relaxation time)
        H(t-T_i)    Heaviside func of i-th exp term (ensuring the exp func is one-sided)

    Parameters: date_list - list of dates in YYYYMMDD format
                exp_dict  - dict of exp func(s) info as:
                            {{onset_time1} : [{char_time11,...,char_time1N}],
                             {onset_time2} : [{char_time21,...,char_time2N}],
                             ...
                             }
                            where onset_time is string  in YYYYMMDD format and
                                  char_time  is float32 in decimal days
    Returns:    A         - 2D np.ndarray of zeros & ones in size of (num_date, num_exp)
    """
    num_date = len(date_list)
    num_exp  = sum(len(val) for key, val in exp_dict.items())
    A = np.zeros((num_date, num_exp), dtype=np.float32)

    t = np.array(ptime.yyyymmdd2years(date_list, seconds=seconds))
    # loop for onset time(s)
    i = 0
    for exp_onset in exp_dict.keys():
        # convert string to float in years
        exp_T = ptime.yyyymmdd2years(exp_onset)

        # loop for charateristic time(s)
        for exp_tau in exp_dict[exp_onset]:
            # convert time from days to years
            exp_tau /= 365.25
            A[:, i] = np.array(t > exp_T).flatten() * (1 - np.exp(-1 * (t - exp_T) / exp_tau))
            i += 1

    return A


def get_design_matrix4log_func(date_list, log_dict, seconds=0):
    """design matrix/function model of logarithmic postseismic relaxation estimation

    Reference: Eq. (4) in Hetland et al. (2012, JGR)
    Note that there is a typo in the paper for this equation, based on the MInTS code, it should be:
        Sum_i{ a_i * H(t-Ti) * [1 + log((t-T_i)/tau_i)] }
    instead of the one below shown in the paper:
        Sum_i{ a_i * H(t-Ti) * [1 + log((t)/tau_i)] }
    where:
        a_i         amplitude      of i-th log term
        T_i         onset time     of i-th log term
        tau_i       char time      of i-th log term (relaxation time)
        H(t-T_i)    Heaviside func of i-th log term (ensuring the log func is one-sided)

    Parameters: date_list - list of dates in YYYYMMDD format
                log_dict  - dict of log func(s) info as:
                            {{onset_time1} : [{char_time11,...,char_time1N}],
                             {onset_time2} : [{char_time21,...,char_time2N}],
                             ...
                             }
                            where onset_time is string  in YYYYMMDD format and
                                  char_time  is float32 in decimal days
    Returns:    A         - 2D np.ndarray of zeros & ones in size of (num_date, num_log)
    """
    num_date = len(date_list)
    num_log  = sum(len(log_dict[x]) for x in log_dict)
    A = np.zeros((num_date, num_log), dtype=np.float32)

    t = np.array(ptime.yyyymmdd2years(date_list, seconds=seconds))
    # loop for onset time(s)
    i = 0
    for log_onset in log_dict.keys():
        # convert string to float in years
        log_T = ptime.yyyymmdd2years(log_onset)

        # loop for charateristic time(s)
        for log_tau in log_dict[log_onset]:
            # convert time from days to years
            log_tau /= 365.25

            olderr = np.seterr(invalid='ignore', divide='ignore')
            A[:, i] = np.array(t > log_T).flatten() * np.nan_to_num(
                np.log(1 + (t - log_T) / log_tau),
                nan=0,
                neginf=0,
            )
            np.seterr(**olderr)
            i += 1

    return A
