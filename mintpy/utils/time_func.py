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


def estimate_time_func(model, date_list, dis_ts, ref_date=None, seconds=0):
    """Deformation model estimator, using a suite of linear, periodic, step, exponential, and logarithmic function(s).

    Problem setup:
        Gm = d

    Parameters: model     - dict of time functions, e.g.:
                            {'polynomial' : 2,                    # int, polynomial degree with 1 (linear), 2 (quadratic), 3 (cubic), etc.
                             'periodic'   : [1.0, 0.5],           # list of float, period(s) in years. 1.0 (annual), 0.5 (semiannual), etc.
                             'step'       : ['20061014'],         # list of str, date(s) in YYYYMMDD.
                             'exp'        : {'20181026': [60],    # dict, key for onset time in YYYYMMDD and value for char times in days.
                                             ...
                                            },
                             'log'        : {'20161231': [80.5],  # dict, key for onset time in YYYYMMDD and value for char times in days.
                                             '20190125': [100, 200],
                                             ...
                                            },
                             ...
                             }
                date_list - list of str, dates in YYYYMMDD format
                dis_ts    - 2D np.ndarray, displacement observation in size of (num_date, num_pixel)
                ref_date  - reference date from date_list
                seconds   - float or str, acquisition time of the day info in seconds.
    Returns:    G         - 2D np.ndarray, design matrix           in size of (num_date, num_param)
                m         - 2D np.ndarray, parameter solution      in size of (num_param, num_pixel)
                e2        - 1D np.ndarray, sum of squared residual in size of (num_pixel,)
    """

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



#################################### Design Matrices ##########################################

def get_design_matrix4time_func(date_list, model=None, ref_date=None, seconds=0):
    """Design matrix (function model) for time functions parameter estimation.
    Parameters: date_list - list of str in YYYYMMDD format, size=num_date
                model     - dict of time functions, e.g.:
                            {'polynomial' : 2,                    # int, polynomial degree with 1 (linear), 2 (quadratic), 3 (cubic), etc.
                             'periodic'   : [1.0, 0.5],           # list of float, period(s) in years. 1.0 (annual), 0.5 (semiannual), etc.
                             'step'       : ['20061014'],         # list of str, date(s) in YYYYMMDD.
                             'exp'        : {'20181026': [60],    # dict, key for onset time in YYYYMMDD and value for char. times in days.
                                             ...
                                             },
                             'log'        : {'20161231': [80.5],  # dict, key for onset time in YYYYMMDD and value for char. times in days.
                                             '20190125': [100, 200],
                                             ...
                                             },
                             ...
                             }
                ref_date  - reference date from date_list
                seconds   - float or str, acquisition time of the day info in seconds.
    Returns:    A         - 2D array of design matrix in size of (num_date, num_param)
                            num_param = (poly_deg + 1) + 2*len(periodic) + len(steps) + len(exp_taus) + len(log_taus)
    """

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
    steps      = model.get('step', [])
    exps       = model.get('exp', dict())
    logs       = model.get('log', dict())
    num_period = len(periods)
    num_step   = len(steps)
    num_exp    = sum([len(val) for key, val in exps.items()])
    num_log    = sum([len(val) for key, val in logs.items()])

    num_param = (poly_deg + 1) + (2 * num_period) + num_step + num_exp + num_log
    if num_param <= 1:
        raise ValueError('NO time functions specified!')

    # initialize the design matrix
    num_date = len(yr_diff)
    A = np.zeros((num_date, num_param), dtype=np.float32)
    c0 = 0

    # update linear/polynomial term(s)
    if poly_deg > 0:
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

    The k! denominator makes the estimated polynomial coefficient (c_k) physically meaningful:
        k=1 makes c_1 the velocity;
        k=2 makes c_2 the acceleration;
        k=3 makes c_3 the acceleration rate;

    Parameters: yr_diff: time difference from ref_date in decimal years
                degree : polynomial models: 1=linear, 2=quadratic, 3=cubic, etc.
    Returns:    A      : 2D array of poly-coeff. in size of (num_date, degree+1)
    """
    A = np.zeros([len(yr_diff), degree + 1], dtype=np.float32)
    for i in range(degree+1):
        A[:,i] = (yr_diff**i) / np.math.factorial(i)

    return A


def get_design_matrix4periodic_func(yr_diff, periods):
    """design matrix/function model of periodic velocity estimation
    Parameters: yr_diff : 1D array of time difference from ref_date in decimal years
                periods : list of period in years: 1=annual, 0.5=semiannual, etc.
    Returns:    A       : 2D array of periodic sine & cosine coeff. in size of (num_date, 2*num_period)
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
    """design matrix/function model of coseismic velocity estimation
    Parameters: date_list      : list of dates in YYYYMMDD format
                step_date_list : Heaviside step function(s) with date in YYYYMMDD
    Returns:    A              : 2D array of zeros & ones in size of (num_date, num_step)
    """
    num_date = len(date_list)
    num_step = len(step_date_list)
    A = np.zeros((num_date, num_step), dtype=np.float32)

    t = np.array(ptime.yyyymmdd2years(date_list, seconds=seconds))
    t_steps = ptime.yyyymmdd2years(step_date_list)
    for i, t_step in enumerate(t_steps):
        A[:, i] = np.array(t > t_step).flatten()

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

    Parameters: date_list      : list of dates in YYYYMMDD format
                exp_dict       : dict of exp func(s) info as:
                                 {{onset_time1} : [{char_time11,...,char_time1N}],
                                  {onset_time2} : [{char_time21,...,char_time2N}],
                                  ...
                                  }
                                 where onset_time is string  in YYYYMMDD format and
                                       char_time  is float32 in decimal days
    Returns:    A              : 2D array of zeros & ones in size of (num_date, num_exp)
    """
    num_date = len(date_list)
    num_exp  = sum([len(val) for key, val in exp_dict.items()])
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

    Parameters: date_list      : list of dates in YYYYMMDD format
                log_dict       : dict of log func(s) info as:
                                 {{onset_time1} : [{char_time11,...,char_time1N}],
                                  {onset_time2} : [{char_time21,...,char_time2N}],
                                  ...
                                  }
                                 where onset_time is string  in YYYYMMDD format and
                                       char_time  is float32 in decimal days
    Returns:    A              : 2D array of zeros & ones in size of (num_date, num_log)
    """
    num_date = len(date_list)
    num_log  = sum([len(log_dict[x]) for x in log_dict])
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
            A[:, i] = np.array(t > log_T).flatten() * np.nan_to_num(np.log(1 + (t - log_T) / log_tau), nan=0, neginf=0)
            np.seterr(**olderr)
            i += 1

    return A

