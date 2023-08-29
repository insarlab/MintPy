"""Miscellaneous utilities for simulation."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2018                               #
############################################################
# Recommend usage:
#   from mintpy.simulation import simulation as sim


import random

import matplotlib.pyplot as plt
import numpy as np

from mintpy.defaults.plot import *
from mintpy.objects import sensor
# load all modules in this sub-directory for easy import
from mintpy.simulation.decorrelation import *
from mintpy.simulation.defo_model import *
from mintpy.simulation.fractal import *
from mintpy.utils import network as pnet, ptime, time_func

SPEED_OF_LIGHT = 299792458   # m/s



############################ Deformation Time-series ############################
def velocity2timeseries(date_list, vel=0.03, display=False):
    '''Simulate displacement time-series from linear velocity
    Inputs:
        date_list - list of string in YYYYMMDD or YYMMDD format
        vel        - float, velocity in meter per year
        display    - bool, display simulation or not
    Output:
        ts         - 2D np.array in size of (date_num,1), displacement time-series in m
    Example:
        date_list = pnet.read_baseline_file('bl_list.txt')[0]
        ts0 = simulate_timeseries(date_list, vel=0.03, display=True)
    '''
    tbase_list = ptime.date_list2tbase(date_list)[0]
    ts = vel / 365.25 * np.array(tbase_list)
    ts = ts.reshape(-1,1)

    if display:
        dates = ptime.date_list2vector(date_list)[0]
        ## Display
        marker_size = 5
        plt.figure()
        plt.scatter(dates, ts*100.0, s=marker_size**2)
        plt.xlabel('time [years]')
        plt.ylabel('LOS displacement [cm]')
        plt.title('displacement time-series with velocity = '+str(vel)+' m/yr')
        plt.show()
    return ts


def sim_variable_timeseries(num_date=100, display=False):
    """Simulate time variable displacement time-series

    Parameters: num_date - int, number of acquisitions
    Returns:    ts_sim   - 1D np.ndarray in size of (num_date,), displacement time-series in meters
    """
    tbase = np.linspace(0, num_date-1, num=num_date) * 12 # days

    ts_sim = np.zeros(num_date, np.float32)

    # comp 1 - exponential increase
    idx1 = int(0.2 * num_date)
    ts_sim[idx1:] = 0.01 * np.log(tbase[idx1:] - tbase[idx1-1])

    # comp 2 - step decrease + linear increase of 5 cm/yr
    idx2 = int(0.5 * num_date)
    ts_sim[idx2:] = 0.03 + 0.05 * (tbase[idx2:] - tbase[idx2-1]) / 365.25

    idx3 = int(0.8 * num_date)
    ts_sim[idx3:] = 0.

    # comp 3 - overall linear descrease of 0.5 cm/yr
    ts_sim += tbase * -0.005 / 365.25

    # scale and reference to the 1st image
    ts_sim *= 3.
    ts_sim -= ts_sim[0]

    if display:
        _, ax = plt.subplots(figsize=[6, 3])
        ax.plot(tbase, ts_sim * 100., '--')
        ax.set_xlabel('time [days]', fontsize=font_size)
        ax.set_ylabel('displacement [cm]', fontsize=font_size)
        ax.tick_params(direction='in', labelsize=font_size)
        plt.show()
    return ts_sim


def sim_variable_timeseries_v1(tbase, scale=3., display=False):
    """Simulate time variable displacement time-series

    Parameters: tbase  - 1D np.ndarray in size of (num_date,), temporal baseline in days
    Returns:    ts_sim - 1D np.ndarray in size of (num_date,), displacement time-series in meters
    """

    # Opt 2 - Time variable
    num_date = len(tbase)
    ts_sim = np.zeros(num_date, np.float32)

    # comp 1 - exponential increase
    idx1 = int(0.2*num_date)
    ts_sim[idx1:] = 0.01 * np.log(tbase[idx1:] - tbase[idx1-1])

    # comp 2 - step decrease + linear increase of 5 cm/yr
    idx2 = int(0.5*num_date)
    ts_sim[idx2:] = 0.03 + 0.05 * (tbase[idx2:] - tbase[idx2-1]) / 365.25

    idx3 = int(0.8*num_date)
    ts_sim[idx3:] = 0.

    # comp 3 - overall linear descrease of 0.5 cm/yr
    ts_sim += tbase * -0.005 / 365.25

    # scale and reference to the 1st image
    ts_sim *= 3.
    ts_sim -= ts_sim[0]

    if display:
        fig, ax = plt.subplots(figsize=[6, 3])
        ax.plot(tbase, ts_sim * 100., '--')
        ax.set_xlabel('time [days]', fontsize=font_size)
        ax.set_ylabel('displacement [cm]', fontsize=font_size)
        ax.tick_params(direction='in', labelsize=font_size)
        plt.show()

    return ts_sim


def timeseries2ifgram(ts_sim, date_list, date12_list, wvl=0.055, display=False, ax=None):
    """re-construct a stack of interferometric phase from a time-series
    Parameters: ts_sim      - 1D np.ndarray in size of (num_date,) in unit of m
                date_list   - list of str in YYYYMMDD in size of (num_date,)
                date12_list - list of str in YYYYMMDD_YYYYMMDD in size of (num_ifg,)
                wvl         - float, wavelength in m
                display     - bool, display the reconstruction result
                ax          - matplotlib.axes object, for display
    Returns:    ifgram_sim  - 1D np.ndarray in size of (num_ifg) in unit of radian
    """
    range2phase = -4.0 * np.pi / wvl
    num_ifgram = len(date12_list)
    ifgram_sim = np.zeros((num_ifgram,1), np.float32)
    for i in range(num_ifgram):
        m_date, s_date = date12_list[i].split('_')
        m_idx = date_list.index(m_date)
        s_idx = date_list.index(s_date)
        ifgram_sim[i] = ts_sim[s_idx] - ts_sim[m_idx]
    ifgram_sim *= range2phase

    if display:
        if not ax:
            fig, ax = plt.subplots(nrows=1, ncols=1)
        else:
            fig = None
        ifgram_sim_mat = pnet.coherence_matrix(date12_list, ifgram_sim)
        im = ax.imshow(ifgram_sim_mat, cmap='jet', interpolation='nearest')
        ax.set_xlabel('image number')
        ax.set_ylabel('image number')
        ax.set_title('interferometric phase')
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('phase [rad]')
        if fig is not None:
            plt.show()
    return ifgram_sim


def simulate_network(ts_sim, date12_list, decor_day, coh_resid, L=75, num_repeat=int(1e4),
                     baseline_file='bl_list.txt', sensor_name='Sen', inc_angle=33.4):
    """Simulate the InSAR stack for one pixel, including:
        simulated coherence --> decorrelation noise
        simulated ifgram phase with / without decorrelation noise
        estimated coherence"""
    # simulated (true) phase
    m_dates = [i.split('_')[0] for i in date12_list]
    s_dates = [i.split('_')[1] for i in date12_list]
    date_list = sorted(list(set(m_dates + s_dates)))
    wvl = SPEED_OF_LIGHT / sensor.SENSOR_DICT[sensor_name.lower()]['carrier_frequency']
    ifgram_sim = timeseries2ifgram(ts_sim, date_list, date12_list, wvl=wvl, display=False)

    # simulated (true) coherence
    coh_sim = pnet.simulate_coherence(date12_list,
                                      baseline_file=baseline_file,
                                      sensor_name=sensor_name,
                                      inc_angle=inc_angle,
                                      decor_time=decor_day,
                                      coh_resid=coh_resid)

    # simulated (estimated) phase
    decor_noise = coherence2decorrelation_phase(coh_sim, L=int(L), num_repeat=num_repeat)
    ifgram_est = decor_noise + np.tile(ifgram_sim.reshape(-1,1), (1, num_repeat))

    # estimated coherence
    coh_est = estimate_coherence(ifgram_est, L=L, win_size=25)

    return ifgram_est, coh_est, ifgram_sim, coh_sim


def estimate_coherence(ifgram, L=20, win_size=25):
    """Estimate coherence based on phase variance
    Reference:
      Rodriguez and Martin, 1992;
      Agram and Simons, 2015.
    Parameters: phase    : 2D np.array in size of (num_ifgram, num_repeat)
                L        : int, number of looks used to determine the phase PDF
                win_size : int, number of samples used to estimate phase variance
    Returns:    coh_est : 1D np.array in size of (num_ifgram,)
    """
    idx = np.random.choice(ifgram.shape[1], size=win_size)
    ifgram_diff = ifgram[:,idx]
    ifgram_diff -= np.tile(np.mean(ifgram_diff, axis=1).reshape(-1, 1), (1, win_size))
    ifgram_std = np.std(ifgram_diff, axis=1)
    coh_est = 1. / np.sqrt(1. + 2. * L * ifgram_std**2)
    return coh_est


def timeseries2velocity(date_list, defo_list):
    # date_list --> design_matrix
    A = time_func.get_design_matrix4time_func(date_list)
    A_inv = np.linalg.pinv(A)

    # least square inversion
    defo = np.array(defo_list, np.float32).reshape(-1,1)
    vel = np.dot(A_inv, defo)[1, :]
    return vel


def check_board(water_mask, grid_step=100, scale=1., display=True):
    length, width = water_mask.shape
    row_mask = np.ones((length, width), np.float32)
    for i in range(np.ceil(length/grid_step).astype(int)):
        i0 = i * grid_step
        i1 = min(length, (i + 1) * grid_step)
        if i % 2 == 0:
            row_mask[i0:i1, :] = -1

    col_mask = np.ones((length, width), np.float32)
    for i in range(np.ceil(width/grid_step).astype(int)):
        i0 = i * grid_step
        i1 = min(width, (i + 1) * grid_step)
        if i % 2 == 0:
            col_mask[:, i0:i1] = -1

    mask = np.multiply(row_mask, col_mask)
    mask[mask == -1] = 0.

    mask[water_mask == 0.] = np.nan
    mask *= scale

    if display:
        fig, ax = plt.subplots()
        im = ax.imshow(mask)
        fig.colorbar(im)
        plt.show()

    return mask


def add_unw_err2ifgram(ifgram, percentage=0.1, Nmax=2, print_msg=True):
    """Add unwrapping error to interferometric phase
    Parameters: ifgram     : 1D / 2D np.array in size of (num_ifgram, num_pixel) in float32
                percentage : float in [0, 1], percentage of interferograms with unwrapping errors
                Nmax       : int, maximum integer numbers of cycles of the added unwrapping errors
    Returns:    ifgram_err : 1D / 2D np.array in size of (num_ifgram, num_pixel) in float32
                idx_ifg_err : list of index, indicating interferograms with unwrapping errors
    """
    Nlist = np.hstack((np.arange(Nmax)+1, -1*np.arange(Nmax)-1))
    num_ifg_err = int(len(ifgram) * percentage)
    idx_ifg_err = random.sample(list(range(len(ifgram))), num_ifg_err)
    if print_msg:
        print(f'ifgram with unwrap error: {percentage}')
        print('unwrap error jump in 2*pi*(-{n}, {n}): '.format(n=Nmax))
        print(f'number of ifgrams with unwrap error: {num_ifg_err}')
    ifgram_err = np.array(ifgram, dtype=np.float32)
    ifgram_err[idx_ifg_err] += 2.*np.pi*np.random.choice(Nlist, size=num_ifg_err)
    return ifgram_err, idx_ifg_err
