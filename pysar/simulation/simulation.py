#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2018, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################
# Recommend usage:
#   import pysar.simulation as psim

import numpy as np
import random
import scipy.stats as stats
import matplotlib.pyplot as plt
from pysar.objects import timeseries
from pysar.utils import ptime, network as pnet, utils as ut
from pysar.simulation.forward_model import mogi
from pysar.simulation.plot import *
from pysar import ifgram_inversion as ifginv


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
        plt.xlabel('Time (years)')
        plt.ylabel('LOS Displacement (cm)')
        plt.title('Displacement time-series with velocity = '+str(vel)+' m/yr')
        plt.show()
    return ts


def sim_variable_timeseries(tbase, scale=3., display=False):
    # Opt 2 - Time variable
    ts_sim = np.zeros(tbase.shape, np.float32)
    idx1 = 20
    ts_sim[idx1:] = 0.01 * np.log(tbase[idx1:] - tbase[idx1-1])
    idx2 = 50
    ts_sim[idx2:] = 0.03 + 0.06 * (tbase[idx2:] - tbase[idx2-1]) / 365.25
    idx3 = 70
    ts_sim[idx3:] = 0.
    ts_sim += tbase * -0.005 / 365.25
    ts_sim *= 3.
    ts_sim -= ts_sim[0]

    if display:
        fig, ax = plt.subplots(figsize=[6, 3])
        ax.plot(tbase, ts_sim * 100., '--')
        ax.set_xlabel('Time (days)', fontsize=font_size)
        ax.set_ylabel('Displacement (cm)', fontsize=font_size)
        ax.tick_params(direction='in', labelsize=font_size)
        plt.show()
    return ts_sim


def timeseries2ifgram(ts_sim, date_list, date12_list, wvl=0.055, display=False):
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
        ifgram_sim_mat = pnet.coherence_matrix(date12_list, ifgram_sim)
        plt.figure()
        plt.imshow(ifgram_sim_mat, cmap='jet')
        plt.xlabel('Image number')
        plt.ylabel('Image number')
        cbar = plt.colorbar()
        cbar.set_label('Phase (radian)')
        plt.title('Interferometric Phase')
        plt.show()
    return ifgram_sim


def sample_decorrelation_phase(L, coherence, size=1, display=False, scale=1.0, font_size=12):
    '''Sample decorrelation phase noise with PDF determined by L and coherence
    Inputs:
        L         - int, multilook number
        coherence - float, spatial coherence
        size      - int, sample number
    Output:
        sample    - 1D np.array in size of (size,), sampled phase
    unw_n = sample_decorrelation_phase(L=1, coherence=0.7, size=100000, display=True)
    '''
    phiNum = 100
    phiMax = np.pi * float(scale)
    pdf = ifginv.phase_pdf_ds(int(L), coherence, phi_num=phiNum)[0].flatten()   #for PS: ifginv.phase_variance_ps()
    phi = np.linspace(-phiMax, phiMax, phiNum+1, endpoint=True)
    phi_dist = stats.rv_histogram((pdf, phi))
    #sample = np.nan
    #while sample is np.nan:
    sample = phi_dist.rvs(size=size)

    if display:
        #size = 10000
        fig, ax = plt.subplots(figsize=[5,3])
        ax.hist(sample, bins=50, density=True, label='Sample\nHistogram\n(norm)')
        ax.plot(phi, phi_dist.pdf(phi), label='PDF')
        ax.plot(phi, phi_dist.cdf(phi), label='CDF')
        ax.set_xlabel('Phase', fontsize=font_size)
        ax.set_ylabel('Probability', fontsize=font_size)
        ax.set_title(r'L = %d, $\gamma$ = %.1f, sample size = %d' % (L, coherence, size), fontsize=font_size)
        ax.set_xlim([-np.pi, np.pi])
        ax.set_xticks([-np.pi, 0, np.pi])
        ax.set_xticklabels([r'-$\pi$', '0', r'$\pi$'], fontsize=font_size)
        ax.tick_params(direction='in', labelsize=font_size)
        ax.legend(fontsize=font_size)
        plt.savefig('DecorNoiseSampling.jpg', bbox_inches='tight', dpi=600)
        plt.show()
    return sample


def simulate_decorrelation_noises(date12_list, cohs, L=20, size:int=1, display=False, scale=1.0):
    '''Simuate decorrelation phase noise for input interferometric pairs
    Inputs:
        date12_list - list of string in YYMMDD-YYMMDD format, indicating pairs configuration
        cohs - 2D np.array in size of (ifgram_num,1)
        L    - int, multilook number
    Output:
        decorNoises - 2D np.array in size of (ifgram_num, 1)
    Example:
        from pysar.utils import network as pnet
        date12_list_all = pnet.get_date12_list('ifgram_list_all.txt')
        cohs = pnet.simulate_coherence(date12_list_all, decor_time=1000, coh_resid=0.2, display=True, inc_angle=22.8)
        decorNoises = simulate_decorrelation_noises(date12_list_all, cohs, L=20, display=True)
    '''
    ifgram_num = len(cohs)
    decorNoises = np.zeros((ifgram_num, size), np.float32)
    for i in range(ifgram_num):
        decorNoises[i, :] = sample_decorrelation_phase(int(L), cohs[i], size=size, scale=scale)

    if display:
        decorNoisesMat = pnet.coherence_matrix(date12_list, decorNoises)
        plt.figure()
        #plt.imshow(decorNoisesMat, vmin=-np.pi, vmax=np.pi, cmap='jet')
        plt.imshow(decorNoisesMat, cmap='jet')
        plt.xlabel('Image number')
        plt.ylabel('Image number')
        cbar = plt.colorbar()
        cbar.set_label('Decorrelation Phase Noise (radian)')
        plt.title('Decorrelation Noise')
        plt.show()
    return decorNoises


def simulate_network(ts_sim, date12_list, decor_day, coh_resid, L=75, num_sample=int(1e4),
                     baseline_file='bl_list.txt', sensor_name='Sen', inc_angle=33.4):
    """Simulate coherence --> decorrelation noise --> ifgram phase and estimated coherence"""
    coh_sim = pnet.simulate_coherence(date12_list,
                                      baseline_file='bl_list.txt',
                                      sensor_name=sensor_name,
                                      inc_angle=inc_angle,
                                      decor_time=decor_day,
                                      coh_resid=coh_resid)

    decor_noise = simulate_decorrelation_noises(date12_list, coh_sim, L=int(L), size=num_sample)

    m_dates = [i.split('_')[0] for i in date12_list]
    s_dates = [i.split('_')[1] for i in date12_list]
    date_list = sorted(list(set(m_dates + s_dates)))
    ifgram_sim = timeseries2ifgram(ts_sim, date_list, date12_list, display=False)
    ifgram_est = decor_noise + np.tile(ifgram_sim.reshape(-1,1), (1, num_sample))
    coh_est = estimate_coherence(ifgram_est, L=L, win_size=25)
    return ifgram_est, coh_est, ifgram_sim, coh_sim


def estimate_coherence(ifgram, L=20, win_size=25):
    """Estimate coherence based on phase variance
    Reference:
      Rodriguez and Martin, 1992;
      Agram and Simons, 2015.
    Parameters: phase    : 2D np.array in size of (num_ifgram, num_sample)
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
    A = timeseries.get_design_matrix4average_velocity(date_list)
    A_inv = np.linalg.pinv(A)

    # least square inversion
    defo = np.array(defo_list, np.float32).reshape(-1,1)
    vel = np.dot(A_inv, defo)[0, :]
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
        plt.colorbar(im)
        plt.show()

    return mask


def mogi_deformation(shape, source_geom, resolution=60., scale=1., display=True):
    """Simulate 2D deformation caused by the overpress of a Mogi source underneath
    
    Parameters: shape: 2-tuple of int in (length, width) or 2D np.ndarray in size of (length, width) in np.bool_
                source_geom : 4-tuple of float, Mogi source geometry: East, North, Depth, Volomn change in SI unit.
    
    """
    if isinstance(shape, np.ndarray):
        mask = np.multiply(np.array(shape != 0), ~np.isnan(shape))
        shape = mask.shape
    else:
        mask = np.ones(shape, np.bool_)

    length, width = shape
    yy, xx = np.mgrid[0:length:length*1j, 0:width:width*1j]
    yy *= resolution
    xx *= resolution
    xloc = np.vstack((xx.reshape(1, -1), yy.reshape(1, -1)))

    dis_map = mogi(source_geom, xloc)[0]
    dis_e = dis_map[0, :].reshape(length, width)
    dis_n = dis_map[1, :].reshape(length, width)
    dis_u = dis_map[2, :].reshape(length, width)
    dis_los = ut.enu2los(dis_e, dis_n, dis_u)

    dis_los[mask == 0.] = np.nan
    dis_los *= scale

    if display:
        fig, ax = plt.subplots(1, 4, figsize=[10, 3])
        dmin = np.nanmin(dis_los)
        dmax = np.nanmax(dis_los)
        for i, fig_title in enumerate(['east','north','vertical']):
            ax[i].imshow(dis_map[i, :].reshape(length, width), vmin=dmin, vmax=dmax)
            ax[i].set_title(fig_title)
        im = ax[3].imshow(dis_los, vmin=dmin, vmax=dmax)
        ax[3].set_title('los - SenD')
        fig.subplots_adjust(right=0.90)
        cax = fig.add_axes([0.92, 0.25, 0.01, 0.5])
        cbar = fig.colorbar(im, cax=cax)
        cbar.set_label('Displacement [m]')
        plt.show()

    return dis_los


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
        print('ifgram with unwrap error: {}'.format(percentage))
        print('unwrap error jump in 2*pi*(-{n}, {n}): '.format(n=Nmax))
        print('number of ifgrams with unwrap error: {}'.format(num_ifg_err))
    ifgram_err = np.array(ifgram, dtype=np.float32)
    ifgram_err[idx_ifg_err] += 2.*np.pi*np.random.choice(Nlist, size=num_ifg_err)
    return ifgram_err, idx_ifg_err

