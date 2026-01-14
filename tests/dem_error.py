#!/usr/bin/env python3
# Author: Zhang Yunjun, Nov 2022
"""Test topographic residuel correction due to DEM errors."""


import argparse
import datetime
import sys

import numpy as np
from matplotlib import pyplot as plt

from mintpy.dem_error import estimate_dem_error
from mintpy.utils import ptime, time_func

################################################################################
# truth
delta_z_sim = 50   # meters

# setting: SAR geometry
range_dist = 800e3  # meters
inc_angle = 34      # degree
max_pbase = 500     # meters

# setting: time / acquisition
revisit_time = datetime.timedelta(days=24)
start_date = datetime.datetime(2018, 1, 1)
num_date = 50
ref_ind = 10

# setting: dem_error.py
phase_velocity = False
cond = 1e-8


################################################################################
EXAMPLE = """example:
  $MINTPY_HOME/tests/dem_error.py
  $MINTPY_HOME/tests/dem_error.py --plot
"""

def cmd_line_parse(iargs=None):
    # create parser
    parser = argparse.ArgumentParser(description='Test dem_error.py',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)
    parser.add_argument('--plot', dest='plot', action='store_true', help='Plot testing results.')

    # parsing
    inps = parser.parse_args(args=iargs)

    return inps


################################################################################
def sim_pbase_and_topo_residual(num_date, delta_z, ref_ind=0,
                                max_pbase=500, inc_angle=34, range_dist=800e3):

    # sim perp baseline TS
    np.random.seed(12138)
    pbase = np.random.rand(num_date) * max_pbase
    pbase -= pbase[ref_ind]

    # sim topo residual phase
    G_geom = pbase / (range_dist * np.sin(np.deg2rad(inc_angle)))
    ts_topo_res = np.dot(G_geom, delta_z_sim)

    return pbase, ts_topo_res


def plot_result(date_list, ts_sim, ts_obs=None, ts_cor=None, model=None):

    dt_list = ptime.date_list2vector(date_list)[0]

    _, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 4))
    ax.plot(dt_list, ts_sim*100, '--', color='C0', lw=3, label='simulation')
    if ts_obs is not None:
        ax.plot(dt_list, ts_obs*100, '.', color='C1', lw=3, label='observation')
    if ts_cor is not None:
        ax.plot(dt_list, ts_cor*100, 'o',  color='C2', ms=12, fillstyle='none', label='corrected')

    # plot step function
    if model is not None and 'stepDate' in model.keys():
        for step_date in model['stepDate']:
            step_dt = datetime.datetime.strptime(step_date, '%Y%m%d')
            ax.axvline(step_dt, ls='--', color='k')

    # axis format
    ax.set_xlabel('Time')
    ax.set_ylabel('Displacement [cm]')
    ax.legend()
    plt.show()


################################################################################
def test_dem_error_with_linear_defo(date_list, tbase, rel_tol=0.05, plot=False):
    print('Test 1: simple time-series with linear displacement.')

    # setting
    vel_sim = 0.05  # m/yr
    model = {'polynomial' : 1}

    # simulate displacement time-series
    ts_sim = vel_sim * tbase
    ts_sim -= ts_sim[ref_ind]

    # add topo residual
    pbase, ts_topo_res = sim_pbase_and_topo_residual(
        num_date,
        delta_z_sim,
        ref_ind,
        max_pbase=max_pbase,
        inc_angle=inc_angle,
        range_dist=range_dist,
    )

    # observed time-series: displacement + topo res
    ts_obs = ts_sim + ts_topo_res

    # estimate DEM error
    G_defo = time_func.get_design_matrix4time_func(date_list, model)
    G_geom = np.reshape(pbase / (range_dist * np.sin(np.deg2rad(inc_angle))), (-1,1))
    G = np.hstack((G_geom, G_defo))

    delta_z_est, ts_cor = estimate_dem_error(
        ts_obs, G, tbase,
        phase_velocity=phase_velocity,
        cond=cond,
    )[:2]

    # plot
    if plot:
        plot_result(date_list, ts_sim, ts_obs, ts_cor, model)

    # validate
    # use np.isclose(), instead of math.isclose() to handle better handle scalar/array types
    print(f'Specified DEM error: {delta_z_sim:.2f} m')
    print(f'Estimated DEM error: {delta_z_est[0]:.2f} m')
    assert np.isclose(delta_z_sim, delta_z_est, rtol=rel_tol)
    print('Pass.')


def test_dem_error_with_complex_defo(date_list, tbase, rel_tol=0.05, plot=False):
    print('Test 2: complex time-series with highly non-linear displacement.')

    # setting
    model = {'polynomial' : 2, 'stepDate' : ['20190818', '20200812']}

    # simulate displacement time-series
    # run the following to re-generate:
    #   from mintpy.simulation import simulation as sim
    #   ts_sim = sim.sim_variable_timeseries(50)
    ts_sim = np.array([
         0.        , -0.00049281, -0.00098563, -0.00147844, -0.00197125,
        -0.00246407, -0.00295688, -0.00344969, -0.0039425 , -0.00443532,
         0.06961907,  0.08992067,  0.10159181,  0.10972946,  0.11593095,
         0.12090778,  0.12503949,  0.12855262,  0.1315933 ,  0.13426131,
         0.13662781,  0.13874532,  0.14065379,  0.14238422,  0.14396119,
         0.0826078 ,  0.08704312,  0.09147844,  0.09591376,  0.10034907,
         0.1047844 ,  0.10921971,  0.11365503,  0.11809035,  0.12252566,
         0.126961  ,  0.1313963 ,  0.13583162,  0.14026694,  0.14470226,
        -0.01971253, -0.02020534, -0.02069815, -0.02119097, -0.02168378,
        -0.02217659, -0.0226694 , -0.02316222, -0.02365503, -0.02414784,
    ], dtype=np.float32)

    # add topo residual
    pbase, ts_topo_res = sim_pbase_and_topo_residual(
        num_date,
        delta_z_sim,
        ref_ind,
        max_pbase=max_pbase,
        inc_angle=inc_angle,
        range_dist=range_dist,
    )

    # observed time-series: displacement + topo res
    ts_obs = ts_sim + ts_topo_res

    # estimate DEM error
    G_defo = time_func.get_design_matrix4time_func(date_list, model)
    G_geom = np.reshape(pbase / (range_dist * np.sin(np.deg2rad(inc_angle))), (-1,1))
    G = np.hstack((G_geom, G_defo))

    delta_z_est, ts_cor = estimate_dem_error(
        ts_obs, G, tbase,
        phase_velocity=phase_velocity,
        cond=cond,
    )[:2]

    # plot
    if plot:
        plot_result(date_list, ts_sim, ts_obs, ts_cor, model)

    # validate
    print(f'Specified DEM error: {delta_z_sim:.2f} m')
    print(f'Estimated DEM error: {delta_z_est[0]:.2f} m')
    assert np.isclose(delta_z_sim, delta_z_est, rtol=rel_tol)
    print('Pass.')


################################################################################
def main(iargs=None):

    print('-'*50)
    print(f'Testing {__file__}')
    inps = cmd_line_parse(iargs)

    # prepare common data: date_list, tbase
    dt_list = [start_date + revisit_time * x for x in range(num_date)]
    date_list = [x.strftime('%Y%m%d') for x in dt_list]

    tbase = np.array(ptime.date_list2tbase(date_list)[0]) / 365.25   # years
    tbase -= tbase[ref_ind]

    # scenario 1 - simple linear deformation
    test_dem_error_with_linear_defo(
        date_list,
        tbase,
        rel_tol=0.05,
        plot=inps.plot,
    )

    # scenario 2 - complex nonlinear deformation (with more relaxed tolerance: 10%)
    test_dem_error_with_complex_defo(
        date_list,
        tbase,
        rel_tol=0.10,
        plot=inps.plot,
    )


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
