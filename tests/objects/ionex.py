#!/usr/bin/env python3
# Author: Zhang Yunjun, Jun 2022
"""Test mintpy.objects.ionex module for the IONEX file reading and interpolation."""


import os

import numpy as np

from mintpy.objects import ionex

# test dataset at $MINTPY_HOME/tests/data directory
tec_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
date_str = '20151115'
sol_code = 'jpl'


def prep_test_data(prep_mode=True):
    """Prepare the test data."""

    os.makedirs(tec_dir, exist_ok=True)
    tec_file = ionex.get_ionex_filename(
        date_str,
        tec_dir=tec_dir,
        sol_code=sol_code,
    )
    print(f'TEC file: {tec_file}')

    if prep_mode:
        # download
        if not os.path.isfile(tec_file):
            print(f'download IONEX at {date_str} from {sol_code} to {tec_file}')
            ionex.dload_ionex(date_str, tec_dir=tec_dir, sol_code=sol_code)

        # plot
        print('plotting IONEX file as animation')
        ionex.plot_ionex(tec_file, save_fig=False)

    return tec_file



def test_read_ionex():
    """Test mintpy.objects.ionex.read_ionex()"""

    print('Test 1: read GIM files via ionex.read_ionex()')
    time_ind = 1
    x0, x1, y0, y1 = 3, 9, 28, 33
    tec_aoi = np.array(
        [[71.8, 64.2, 55.9, 47.1, 38.6, 31.2],
         [80.2, 73.9, 66.6, 58.2, 49.5, 41.1],
         [83.2, 79.6, 74.6, 68. , 60.1, 51.6],
         [79.6, 79.5, 78.1, 74.5, 68.5, 60.9],
         [71.9, 74.5, 76.5, 76.2, 73.1, 67.3]],
    )

    # get IONEX file path
    tec_file = ionex.get_ionex_filename(
        date_str,
        tec_dir=tec_dir,
        sol_code=sol_code,
    )

    # read IONEX
    mins, lats, lons, tec_maps = ionex.read_ionex(tec_file)[:4]

    # compare
    assert np.allclose(tec_maps[time_ind, y0:y1, x0:x1], tec_aoi)
    print('Pass.')


def test_get_ionex_value():
    """Test mintpy.objects.ionex.get_ionex_value()"""

    print('Test 2: Get the ionex for a given location/time via ionex.get_ionex_value()')
    lat, lon = -21.3, -67.4         # northern Chile
    utc_sec = 23 * 3600 + 7 * 60    # 23:07 UTC

    methods = ['nearest', 'linear2d', 'linear3d', 'linear3d']
    rotates = [False, False, False, True]
    values = [60.8, 58.90687978, 64.96605174, 65.15525905]

    # get IONEX file path
    tec_file = ionex.get_ionex_filename(
        date_str,
        tec_dir=tec_dir,
        sol_code=sol_code,
    )

    # compare
    for method, rotate, value in zip(methods, rotates, values):
        tec_val = ionex.get_ionex_value(
            tec_file, utc_sec, lat, lon,
            interp_method=method,
            rotate_tec_map=rotate,
        )
        assert np.allclose(tec_val, value, atol=1e-05, rtol=1e-05)
        print(f'Interpolation method ({method:8s}) with rotation ({str(rotate):5s}): Pass.')


def test_plot_ionex():
    """Test mintpy.objects.ionex.plot_ionex()"""
    # get IONEX file path
    tec_file = ionex.get_ionex_filename(
        date_str,
        tec_dir=tec_dir,
        sol_code=sol_code,
    )

    # plot IONEX file
    ionex.plot_ionex(tec_file, save_fig=True)


if __name__ == '__main__':

    print('-'*50)
    print(f'Testing {__file__}')

    prep_test_data(prep_mode=False)

    test_read_ionex()

    test_get_ionex_value()

    #test_plot_ionex()
