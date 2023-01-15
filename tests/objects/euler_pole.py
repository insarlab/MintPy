#!/usr/bin/env python3
# Author: Yuan-Kai Liu, Zhang Yunjun, Oct 2022
"""Test mintpy.objects.euler module for the Euler pole and velocity computation."""


import collections
import math

from mintpy.objects.euler_pole import ITRF2014_PMM, MASY2DMY, EulerPole

# validation against the UNAVCO Plate Motion Calculator
# https://www.unavco.org/software/geodetic-utilities/plate-motion-calculator/plate-motion-calculator.html
# accessed on Oct 29, 2022.
# config: height=0, model=ITRF2014, reference=NNR no-net-rotation
# note: azimuth here is measured from the north with positive for clock-wise direction
# unit:                              deg deg mm/yr   deg   mm/yr mm/yr
Tag = collections.namedtuple('Tag', 'lat lon speed azimuth vel_n vel_e')
POINT_PM = {
    'Australia' : Tag(-24,  132,  67.56,  28.94,  59.12,  32.69),
    'Eurasia'   : Tag( 27,   62,  28.85,  79.21,   5.40,  28.34),
    'Arabia'    : Tag( 18,   48,  46.51,  50.92,  29.32,  36.11),
}
PLATE_NAMES = POINT_PM.keys()


def test_euler_pole_initiation():
    print('Test 1: EulerPole object initiation and vector2pole conversion.')

    for plate_name in PLATE_NAMES:
        print(f'Plate name: ITRF2014-PMM {plate_name}')

        # get PMM info from Table 1 in Altamimi et al. (2017) as reference
        plate_pmm = ITRF2014_PMM[plate_name]

        # build EulerPole obj
        pole_obj = EulerPole(wx=plate_pmm.omega_x, wy=plate_pmm.omega_y, wz=plate_pmm.omega_z)

        # compare rotation rate: ITRF2014_PMM vs. Euler vector to pole conversion
        print(f'Reference  rotation rate from Altamimi et al. (2017): {plate_pmm.omega:.4f} deg/Ma')
        print(f'Calculated rotation rate from pole2vector conversion: {plate_pmm.omega:.4f} deg/Ma')
        assert math.isclose(plate_pmm.omega, pole_obj.rotRate*MASY2DMY, abs_tol=5e-4)
        print('Pass.')


def test_plate_motion_calc():
    print('Test 2: Plate motion calculation and validation against UNAVCO website.')

    for plate_name in PLATE_NAMES:
        # get UNAVCO result as reference
        point_pm = POINT_PM[plate_name]
        print(f'Plate = ITRF2014-PMM {plate_name}, point lat/lon = {point_pm.lat}/{point_pm.lon}')

        # calculate using EulerPole in m/yr
        plate_pmm = ITRF2014_PMM[plate_name]
        pole_obj = EulerPole(wx=plate_pmm.omega_x, wy=plate_pmm.omega_y, wz=plate_pmm.omega_z)
        ve, vn = pole_obj.get_velocity_enu(point_pm.lat, point_pm.lon, print_msg=False)[:2]

        print(f'Reference   (UNAVCO): vel_e={point_pm.vel_e:.2f}, vel_n={point_pm.vel_n:.2f} mm/yr')
        print(f'Calculation (MintPy): vel_e={ve*1e3:.2f}, vel_n={vn*1e3:.2f} mm/yr')
        assert math.isclose(point_pm.vel_e, ve*1e3, abs_tol=0.05)
        assert math.isclose(point_pm.vel_n, vn*1e3, abs_tol=0.05)
        print('Pass.')


if __name__ == '__main__':

    print('-'*50)
    print(f'Testing {__file__}')

    test_euler_pole_initiation()

    test_plate_motion_calc()
