#!/usr/bin/env python3
# Author: Yuan-Kai Liu, Oct 2022
"""Test mintpy.objects.euler module for the Euler pole and velocity computation."""


import numpy as np

from mintpy.objects.euler_pole import EulerPole
from mintpy.plate_motion import ITRF2014_PMM


def test_build_euler_pole():
    # Read the Poles
    print('ITRF2014 No-net-rotation Plate Motion Model (Altamimi et al., 2017)')

    Tags = ['Tag','name','num_site','omega_x','omega_y','omega_z','omega','wrms_e','wrms_n']
    print('{:8s}{:15s}{:10s}{:12s}{:12s}{:12s}{:10s}{:10s}{:10s}'.format(*Tags))
    for k, v in ITRF2014_PMM.items():
        print('{:8s}{:15s}{:3d}{:12.3f}{:12.3f}{:12.3f}{:12.3f}{:12.3f}{:12.3f}'.format(
            k, v.name, v.num_site, v.omega_x, v.omega_y, v.omega_z, v.omega, v.wrms_e, v.wrms_n))

    # Inititalized the Euler poles for the plates
    aust = EulerPole(wx=ITRF2014_PMM['AUST'].omega_x, wy=ITRF2014_PMM['AUST'].omega_y, wz=ITRF2014_PMM['AUST'].omega_z)
    eura = EulerPole(wx=ITRF2014_PMM['EURA'].omega_x, wy=ITRF2014_PMM['EURA'].omega_y, wz=ITRF2014_PMM['EURA'].omega_z)
    arab = EulerPole(wx=ITRF2014_PMM['ARAB'].omega_x, wy=ITRF2014_PMM['ARAB'].omega_y, wz=ITRF2014_PMM['ARAB'].omega_z)
    return aust, eura, arab


def read_validation_truth():
    # Results from UNAVCO calculator on the Eurasian plate:
    # (https://www.unavco.org/software/geodetic-utilities/plate-motion-calculator/plate-motion-calculator.html)
    truths = {
        'Latitude'       : [47, 43, 39, 35, 31],
        'Longitude'      : [109, 105, 101, 97, 93],
        #'Speed [mm/yr]'  : [28.07, 28.57, 28.89, 29.02, 28.98],
        #'Azimuth [deg]'  : [106.15, 103.70, 101.37, 99.11, 96.88],
        'E vel [mm/yr]'  : [26.96, 27.76, 28.32, 28.66, 28.77],
        'N vel [mm/yr]'  : [-7.81, -6.77, -5.70, -4.59, -3.47]}

    print('\nGround truth from UNAVCO calculator:')
    print('Model: ITRF2014 Eurasian    Reference: NNR no-net-rotation')
    tests_res = []
    for i in range(len(truths['Latitude'])):
        tests_res.append([])
        for v in truths.values():
            tests_res[i].append(v[i])
    print('\t'.join(truths.keys()))
    for i, res in enumerate(tests_res):
        print('\t\t'.join(np.array(res).astype('str')))
    print('\n\n')
    return truths


def test_validate_euler_pole():
    # Build plate models
    aust, eura, arab = test_build_euler_pole()
    print('\n\nExample poles:')
    print('\n < Australian plate >')
    aust.print_info()
    print('\n < Eurasian plate >')
    eura.print_info()
    print('\n < Arabian plate >')
    arab.print_info()

    # Read UNAVCO as ground truth
    truths = read_validation_truth()

    # Test our MintPy calculation from mintpy/objects/euler_pole.py
    lats = truths['Latitude']
    lons = truths['Longitude']

    print('Compute linear velocity from the Euler pole...')
    v = 1e3 * np.array(eura.get_velocity_enu(lats, lons))
    #azi, speed = eura.getVelocityAzi(lats, lons)

    print('\nComputed results by euler_pole.py:')
    print('\t'.join(truths.keys())+'\t U vel [mm/yr]')
    for i, (lat, lon) in enumerate(zip(lats, lons)):
        speed = np.sqrt(v[0,i]**2 + v[1,i]**2 + v[2,i]**2)
        print(f'{lat:.2f}\t\t{lon:.2f}\t\t{v[0,i]:.4f}\t\t{v[1,i]:.2f}\t\t{v[2,i]:.4f}')

    # compare
    assert np.allclose([v[0], v[1]], [truths['E vel [mm/yr]'], truths['N vel [mm/yr]']], rtol=1e-02)



if __name__ == '__main__':

    print(f'Testing {__file__}')

    test_validate_euler_pole()

    print('\nTesting is successful!')
