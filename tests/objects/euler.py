#!/usr/bin/env python3
# Author: Yuan-Kai Liu, Oct 2022
"""Test mintpy.objects.euler module for the Euler pole and velocity computation."""

############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Yuan-Kai Liu, Oct 2022                           #
############################################################

import numpy as np
import pandas as pd

from mintpy.objects.euler import EulerPole
from mintpy.plate_motion import ITRF2014_PMM

# Read the Poles
pmm = pd.DataFrame.from_dict(data=ITRF2014_PMM, orient='index')
print('ITRF2014 No-net-rotation Plate Motion Model (11 plates)')
print(pmm)

# Inititalized the Euler poles for the plates
aust = EulerPole([ITRF2014_PMM['AUST'].omega_x, ITRF2014_PMM['AUST'].omega_y, ITRF2014_PMM['AUST'].omega_z])
eura = EulerPole([ITRF2014_PMM['EURA'].omega_x, ITRF2014_PMM['EURA'].omega_y, ITRF2014_PMM['EURA'].omega_z])
arab = EulerPole([ITRF2014_PMM['ARAB'].omega_x, ITRF2014_PMM['ARAB'].omega_y, ITRF2014_PMM['ARAB'].omega_z])

# Test printing the Euler pole information
print('\nAustralian plate')
aust.print_msg()
print('\nEurasian plate')
eura.print_msg()
print('\nArabian plate')
arab.print_msg()


# Results from UNAVCO calculator on the Eurasian plate:
# (https://www.unavco.org/software/geodetic-utilities/plate-motion-calculator/plate-motion-calculator.html)
truths = {
    'Latitude'       : [47, 43, 39, 35, 31],
    'Longitude'      : [109, 105, 101, 97, 93],
    'Speed [mm/yr]'  : [28.07, 28.57, 28.89, 29.02, 28.98],
    'Azimuth [deg]'  : [106.15, 103.70, 101.37, 99.11, 96.88],
    'E vel [mm/yr]'  : [26.96, 27.76, 28.32, 28.66, 28.77],
    'N vel [mm/yr]'  : [-7.81, -6.77, -5.70, -4.59, -3.47]}

TruthDF = pd.DataFrame.from_dict(truths)
print('\nGround truth from UNAVCO calculator')
print('Model: ITRF2014    Reference: NNR no-net-rotation')
print(TruthDF)
print('\n\n')



## Test our MintPy calculation from mintpy/objects/euler.py
lats = truths['Latitude']
lons = truths['Longitude']

v = 1e3 * np.array(eura.velocity_enu(lats, lons))
azi, speed = eura.velocity(lats, lons)

print('\nMintPy results computed by euler.py')
for i, (lat, lon) in enumerate(zip(lats, lons)):
    print(f'(lat,lon)=({lat:.2f}, {lon:.2f}) \t E={v[0,i]:.4f} \t N={v[1,i]:.4f} \t U={v[2,i]:.4f} \t \
    Speed={speed[i]*1e3:.2f} mm/yr \t Azi={azi[i]:.2f} deg')
