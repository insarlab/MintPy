#!/usr/bin/env python3
# Author: Zhang Yunjun, Feb 2022
"""Test ascending/descending LOS decomposition into horizontal/vertical."""

import argparse
import os
import sys

import numpy as np
from matplotlib import pyplot as plt

from mintpy.asc_desc2horz_vert import asc_desc2horz_vert
from mintpy.utils import utils as ut

plt.rcParams.update({'font.size': 12})


################################################################################
EXAMPLE = """example:
  $MINTPY_HOME/tests/asc_desc2horz_vert.py
  $MINTPY_HOME/tests/asc_desc2horz_vert.py --plot
"""

def cmd_line_parse(iargs=None):
    # create parser
    parser = argparse.ArgumentParser(description='Test asc_desc2horz_vert.py',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)
    parser.add_argument('--plot', dest='plot', action='store_true', help='Plot testing results.')

    # parsing
    inps = parser.parse_args(args=iargs)

    return inps


################################################################################
def main(iargs=None):

    inps = cmd_line_parse(iargs)
    print('-'*50)
    print(os.path.abspath(__file__))

    ## Setup
    # LOS incidence / azimuth angles
    los_inc_angle = np.array([ 30,   30], dtype=np.float32)
    los_az_angle  = np.array([102, -102], dtype=np.float32)

    # specify horz / vert [truth]
    length, width = 5, 5
    simH = np.ones((length, width), dtype=np.float32) * 0.5
    simV = np.ones((length, width), dtype=np.float32) * 1.0
    # azimuth angle in horizontal direction in degrees
    # measured from the north with anti-clockwise as positive
    # [0 for north, -90 for east]
    horz_az_angle = 30


    ## Testing
    # horz / vert --> east / north / up
    dE = simH * np.sin(np.deg2rad(horz_az_angle)) * -1
    dN = simH * np.cos(np.deg2rad(horz_az_angle))
    dU = simV

    # east / north / up --> asc / desc LOS
    dlos0 = ut.enu2los(dE, dN, dU, inc_angle=los_inc_angle[0], az_angle=los_az_angle[0])
    dlos1 = ut.enu2los(dE, dN, dU, inc_angle=los_inc_angle[1], az_angle=los_az_angle[1])

    # asc / desc LOS --> horz / vert [estimation]
    dlos = np.vstack((dlos0.reshape(1, length, width), dlos1.reshape(1, length, width)))
    estH, estV = asc_desc2horz_vert(dlos, los_inc_angle, los_az_angle, horz_az_angle)

    # check difference between the trush and estimation
    print(f'mean difference for horz / vert: {np.nanmean(estH - simH)} / {np.nanmean(estH - simH)}')
    assert np.allclose(simH, estH)
    assert np.allclose(simV, estV)


    # Plotting
    if inps.plot:
        print(f'plot test result of {os.path.basename(__file__)}')
        fig, axs = plt.subplots(nrows=3, ncols=4, figsize=[8, 6], sharex=True, sharey=True)
        kwargs = dict(vmin=-1.5, vmax=1.5, cmap='RdBu', interpolation='nearest')

        # horz / vert [truth]
        ax = axs[0, 0];  im = ax.imshow(simH,  **kwargs);  ax.set_title('Horz [sim]')
        ax = axs[1, 0];  im = ax.imshow(simV,  **kwargs);  ax.set_title('Vert [sim]')
        ax = axs[2, 0];  ax.axis('off')
        # east / north / up
        ax = axs[0, 1];  im = ax.imshow(dE,    **kwargs);  ax.set_title('East [sim]')
        ax = axs[1, 1];  im = ax.imshow(dN,    **kwargs);  ax.set_title('North [sim]')
        ax = axs[2, 1];  im = ax.imshow(dU,    **kwargs);  ax.set_title('Up [sim]')
        # asc / desc
        ax = axs[0, 2];  im = ax.imshow(dlos0, **kwargs);  ax.set_title('Asc [obs]')
        ax = axs[1, 2];  im = ax.imshow(dlos1, **kwargs);  ax.set_title('Desc [obs]')
        ax = axs[2, 2];  ax.axis('off')
        # horz / vert [estimation]
        ax = axs[0, 3];  im = ax.imshow(estH,  **kwargs);  ax.set_title('Horz [est]')
        ax = axs[1, 3];  im = ax.imshow(estV,  **kwargs);  ax.set_title('Vert [est]')
        ax = axs[2, 3];  ax.axis('off')

        # axis format
        fig.tight_layout()
        # colorbar
        cax = fig.add_axes([0.6, 0.2, 0.3, 0.02])
        fig.colorbar(im, cax=cax, orientation='horizontal')
        plt.show()


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
