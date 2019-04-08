#!/usr/bin/env python3
# Author: Emre Havazli

import os
import sys
import numpy as np
import h5py
from pysar.utils import readfile, writefile


################################################################################
USAGE = """
usage: correct_dem.py dem_file geo_demError_file

example:
  correct_dem.py $DEMDIR/Socorro-30/Socorro_30.dem geo_DEM_error.h5
  correct_dem.py $DEMDIR/Socorro-30/Socorro_30.dem geo_DEM_error.h5
"""


def usage():
    print(USAGE)


def main(argv):
    try:
        dem_file = argv[1]
        dem_error_file = argv[2]
    except:
        usage()
        sys.exit(1)

    print('Correcting the DEM')

    dem, demrsc = readfile.read(dem_file)
    dem_error = readfile.read(dem_error_file)

    dem_out = dem + dem_error
    writefile.write(dem_out, out_file='DEM_w_error.dem', metadata=demrsc)

    date12_file = open('111111-222222_baseline.rsc', 'w')
    date12_file.write('P_BASELINE_TOP_ODR'+'     '+'000')
    date12_file.close()
    return


################################################################################
if __name__ == '__main__':
    main(sys.argv[:])
