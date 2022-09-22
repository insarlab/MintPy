#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, 2013                             #
############################################################


import sys

import numpy as np

import mintpy.utils.readfile as readfile

################################################################################
USAGE = """
usage: correlation_with_dem.py data_file dem_file [y0:y1 x0:x1]

Calculates the correlation of the DEM and file

example:
  correlation_with_dem.py velocity_masked.h5 inputs/geometryRadar.h5
  correlation_with_dem.py velocity_masked.h5 radar_8rlks.hgt
"""


def usage():
    print(USAGE)
    return


################################################################################
def main(argv):
    try:
        File = argv[0]
    except:
        usage()
        sys.exit(1)

    try:
        demFile = argv[1]
    except:
        demFile = 'inputs/geometryRadar.h5'

    dem = readfile.read(demFile, datasetName='height')[0]
    data, atr = readfile.read(File)
    print('input file is '+atr['FILE_TYPE'])

    # Subset
    try:
        y0, y1 = (int(i) for i in argv[2].split(':'))
        x0, x1 = (int(i) for i in argv[3].split(':'))
        data = data[y0:y1, x0:x1]
        dem = dem[y0:y1, x0:x1]
    except:
        pass

    # Calculation
    dem = dem.flatten(1)
    data = data.flatten(1)
    ndx = ~np.isnan(data)
    C1 = np.zeros([2, len(dem[ndx])])
    C1[0][:] = dem[ndx]
    C1[1][:] = data[ndx]

    # Display
    print('-------------------------------------------')
    print('Correlation with the DEM:  %.2f' % np.corrcoef(C1)[0][1])
    print('-------------------------------------------')
    print('DEM info:')
    print('    Max height difference: %.2f m' % (np.max(dem[ndx])-np.min(dem[ndx])))
    print('    Average        height: %.2f m' % np.mean(dem[ndx]))
    print('    Height            Std: %.2f m' % np.std(dem[ndx]))
    return


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
