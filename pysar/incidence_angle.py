#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013, Heresh Fattahi, Zhang Yunjun          #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import os
import sys
import h5py
import numpy as np
from pysar.utils import readfile, writefile, utils as ut


############################################################
USAGE = """
usage:  incidence_angle.py  file  [outfile]

Generates incidence angles (in Radar Coordinate) for each pixel
  with required attributes read from the h5 file

input arguments:
  file    : string, input file name/path
  outfile : string, output file name/path for 2D incidence angle 
            calculated from file in radar coord

example:
  incidence_angle.py  velocity.h5
  incidence_angle.py  timeseries.h5
  incidence_angle.py  temporal_coherence.h5
"""


def usage():
    print(USAGE)
    return


############################################################
def main(argv):
    try:
        File = argv[0]
        atr = readfile.read_attribute(File)
    except:
        usage()
        sys.exit(1)

    try:
        outFile = argv[1]
    except:
        outFile = 'incidenceAngle.h5'

    # Calculate look angle
    angle = ut.incidence_angle(atr, dimension=2)

    # Geo coord
    if 'Y_FIRST' in atr.keys():
        print('Input file is geocoded, only center incident angle is calculated: ')
        print(angle)
        length = int(atr['LENGTH'])
        width = int(atr['WIDTH'])
        angle_mat = np.zeros((length, width), np.float32)
        angle_mat[:] = angle
        angle = angle_mat

    atr['FILE_TYPE'] = 'mask'
    atr['UNIT'] = 'degree'
    if 'REF_DATE' in atr.keys():
        atr.pop('REF_DATE')
    writefile.write(angle, out_file=outFile, metadata=atr)
    return outFile


############################################################
if __name__ == '__main__':
    main(sys.argv[1:])
