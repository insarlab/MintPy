#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Jun 2016: use read() and write() for IO
#                   add incidence_angle()
# Yunjun, Feb 2017: move incidence_angle() to _pysar_utilities


import sys
import os

import h5py
import numpy as np

import pysar.utils.readfile as readfile
import pysar.utils.writefile as writefile
import pysar.utils.utils as ut


def usage():
    print('''usage:  incidence_angle.py  file  [outfile]

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
    ''')
    return

def main(argv):
    try:
        File = argv[0]
        atr = readfile.read_attribute(File)
    except:
        usage();  sys.exit(1)
    
    try:    outFile = argv[1]
    except: outFile = 'incidenceAngle.h5'
    
    # Calculate look angle
    angle = ut.incidence_angle(atr, dimension=2)
    
    # Geo coord
    if 'Y_FIRST' in atr.keys():
        print('Input file is geocoded, only center incident angle is calculated: ')
        print(angle)
        length = int(atr['FILE_LENGTH'])
        width = int(atr['WIDTH'])
        angle_mat = np.zeros((length, width), np.float32)
        angle_mat[:] = angle
        angle = angle_mat

    print('writing >>> '+outFile)
    atr['FILE_TYPE'] = 'mask'
    atr['UNIT'] = 'degree'
    try: atr.pop('ref_date')
    except: pass
    writefile.write(angle, atr, outFile)
    return outFile

############################################################
if __name__ == '__main__':
    main(sys.argv[1:])





