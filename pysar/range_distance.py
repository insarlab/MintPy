#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2017, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################


import sys
import os

import h5py
import numpy as np

import pysar.utils.readfile as readfile
import pysar.utils.writefile as writefile
import pysar.utils.utils as ut


def usage():
    print('''
usage:  range_distance.py  file  [outfile]

Generates range distance (in Radar Coordinate) for each pixel
  with required attributes read from the h5 file

input arguments:
  file    : string, input file name/path
  outfile : string, output file name/path for 2D incidence angle 
            calculated from file in radar coord

example:
  range_distance.py  velocity.h5
  range_distance.py  timeseries.h5
  range_distance.py  temporal_coherence.h5
    ''')
    return

def main(argv):
    try:
        File = argv[0]
        atr = readfile.read_attribute(File)
    except:
        usage();  sys.exit(1)
    
    try:    outFile = argv[1]
    except: outFile = 'rangeDistance.h5'
    
    # Calculate look angle
    range_dis = ut.range_distance(atr, dimension=2)
    
    # Geo coord
    if 'Y_FIRST' in list(atr.keys()):
        print('Input file is geocoded, only center range distance is calculated: ')
        print(range_dis)
        length = int(atr['FILE_LENGTH'])
        width = int(atr['WIDTH'])
        range_dis_mat = np.zeros((length, width), np.float32)
        range_dis_mat[:] = range_dis
        range_dis = range_dis_mat

    print('writing >>> '+outFile)
    atr['FILE_TYPE'] = 'mask'
    atr['UNIT'] = 'm'
    try: atr.pop('ref_date')
    except: pass
    writefile.write(range_dis, atr, outFile)
    return outFile

############################################################
if __name__ == '__main__':
    main(sys.argv[1:])





