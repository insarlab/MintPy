#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2017, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################


import os
import sys

import h5py
import numpy as np

import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar._pysar_utilities as ut


def usage():
    print '''usage: coord_radar2glob.py az rg [trans_file] [hdf5_file_radarCoord]

Generates the sum of two input files.

example:
  coord_radar2glob.py 400 800
  coord_radar2glob.py 400 800 geomap_4rlks.trans velocity.h5
    '''

################################################################################

def main(argv):
    if len(sys.argv) < 3:
        usage(); sys.exit(1)

    y = int(argv[0])
    x = int(argv[1])

    try:    trans_file = argv[2]
    except: trans_file = ut.get_file_list('geomap*.trans')[0]

    try:    radar_file = argv[3]
    except: radar_file = 'unwrapIfgram.h5'

    atr_rdr = readfile.read_attribute(radar_file)
    
    print 'input radar coord: y/azimuth=%d, x/range=%d' % (y, x)
     
    lat, lon = ut.radar2glob(np.array(y), np.array(x), trans_file, atr_rdr)[0:2]
    print 'corresponding geo coord: lat=%.4f, lon=%.4f' % (lat, lon)

    return

################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])  

