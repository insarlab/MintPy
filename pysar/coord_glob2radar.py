#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2017, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################


import os
import sys

import h5py
import numpy as np

import _readfile as readfile
import _pysar_utilities as ut


def usage():
    print('''usage: coord_glob2radar.py lat lon [trans_file] [hdf5_file_radarCoord]

Generates the sum of two input files.

example:
  coord_glob2radar.py 33.5 130.8
  coord_glob2radar.py 33.5 130.8 geomap_4rlks.trans velocity.h5
    ''')

################################################################################

def main(argv):
    if len(sys.argv) < 3:
        usage(); sys.exit(1)

    lat = float(argv[0])
    lon = float(argv[1])

    try:    trans_file = argv[2]
    except: trans_file = ut.get_lookup_file()

    try:    radar_file = argv[3]
    except: radar_file = 'unwrapIfgram.h5'

    atr_rdr = readfile.read_attribute(radar_file)
    
    print('input geo coord: lat=%.4f, lon=%.4f' % (lat, lon))
     
    y, x = ut.glob2radar(np.array(lat), np.array(lon), trans_file, atr_rdr)[0:2]
    print('corresponding radar coord: y=%d, x=%d' % (y, x))

    return

################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])  

