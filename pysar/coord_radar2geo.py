#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2017-2018, Zhang Yunjun                     #
# Author:  Zhang Yunjun                                    #
############################################################


import sys
import numpy as np
from pysar.utils import readfile, utils as ut


USAGE = """
usage: coord_radar2geo.py az rg [trans_file] [hdf5_file_radarCoord]

Generates the sum of two input files.

example:
  coord_radar2geo.py 400 800
  coord_radar2geo.py 400 800 geomap_4rlks.trans velocity.h5
"""


def usage():
    print(USAGE)
    return


################################################################################
def main(argv):
    if len(sys.argv) < 3:
        usage()
        sys.exit(1)

    y = int(argv[0])
    x = int(argv[1])
    print('input radar coord: y/azimuth=%d, x/range=%d' % (y, x))

    try:
        trans_file = argv[2]
    except:
        trans_file = ut.get_lookup_file()

    try:
        radar_file = argv[3]
    except:
        radar_file = 'INPUTS/ifgramStack.h5'
    atr_rdr = readfile.read_attribute(radar_file)

    coord = ut.coordinate(atr_rdr, lookup_file=trans_file)
    lat, lon = coord.radar2geo(np.array(y), np.array(x))[0:2]
    print('corresponding geo coord: lat=%.4f, lon=%.4f' % (lat, lon))
    return


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
