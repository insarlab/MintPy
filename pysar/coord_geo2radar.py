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
usage: coord_geo2radar.py lat lon [trans_file] [hdf5_file_radarCoord]

Generates the sum of two input files.

example:
  coord_geo2radar.py 33.5 130.8
  coord_geo2radar.py 33.5 130.8 INPUTS/geometryRadar.h5
  coord_geo2radar.py 33.5 130.8 geomap_4rlks.trans velocity.h5
"""


def usage():
    print(USAGE)
    return


################################################################################
def main(argv):
    if len(sys.argv) < 3:
        usage()
        sys.exit(1)

    lat = float(argv[0])
    lon = float(argv[1])

    try:
        trans_file = argv[2]
    except:
        trans_file = ut.get_lookup_file()

    try:
        radar_file = argv[3]
    except:
        radar_file = 'INPUTS/ifgramStack.h5'

    atr_rdr = readfile.read_attribute(radar_file)

    print('input geo coord: lat=%.4f, lon=%.4f' % (lat, lon))

    coord = ut.coordinate(atr_rdr, lookup_file=trans_file)
    y, x = coord.geo2radar(np.array(lat), np.array(lon))[0:2]
    print('corresponding radar coord: y=%d, x=%d' % (y, x))
    return


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
