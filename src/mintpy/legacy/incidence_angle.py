#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, Zhang Yunjun, 2013               #
############################################################


import argparse

import numpy as np

from mintpy.utils import readfile, utils as ut, writefile

############################################################
EXAMPLE = """example:
  incidence_angle.py  velocity.h5
  incidence_angle.py  timeseries.h5         -d inputs/geometryRadar.h5
  incidence_angle.py  temporalCoherence.h5
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Generate incidence angle for each pixel',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='file containing STARTING_RANGE / RANGE_PIXEL_SIZE metadata')
    parser.add_argument('-d', '--dem', dest='dem_file', help='DEM file')
    parser.add_argument('-o', '--output', dest='outfile',
                        help='output file name/path for 2D incidence angle ')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


############################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    # Calculate look angle
    atr = readfile.read_attribute(inps.file)
    dem = None
    if inps.dem_file:
        dem = readfile.read(inps.dem_file, datasetName='height')[0]
    angle = ut.incidence_angle(atr, dem=dem, dimension=2)

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

    if not inps.outfile:
        inps.outfile = 'incidenceAngle.h5'
    writefile.write(angle, out_file=inps.outfile, metadata=atr)
    return inps.outfile


############################################################
if __name__ == '__main__':
    main()
