#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2018, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################


import os
import argparse
from pysar.utils import readfile, writefile, utils as ut


############################################################
EXAMPLE = """example:
  load2hdf5.py SanAndreas.dem      -o demGeo.h5
  load2hdf5.py geomap_4rlks.trans  -o geometryGeo.h5
  load2hdf5.py filt_fine_101120-110220.unw  -o 20101120_20110220_unw.h5
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Load binary data file(s) into an HDF5 file.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='file to be loaded.')
    parser.add_argument('-o', '--output', dest='outfile',
                        help='output HDF5 file name')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


############################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    if not inps.outfile:
        inps.outfile = os.path.splitext(inps.file)[0]+'.h5'

    atr = readfile.read_attribute(inps.file)
    dsNames = readfile.get_dataset_list(inps.file)
    dsDict = {}
    for dsName in dsNames:
        data = readfile.read(inps.file, datasetName=dsName)[0]
        dsDict[dsName] = data
    writefile.write(dsDict, out_file=inps.outfile, metadata=atr)
    return inps.outfile


############################################################
if __name__ == '__main__':
    main()
