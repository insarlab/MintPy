#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2018, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################


import os
import argparse
import numpy as np
from pysar.utils import readfile, writefile


############################################################
EXAMPLE = """example:
  load2hdf5.py SanAndreas.dem  --data-type float32  -o demGeo.h5
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
    parser.add_argument('--data-type', dest='data_type', help='output data type')
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

    if inps.data_type:
        if inps.data_type in ['float', 'float32', 'np.float32']:
            inps.data_type = np.float32
        elif inps.data_type in ['float64', 'np.float64']:
            inps.data_type = np.float64
        elif inps.data_type in ['int', 'int16', 'np.int16']:
            inps.data_type = np.int16
        elif inps.data_type in ['bool', 'np.bool_']:
            inps.data_type = np.bool_
        elif inps.data_type in ['complex', 'np.complex64']:
            inps.data_type = np.complex64
        elif inps.data_type in ['complex128', 'np.complex128']:
            inps.data_type = np.complex128
        else:
            raise ValueError('un-recognized input data type: {}'.format(inps.data_type))

    atr = readfile.read_attribute(inps.file)
    dsNames = readfile.get_dataset_list(inps.file)
    dsDict = {}
    for dsName in dsNames:
        data = readfile.read(inps.file, datasetName=dsName)[0]
        if inps.data_type:
            data = np.array(data, inps.data_type)
        dsDict[dsName] = data
    writefile.write(dsDict, out_file=inps.outfile, metadata=atr)
    return inps.outfile


############################################################
if __name__ == '__main__':
    main()
