#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os
import sys
import argparse
import numpy as np
from mintpy.utils import readfile, writefile


################################################################################
EXAMPLE = """example:
  add.py mask_1.h5 mask_2.h5 mask_3.h5 -o mask_all.h5
  add.py 081008_100220.unw 100220_110417.unw -o 081008_110417.unw
  add.py timeseries_ERA5.h5 inputs/ERA5.h5 -o timeseries.h5
"""


def create_parser():
    """ Command line parser """
    parser = argparse.ArgumentParser(description='Generate sum of multiple input files.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs='+', help='files (2 or more) to be added')
    parser.add_argument('-o', '--output', dest='outfile', help='output file name')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    if len(inps.file) < 2:
        parser.print_usage()
        sys.exit('ERROR: At least 2 input files needed!')
    return inps


################################################################################
def add_matrix(data1, data2):
    """Sum of 2 input matrix"""
    data1 = np.array(data1, dtype=np.float32)
    data2 = np.array(data2, dtype=np.float32)
    data = data1 + data2
    data[np.isnan(data1)] = data2[np.isnan(data1)]
    data[np.isnan(data2)] = data1[np.isnan(data2)]
    return data


def add_file(fnames, out_file=None):
    """Generate sum of all input files
    Parameters: fnames : list of str, path/name of input files to be added
                out_file : str, optional, path/name of output file
    Returns:    out_file : str, path/name of output file
    Example:    'mask_all.h5' = add_file(['mask_1.h5','mask_2.h5','mask_3.h5'], 'mask_all.h5')
    """
    # Default output file name
    ext = os.path.splitext(fnames[0])[1]
    if not out_file:
        out_file = os.path.splitext(fnames[0])[0]
        for i in range(1, len(fnames)):
            out_file += '_plus_' + os.path.splitext(os.path.basename(fnames[i]))[0]
        out_file += ext

    dsDict = {}
    dsNames = readfile.get_dataset_list(fnames[0])
    for dsName in dsNames:
        # ignore dsName if input file has single dataset
        if len(dsNames) == 1:
            dsName2read = None
        else:
            dsName2read = dsName

        print('adding {} ...'.format(dsName))
        data = readfile.read(fnames[0], datasetName=dsName2read)[0]
        for i in range(1, len(fnames)):
            data2 = readfile.read(fnames[i], datasetName=dsName2read)[0]
            data = add_matrix(data, data2)
        dsDict[dsName] = data

    # output
    atr = readfile.read_attribute(fnames[0])
    print('use metadata from the 1st file: {}'.format(fnames[0]))
    writefile.write(dsDict, out_file=out_file, metadata=atr, ref_file=fnames[0])
    return out_file


################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    print('input files to be added: ({})\n{}'.format(len(inps.file), inps.file))

    inps.outfile = add_file(inps.file, inps.outfile)

    print('Done.')
    return inps.outfile


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
