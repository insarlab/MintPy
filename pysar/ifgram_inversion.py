#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################


import os
import sys
import argparse

import _readfile as readfile
import _pysar_utilities as ut


################################################################################################
EXAMPLE='''example:
  ifgram_inversion.py  unwrapIfgram.h5
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Inverse network of interferograms into timeseries.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('ifgram_file', help='interferograms file to be inversed')
    parser.add_argument('--method', dest='inverse_method', default='L2', choices=['L1','L2'],\
                        help='Inverse method, L1 or L2 norm minimization. Default: L2')
    parser.add_argument('-o','--output', dest='timeseries_file', default='timeseries.h5',\
                        help='output file name. Default: timeseries.h5')

    inps = parser.parse_args()
    return inps


################################################################################################
def main(argv):
    inps = cmdLineParse()

    # Input file info
    atr = readfile.read_attribute(inps.ifgram_file)
    k = atr['FILE_TYPE']
    if not k == 'interferograms':
        sys.exit('ERROR: only interferograms file supported, input is '+k+' file!')

    # Network Inversion
    if not inps.inverse_method == 'L1':
        print('Inverse time-series using L2 norm minimization')
        ut.timeseries_inversion(inps.ifgram_file, inps.timeseries_file)
    else:
        print('Inverse time-series using L1 norm minimization')
        ut.timeseries_inversion_L1(inps.ifgram_file, inps.timeseries_file)

    return inps.timeseries_file


################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])

