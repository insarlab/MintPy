#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import argparse
import sys

import numpy as np

from mintpy.objects import ifgramStack
from mintpy.utils import readfile, writefile

#####################################################################################
EXAMPLE = """example:
  ifgram_reconstruction.py timeseries_ERA5_ramp_demErr.h5
  ifgram_reconstruction.py timeseries_ERA5_ramp_demErr.h5 -r inputs/ifgramStack.h5 -o ifgramStackRecon.h5
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Reconstruct network of interferograms from time-series',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('timeseries_file', type=str, help='time-series file.')
    parser.add_argument('-r', dest='ifgram_file', type=str, default='./inputs/ifgramStack.h5',
                        help='reference interferograms stack file')
    parser.add_argument('-o','--output', dest='out_file', default='ifgramStackRecon.h5',
                        help='output filename for the reconstructed interferograms.')
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


#####################################################################################
def timeseries2ifgram(ts_file, ifgram_file, out_file='reconUnwrapIfgram.h5'):
    # read time-series
    atr = readfile.read_attribute(ts_file)
    range2phase = -4.*np.pi / float(atr['WAVELENGTH'])
    print(f'reading timeseries data from file {ts_file} ...')
    ts_data = readfile.read(ts_file)[0] * range2phase
    num_date, length, width = ts_data.shape
    ts_data = ts_data.reshape(num_date, -1)

    # reconstruct unwrapPhase
    print('reconstructing the interferograms from timeseries')
    stack_obj = ifgramStack(ifgram_file)
    stack_obj.open(print_msg=False)
    date12_list = stack_obj.get_date12_list(dropIfgram=False)
    A = stack_obj.get_design_matrix4timeseries(date12_list, refDate='no')[0]
    ifgram_est = np.dot(A, ts_data).reshape(A.shape[0], length, width)
    ifgram_est = np.array(ifgram_est, dtype=ts_data.dtype)
    del ts_data

    # write to ifgram file
    dsDict = {}
    dsDict['unwrapPhase'] = ifgram_est
    writefile.write(dsDict, out_file=out_file, ref_file=ifgram_file)
    return ifgram_file


def main(iargs=None):
    inps = cmd_line_parse(iargs)
    timeseries2ifgram(inps.timeseries_file, inps.ifgram_file, inps.out_file)
    return


#####################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
