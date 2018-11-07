#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013-2018, Zhang Yunjun, Heresh Fattahi     #
# Author:  Zhang Yunjun, Heresh Fattahi                    #
############################################################


import os
import sys
import argparse
import numpy as np
from pysar.utils import readfile, writefile, utils as ut


############################################################
EXAMPLE = """example:
  mask.py  velocity.h5     -m Mask.h5
  mask.py  timeseries.h5   -m temporal_coherence.h5  -t 0.7
  mask.py  ifgramStack.h5  -m 100102_101120.cor      -t 0.9  -y  200 300  -x 300 400
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Mask file',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='File to be masked')
    parser.add_argument('-m', '--mask', dest='mask_file',
                        help='mask for pixels used in ramp estimation')
    parser.add_argument(
        '-o', '--outfile', help='Output file name. Disabled when more than 1 input files')

    # modify input mask
    parser.add_argument('-t', dest='threshold', type=float,
                        help='threshold value used for masking.\n' +
                        'if not specified, only pixels with mask value equal to zero is masked out.')
    parser.add_argument('--fill', dest='fill_value', type=float,
                        help="fill masked out area with input value. i.e. \n"
                             "np.nan, 0, 1000, ... \n"
                             "By default, it's np.ma.masked for int16 type and np.nan for all the others.")
    parser.add_argument('-x', dest='subset_x', type=int, nargs=2,
                        help='subset range in x/cross-track/column direction')
    parser.add_argument('-y', dest='subset_y', type=int, nargs=2,
                        help='subset range in y/along-track/row direction')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


############################################################
def mask_matrix(data, mask, fill_value=None):
    """mask a 2D matrxi data with mask
    Parameters: data : 2D / 3D np.array
                mask : 2D np.array of bool
                fill_value : number
    Returns:    data : same shape of array as input
    """
    # Masked Value
    if fill_value is None:
        if data.dtype == np.float32:
            fill_value = np.nan
        else:
            raise ValueError('specify fill_value for input data type: {}'.format(data.dtype))

    if len(data.shape) == 2:
        data[mask == 0] = fill_value
    elif len(data.shape) == 3:
        data[:, mask == 0] = fill_value
    return data


def update_mask_with_inps(mask, inps=None, print_msg=True):
    """Update mask matrix from input options: subset_x/y and threshold"""
    if not inps:
        inps = cmd_line_parse()

    if inps.subset_x:
        mask[:, 0:inps.subset_x[0]] = 0
        mask[:, inps.subset_x[1]:] = 0
        if print_msg:
            print('mask out area not in x: {}'.format(inps_dict['subset_x']))

    if inps.subset_y:
        mask[0:inps.subset_y[0], :] = 0
        mask[inps.subset_y[1]:, :] = 0
        if print_msg:
            print('mask out area not in y: {}'.format(inps_dict['subset_y']))

    if inps.threshold:
        mask[mask < inps.threshold] = 0
        if print_msg:
            print('mask out pixels with value < {} in mask file'.format(inps.threshold))
    return mask


def mask_file(fname, mask_file, out_file=None, inps=None):
    """ Mask input fname with mask_file
    Inputs:
        fname/mask_file - string, 
        inps_dict - dictionary including the following options:
                    subset_x/y - list of 2 ints, subset in x/y direction
                    threshold - float, threshold/minValue to generate mask
    Output:
        out_file - string
    """
    if not inps:
        inps = cmd_line_parse()

    if not out_file:
        out_file = '{}_masked{}'.format(os.path.splitext(fname)[0], os.path.splitext(fname)[1])

    # read mask_file
    mask = readfile.read(mask_file)[0]
    mask = update_mask_with_inps(mask, inps)

    # masking input file
    dsNames = readfile.get_dataset_list(fname)
    maxDigit = max([len(i) for i in dsNames])
    dsDict = {}
    for dsName in dsNames:
        if dsName not in ['coherence']:
            print('masking {d:<{w}} from {f} ...'.format(d=dsName, w=maxDigit, f=fname))
            data = readfile.read(fname, datasetName=dsName, print_msg=False)[0]
            data = mask_matrix(data, mask, fill_value=inps.fill_value)
        dsDict[dsName] = data
    writefile.write(dsDict, out_file=out_file, ref_file=fname)


############################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    inps.outfile = mask_file(inps.file, inps.mask_file, inps.outfile, inps)

    print('Done.')
    return


############################################################
if __name__ == '__main__':
    main()
