#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, Zhang Yunjun, 2013               #
############################################################


import os
import sys
import argparse

try:
    from skimage import filters, feature
except ImportError:
    raise ImportError('Could not import skimage!')

import numpy as np
from scipy import ndimage
from mintpy.utils import readfile, writefile


################################################################################################
EXAMPLE = """example:
  spatial_filter.py  velocity.h5
  spatial_filter.py  timeseries.h5  lowpass_avg        5
  spatial_filter.py  velocity.h5    lowpass_avg        5
  spatial_filter.py  velocity.h5    highpass_gaussian  3
  spatial_filter.py  velocity.h5    sobel
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Spatial filtering of 2D image.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='File to be filtered')
    parser.add_argument('filter_type', nargs='?', default='lowpass_gaussian',
                        choices=['lowpass_gaussian', 'highpass_gaussian',
                                 'lowpass_avg', 'highpass_avg',
                                 'sobel', 'roberts', 'canny'],
                        help='Type of filter. Default: lowpass_gaussian.\n' +
                             'For more filters, check the link below:\n' +
                             'http://scikit-image.org/docs/dev/api/skimage.filters.html')
    parser.add_argument('filter_par', nargs='?', type=float,
                        help='Filter parameter for low/high pass filter. Default=\n' +
                             'Sigma       for low/high pass gaussian filter, default: 3.0\n' +
                             'Kernel Size for low/high pass average filter, default: 5')
    parser.add_argument('-o', '--outfile',default=None, help='Output file name.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    inps.filter_type = inps.filter_type.lower()
    return inps


################################################################################################
def filter_data(data, filter_type, filter_par=None):
    """Filter 2D matrix with selected filter
    Inputs:
        data        : 2D np.array, matrix to be filtered
        filter_type : string, filter type
        filter_par  : string, optional, parameter for low/high pass filter
                      for low/highpass_avg, it's kernel size in int
                      for low/highpass_gaussain, it's sigma in float
    Output:
        data_filt   : 2D np.array, matrix after filtering.
    """

    if filter_type == "sobel":
        data_filt = filters.sobel(data)
    elif filter_type == "roberts":
        data_filt = filters.roberts(data)
    elif filter_type == "canny":
        data_filt = feature.canny(data)

    elif filter_type == "lowpass_avg":
        p = int(filter_par)
        kernel = np.ones((p, p), np.float32)/(p*p)
        data_filt = ndimage.convolve(data, kernel)
    elif filter_type == "highpass_avg":
        p = int(filter_par)
        kernel = np.ones((p, p), np.float32)/(p*p)
        lp_data = ndimage.convolve(data, kernel)
        data_filt = data - lp_data

    elif filter_type == "lowpass_gaussian":
        data_filt = filters.gaussian(data, sigma=filter_par)
    elif filter_type == "highpass_gaussian":
        lp_data = filters.gaussian(data, sigma=filter_par)
        data_filt = data - lp_data

    else:
        raise Exception('Un-recognized filter type: '+filter_type)

    return data_filt


############################################################
def filter_file(fname, filter_type, filter_par=None, fname_out=None):
    """Filter 2D matrix with selected filter
    Inputs:
        fname       : string, name/path of file to be filtered
        filter_type : string, filter type
        filter_par  : string, optional, parameter for low/high pass filter
                      for low/highpass_avg, it's kernel size in int
                      for low/highpass_gaussain, it's sigma in float
    Output:
        fname_out   : string, optional, output file name/path
    """
    # Info
    filter_type = filter_type.lower()
    atr = readfile.read_attribute(fname)
    k = atr['FILE_TYPE']
    msg = 'filtering {} file: {} using {} filter'.format(k, fname, filter_type)
    if filter_type.endswith('avg'):
        if not filter_par:
            filter_par = 5
        msg += ' with kernel size of {}'.format(filter_par)
    elif filter_type.endswith('gaussian'):
        if not filter_par:
            filter_par = 3.0
        msg += ' with sigma of {:.1f}'.format(filter_par)
    print(msg)

    # output filename
    if not fname_out:
        fname_out = '{}_{}{}'.format(os.path.splitext(fname)[0], filter_type,
                                     os.path.splitext(fname)[1])

    # filtering file
    dsNames = readfile.get_dataset_list(fname)
    maxDigit = max([len(i) for i in dsNames])
    dsDict = dict()
    for dsName in dsNames:
        msg = 'filtering {d:<{w}} from {f} '.format(
            d=dsName, w=maxDigit, f=os.path.basename(fname))
        data = readfile.read(fname, datasetName=dsName, print_msg=False)[0]
        if len(data.shape) == 3:
            num_loop = data.shape[0]
            for i in range(num_loop):
                data[i, :, :] = filter_data(data[i, :, :], filter_type, filter_par)
                sys.stdout.write('\r{} {}/{} ...'.format(msg, i+1, num_loop))
                sys.stdout.flush()
            print('')
        else:
            data = filter_data(data, filter_type, filter_par)
        dsDict[dsName] = data
    writefile.write(dsDict, out_file=fname_out, metadata=atr, ref_file=fname)
    return fname_out


################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    inps.outfile = filter_file(inps.file,
                               filter_type=inps.filter_type,
                               filter_par=inps.filter_par,
                               fname_out=inps.outfile)
    print('Done.')
    return inps.outfile


################################################################################################
if __name__ == '__main__':
    main()
