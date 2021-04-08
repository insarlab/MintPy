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
    from skimage import filters, feature, morphology
except ImportError:
    raise ImportError('Could not import skimage!')

import numpy as np
from scipy import ndimage
from mintpy.utils import readfile, writefile


################################################################################################
EXAMPLE = """example:
  spatial_filter.py  velocity.h5
  spatial_filter.py  timeseries.h5 -f lowpass_avg       -p 5
  spatial_filter.py  velocity.h5   -f lowpass_avg       -p 5
  spatial_filter.py  velocity.h5   -f highpass_gaussian -p 3
  spatial_filter.py  velocity.h5   -f sobel
  spatial_filter.py  ifgramStack.h5 unwrapPhase
  spatial_filter.py  ifgramStack.h5 unwrapPhase -f lowpass_avg -p 5
  spatial_filter.py  ifgramStack.h5 unwrapPhase -f double_difference -p 1 10
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Spatial filtering of 2D image.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='File to be filtered')
    parser.add_argument('dset', type=str, nargs='*', default=[],
                        help='optional - dataset(s) to filter (default: %(default)s).')
    parser.add_argument('-f', '--filter_type', dest='filter_type', nargs='?', default='lowpass_gaussian',
                        choices=['lowpass_gaussian', 'highpass_gaussian',
                                 'lowpass_avg', 'highpass_avg',
                                 'sobel', 'roberts', 'canny', 'double_difference'],
                        help='Type of filter. Default: lowpass_gaussian.\n' +
                             'For more filters, check the link below:\n' +
                             'http://scikit-image.org/docs/dev/api/skimage.filters.html')
    parser.add_argument('-p', '--filter_par', dest='filter_par', nargs='*', type=float,
                        help='Filter parameters for filters. Default=\n' +
                             'Sigma       for low/high pass gaussian filter, default: 3.0\n' +
                             'Kernel Size for low/high pass average filter, default: 5\n' +
                             'Kernel Radius for double difference local and regional filters, default: 1 10\n')
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
    
    elif filter_type == "double_difference":

        kernel = morphology.disk(filter_par[0], np.float32)
        kernel = kernel / kernel.flatten().sum()
        local_filt = ndimage.convolve(data, kernel)

        kernel = morphology.disk(filter_par[1], np.float32)
        kernel = kernel / kernel.flatten().sum()
        regional_filt = ndimage.convolve(data, kernel)

        data_filt = regional_filt - local_filt

    else:
        raise Exception('Un-recognized filter type: '+filter_type)

    return data_filt


############################################################
def filter_file(fname, ds_names=None, filter_type='lowpass_gaussian', filter_par=None, fname_out=None):
    """Filter 2D matrix with selected filter
    Inputs:
        fname       : string, name/path of file to be filtered
        ds_names    : list of string, datasets of interest
        filter_type : string, filter type
        filter_par  : string, optional, parameter for low/high pass filter
                      for low/highpass_avg, it's kernel size in int
                      for low/highpass_gaussain, it's sigma in float
                      for double_difference, it's local and regional kernel sizes in int
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
        else:
            filter_par = filter_par[0]
        msg += ' with kernel size of {}'.format(filter_par)
    elif filter_type.endswith('gaussian'):
        if not filter_par:
            filter_par = 3.0
        else:
            filter_par = filter_par[0]
        msg += ' with sigma of {:.1f}'.format(filter_par)
    elif filter_type == 'double_difference':
        if not filter_par:
            filter_par = [1,10]
        local, regional = filter_par
        msg += ' with local/regional kernel sizes of {}/{}'.format(local, regional)
    print(msg)

    # output filename
    if not fname_out:
        fname_out = '{}_{}{}'.format(os.path.splitext(fname)[0], filter_type,
                                     os.path.splitext(fname)[1])

    # filtering file
    ds_all = readfile.get_dataset_list(fname)
    if not ds_names:
        ds_names = ds_all
    ds_skips = set(ds_all)-set(ds_names)

    maxDigit = max([len(i) for i in ds_names])
    dsDict = dict()

    for ds_name in ds_skips:
        dsDict[ds_name] = readfile.read(fname, datasetName=ds_name, print_msg=False)[0]

    for ds_name in ds_names:
        msg = 'filtering {d:<{w}} from {f} '.format(d=ds_name, w=maxDigit, f=os.path.basename(fname))
        # read
        data = readfile.read(fname, datasetName=ds_name, print_msg=False)[0]
        # filter
        if len(data.shape) == 3:
            num_loop = data.shape[0]
            for i in range(num_loop):
                data[i, :, :] = filter_data(data[i, :, :], filter_type, filter_par)
                sys.stdout.write('\r{} {}/{} ...'.format(msg, i+1, num_loop))
                sys.stdout.flush()
            print('')
        else:
            data = filter_data(data, filter_type, filter_par)
        # write
        dsDict[ds_name] = data
    writefile.write(dsDict, out_file=fname_out, metadata=atr, ref_file=fname)
    return fname_out


################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    inps.outfile = filter_file(inps.file,
                               ds_names=inps.dset,
                               filter_type=inps.filter_type,
                               filter_par=inps.filter_par,
                               fname_out=inps.outfile)
    print('Done.')
    return inps.outfile


################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
