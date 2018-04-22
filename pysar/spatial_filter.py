#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Jul 2017: rewrite using pysar module


import os
import sys
import argparse

try:
    from skimage import filters, feature
except:
    print('++++++++++++++++++++++++++++++++++++++++++++')
    print('Could not import skimage')
    print('To use filter.py you must install skimage')
    print('See: http://scikit-image.org/')
    print('++++++++++++++++++++++++++++++++++++++++++++')
    sys.exit(1)

import h5py
import numpy as np
#from PIL import Image
from scipy import ndimage

from pysar.utils import readfile, writefile, datetime as ptime
from pysar.utils.readfile import multi_group_hdf5_file, multi_dataset_hdf5_file, single_dataset_hdf5_file


################################################################################################
def filter_data(data, filter_type, filter_par=None):
    '''Filter 2D matrix with selected filter
    Inputs:
        data        : 2D np.array, matrix to be filtered
        filter_type : string, filter type
        filter_par  : string, optional, parameter for low/high pass filter
                      for low/highpass_avg, it's kernel size in int
                      for low/highpass_gaussain, it's sigma in float
    Output:
        data_filt   : 2D np.array, matrix after filtering.
    '''

    if   filter_type == "sobel":       data_filt = filters.sobel(data)
    elif filter_type == "roberts":     data_filt = filters.roberts(data)
    elif filter_type == "canny":       data_filt = feature.canny(data)

    elif filter_type == "lowpass_avg":
        p = int(filter_par)
        kernel = np.ones((p,p),np.float32)/(p*p)
        data_filt = ndimage.convolve(data, kernel)
    elif filter_type == "highpass_avg":
        p = int(filter_par)
        kernel = np.ones((p,p),np.float32)/(p*p)
        lp_data = ndimage.convolve(data, kernel)
        data_filt = data - lp_data

    elif filter_type == "lowpass_gaussian":
        data_filt = filters.gaussian(data, sigma=filter_par)
    elif filter_type == "highpass_gaussian":
        lp_data   = filters.gaussian(data, sigma=filter_par)
        data_filt = data - lp_data

    else:
        print('Un-recognized filter type: '+filter_type)
        sys.exit(1)

    return data_filt


############################################################
def filter_file(fname, filter_type, filter_par=None, fname_out=None):
    '''Filter 2D matrix with selected filter
    Inputs:
        fname       : string, name/path of file to be filtered
        filter_type : string, filter type
        filter_par  : string, optional, parameter for low/high pass filter
                      for low/highpass_avg, it's kernel size in int
                      for low/highpass_gaussain, it's sigma in float
    Output:
        fname_out   : string, optional, output file name/path
    '''

    # Basic info
    atr = readfile.read_attribute(fname)
    k = atr['FILE_TYPE']
    try:    ref_yx = [int(atr['REF_Y']), int(atr['REF_X'])]
    except: ref_yx = None

    filter_type = filter_type.lower()
    MSG = 'filtering '+k+' file: '+fname+' using '+filter_type+' filter'
    if filter_type.endswith('avg'):
        if not filter_par:
            filter_par = 5
        MSG += ' with kernel size of %d' % int(filter_par)
    elif filter_type.endswith('gaussian'):
        if not filter_par:
            filter_par = 3.0
        MSG += ' with sigma of %.1f' % filter_par
    print(MSG)

    if not fname_out:
        ext = os.path.splitext(fname)[1]
        fname_out = os.path.splitext(fname)[0]+'_'+filter_type+ext

    ##### Multiple Dataset File
    if k in multi_group_hdf5_file+multi_dataset_hdf5_file:
        h5 = h5py.File(fname,'r')
        epoch_list = sorted(h5[k].keys())
        epoch_num = len(epoch_list)
        prog_bar = ptime.progressBar(maxValue=epoch_num)

        h5out = h5py.File(fname_out,'w')
        group = h5out.create_group(k)
        print('writing >>> '+fname_out)

        if k == 'timeseries':
            print('number of acquisitions: '+str(epoch_num))
            for i in range(epoch_num):
                date = epoch_list[i]
                data = h5[k].get(date)[:]

                data_filt = filter_data(data, filter_type, filter_par)
                if ref_yx:
                    data_filt -= data_filt[ref_yx[0], ref_yx[1]]

                dset = group.create_dataset(date, data=data_filt)
                prog_bar.update(i+1, suffix=date)
            for key,value in iter(atr.items()):
                group.attrs[key] = value

        elif k in ['interferograms','wrapped','coherence']:
            print('number of interferograms: '+str(epoch_num))
            date12_list = ptime.list_ifgram2date12(epoch_list)
            for i in range(epoch_num):
                ifgram = epoch_list[i]
                data = h5[k][ifgram].get(ifgram)[:]

                data_filt = filter_data(data, filter_type, filter_par)
                if ref_yx and k in ['interferograms']:
                    data_filt -= data_filt[ref_yx[0], ref_yx[1]]

                gg = group.create_group(ifgram)
                dset = gg.create_dataset(ifgram, data=data_filt)
                for key, value in h5[k][ifgram].attrs.items():
                    gg.attrs[key] = value
                prog_bar.update(i+1, suffix=date12_list[i])

        h5.close()
        h5out.close()
        prog_bar.close()

    ##### Single Dataset File
    else:
        data, atr = readfile.read(fname)
        data_filt = filter_data(data, filter_type, filter_par)
        if ref_yx and k in ['.unw','velocity']:
            data_filt -= data_filt[ref_yx[0], ref_yx[1]]
        print('writing >>> '+fname_out)
        writefile.write(data_filt, out_file=fname_out, metadata=atr)

    return fname_out


################################################################################################
EXAMPLE='''example:
  spatial_filter.py  velocity.h5
  spatial_filter.py  timeseries.h5  lowpass_avg        5
  spatial_filter.py  velocity.h5    lowpass_avg        5
  spatial_filter.py  velocity.h5    highpass_gaussian  3
  spatial_filter.py  velocity.h5    sobel
'''

def create_parser():
    parser = argparse.ArgumentParser(description='Spatial filtering of 2D image.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='File to be filtered')
    parser.add_argument('filter_type', nargs='?', default='lowpass_gaussian',\
                        choices=['lowpass_gaussian','highpass_gaussian','lowpass_avg','highpass_avg',\
                                 'sobel','roberts','canny'],\
                        help='Type of filter. Default: lowpass_gaussian.\n'+\
                             'For more filters, check the link below:\n'+\
                             'http://scikit-image.org/docs/dev/api/skimage.filters.html')
    parser.add_argument('filter_par', nargs='?', type=float,\
                        help='Filter parameter for low/high pass filter. Default=\n'+\
                             'Sigma       for low/high pass gaussian filter, default: 3.0\n'+\
                             'Kernel Size for low/high pass average filter, default: 5')
    parser.add_argument('-o','--outfile', help='Output file name.')
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    inps.filter_type = inps.filter_type.lower()
    return inps


################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    print('Filter type: '+inps.filter_type)
    if inps.filter_type.startswith(('lowpass','highpass')):
        print('parameters: '+str(inps.filter_par))

    inps.outfile = filter_file(inps.file, inps.filter_type, inps.filter_par)
    print('Done.')
    return inps.outfile


################################################################################################
if __name__ == '__main__':
    main()



