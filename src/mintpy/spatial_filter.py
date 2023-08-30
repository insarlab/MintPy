############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, Zhang Yunjun, 2013               #
############################################################
# Add "double_difference" filter, Forrest Williams, May 2021


import os
import sys

import numpy as np
from scipy import ndimage
from skimage import feature, filters, morphology

from mintpy.utils import readfile, writefile


################################################################################################
def filter_data(data, filter_type, filter_par=None):
    """Filter 2D matrix with selected filter.

    Parameters: data        - 2D np.array, matrix to be filtered
                filter_type - string, filter type
                filter_par  - (list of) int/float, optional, parameter for low/high pass filter
                              for low/highpass_avg, it's kernel size in int
                              for low/highpass_gaussain, it's sigma in float
                              for double_difference, it's local and regional kernel sizes in int
    Returns:    data_filt   - 2D np.array, matrix after filtering.
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
        # ORIGINAL: data_filt = filters.gaussian(data, sigma=filter_par)
        #   nan pixels can enlarge to big holes depending on the size of your gaussian kernel
        #   we can do normalized convolution (https://stackoverflow.com/a/36307291/7128154) as below:
        V=np.array(data)
        V[np.isnan(data)]=0
        VV=filters.gaussian(V, sigma=filter_par)

        W=np.ones_like(data)
        W[np.isnan(data)]=0
        WW=filters.gaussian(W, sigma=filter_par)
        WW[WW==0]=np.nan

        data_filt = VV/WW
        data_filt[np.isnan(data)] = np.nan

    elif filter_type == "highpass_gaussian":
        # Use the existing logic for lowpass_gaussian, then subtract from original data
        return data - filter_data(data, "lowpass_gaussian", filter_par=filter_par)

    elif filter_type == "median":
        data_filt = filters.median(data, morphology.disk(filter_par))

    elif filter_type == "double_difference":
        """Amplifies the local deformation signal by reducing the influence
        of regional deformation trends from atmospheric artifacts, tectonic
        deformation, and other sources. Intended use is to identify landslide-related
        deformation. Filter has the form:

            result =  regional mean of data - local mean of data

        where both the regional and local kernel size can be set using the
        filter_par argument. The regional kernel has a doughnut shape
        because the pixels used to calculate the local mean are not
        included in the regional mean. This ensures that the local mean
        and regional mean results are separate.
        """

        local_kernel = morphology.disk(filter_par[0], dtype=np.float32)
        local_kernel = np.pad(local_kernel, pad_width=int(filter_par[1]-filter_par[0]), mode='constant')

        regional_kernel = morphology.disk(filter_par[1], dtype=np.float32)
        regional_kernel[local_kernel == 1] = 0

        # normalize
        local_kernel /= local_kernel.sum(axis=(0,1))
        regional_kernel /= regional_kernel.sum(axis=(0,1))

        # double-difference kernel
        combined_kernel = regional_kernel - local_kernel

        data_filt = ndimage.convolve(data, combined_kernel)

    else:
        raise Exception('Un-recognized filter type: '+filter_type)

    return data_filt


################################################################################################
def filter_file(fname, ds_names=None, filter_type='lowpass_gaussian', filter_par=None, fname_out=None):
    """Filter 2D matrix with selected filter.

    Parameters: fname       - string, name/path of file to be filtered
                ds_names    - list of string, datasets of interest
                filter_type - string, filter type
                filter_par  - (list of) int/float, optional, parameter for low/high pass filter
                              for low/highpass_avg, it's kernel size in int
                              for low/highpass_gaussain, it's sigma in float
                              for double_difference, it's local and regional kernel sizes in int
    Returns:    fname_out   - string, optional, output file name/path
    """

    # Info
    filter_type = filter_type.lower()
    atr = readfile.read_attribute(fname)
    print(f'input file: {fname}, file type: {atr["FILE_TYPE"]}')
    print(f'filter type: {filter_type}')

    if filter_type.endswith('avg'):
        if not filter_par:
            filter_par = 5
        elif isinstance(filter_par, list):
            filter_par = filter_par[0]
        filter_par = int(filter_par)
        print(f'filter parameter: kernel size = {filter_par}')

    elif filter_type.endswith('gaussian'):
        if not filter_par:
            filter_par = 3.0
        elif isinstance(filter_par, list):
            filter_par = filter_par[0]
        filter_par = float(filter_par)
        print(f'filter parameter: sigma = {filter_par:.1f}')

    elif filter_type == 'double_difference':
        if not filter_par:
            filter_par = [1, 10]
        local, regional = int(filter_par[0]), int(filter_par[1])
        print(f'filter parameter: local / regional kernel sizes = {local} / {regional}')

    elif filter_type == 'median':
        if not filter_par:
            filter_par = 5
        elif isinstance(filter_par, list):
            filter_par = filter_par[0]
        print(f'filter parameter:  median radius of {filter_par} pixels')


    # output filename
    if not fname_out:
        fbase, fext = os.path.splitext(fname)
        fname_out = f'{fbase}_{filter_type}{fext}'

    # filtering file
    ds_all = readfile.get_dataset_list(fname)
    ds_names = ds_names if ds_names else ds_all
    ds_skips = list(set(ds_all) - set(ds_names))

    maxDigit = max(len(i) for i in ds_names)
    dsDict = dict()

    for ds_name in ds_skips:
        dsDict[ds_name] = readfile.read(fname, datasetName=ds_name, print_msg=False)[0]

    # loop over each dataset
    for ds_name in ds_names:
        msg = 'filtering {d:<{w}} from {f} '.format(d=ds_name, w=maxDigit, f=os.path.basename(fname))
        # read
        data = readfile.read(fname, datasetName=ds_name, print_msg=False)[0]

        # filter
        if len(data.shape) == 3:
            # 3D matrix
            num_loop = data.shape[0]
            for i in range(num_loop):
                data[i, :, :] = filter_data(data[i, :, :], filter_type, filter_par)
                sys.stdout.write(f'\r{msg} {i+1}/{num_loop} ...')
                sys.stdout.flush()
            print('')

        else:
            # 2D matrix
            data = filter_data(data, filter_type, filter_par)

        # save
        dsDict[ds_name] = data

    # write to file
    writefile.write(
        dsDict,
        out_file=fname_out,
        metadata=atr,
        ref_file=fname,
    )

    print('Done.')
    return fname_out
