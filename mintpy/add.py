############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os

import numpy as np

from mintpy.diff import check_reference
from mintpy.objects import timeseries
from mintpy.utils import readfile, writefile


################################################################################
def add_matrix(data1, data2):
    """Sum of 2 input matrix"""
    data1 = np.array(data1, dtype=np.float32)
    data2 = np.array(data2, dtype=np.float32)
    data = data1 + data2
    data[np.isnan(data1)] = data2[np.isnan(data1)]
    data[np.isnan(data2)] = data1[np.isnan(data2)]
    return data


def add_file(fnames, out_file=None, force=False):
    """Generate sum of all input files
    Parameters: fnames   - list of str, path/name of input files to be added
                out_file - str, optional, path/name of output file
    Returns:    out_file - str, path/name of output file
    Example:    'mask_all.h5' = add_file(['mask_1.h5','mask_2.h5','mask_3.h5'], 'mask_all.h5')
    """
    num_file = len(fnames)
    print(f'input files to be added: ({num_file})\n{fnames}')

    # Default output file name
    ext = os.path.splitext(fnames[0])[1]
    if not out_file:
        out_file = os.path.splitext(fnames[0])[0]
        for i in range(1, num_file):
            out_file += '_plus_' + os.path.splitext(os.path.basename(fnames[i]))[0]
        out_file += ext

    # read FILE_TYPE
    ftypes = [readfile.read_attribute(x)['FILE_TYPE'] for x in fnames]
    print(f'input file types: {ftypes}')

    if ftypes[0] == 'timeseries':
        # check dates shared by two timeseries files
        file1, file2 = fnames[0], fnames[1]
        atr1 = readfile.read_attribute(file1)
        atr2 = readfile.read_attribute(file2)
        dateList1 = timeseries(file1).get_date_list()
        dateList2 = timeseries(file2).get_date_list()
        dateListShared = [i for i in dateList1 if i in dateList2]
        dateShared = np.ones((len(dateList1)), dtype=np.bool_)
        if dateListShared != dateList1:
            print(f'WARNING: {file2} does not contain all dates in {file1}')
            if force:
                dateListEx = list(set(dateList1) - set(dateListShared))
                print('Continue and enforce the differencing for their shared dates only.')
                print(f'\twith following dates are ignored for differencing:\n{dateListEx}')
                dateShared[np.array([dateList1.index(i) for i in dateListEx])] = 0
            else:
                raise Exception('To enforce the differencing anyway, use --force option.')

        # check reference point
        ref_date, ref_y, ref_x = check_reference(atr1, atr2)

        # read data2 (consider different reference_date/pixel)
        print(f'read from file: {file2}')
        data2 = readfile.read(file2, datasetName=dateListShared)[0]

        if ref_y and ref_x:
            print(f'* referencing data from {os.path.basename(file2)} to y/x: {ref_y}/{ref_x}')
            ref_box = (ref_x, ref_y, ref_x + 1, ref_y + 1)
            ref_val = readfile.read(file2, datasetName=dateListShared, box=ref_box)[0]
            data2 -= np.tile(ref_val.reshape(-1, 1, 1), (1, data2.shape[1], data2.shape[2]))

        if ref_date:
            print(f'* referencing data from {os.path.basename(file2)} to date: {ref_date}')
            ref_ind = dateListShared.index(ref_date)
            data2 -= np.tile(data2[ref_ind, :, :], (data2.shape[0], 1, 1))

        # read data1
        print(f'read from file: {file1}')
        data = readfile.read(file1)[0]

        # apply adding
        mask = data == 0.
        data[dateShared] += data2
        data[mask] = 0.               # Do not change zero phase value
        del data2

        # write file
        writefile.write(data, out_file=out_file, metadata=atr1, ref_file=file1)

    else:
        # get common dataset list
        ds_names_list = [readfile.get_dataset_list(x) for x in fnames]
        ds_names = list(set.intersection(*map(set, ds_names_list)))
        # if all files have one dataset, ignore dataset name variation and take the 1st one as reference
        if all(len(x) == 1 for x in ds_names_list):
            ds_names = ds_names_list[0]
        print('List of common datasets across files: ', ds_names)
        if len(ds_names) < 1:
            raise ValueError(f'No common datasets found among files:\n{fnames}')

        # loop over each file
        dsDict = {}
        for ds_name in ds_names:
            print(f'adding {ds_name} ...')
            data, atr = readfile.read(fnames[0], datasetName=ds_name)

            for i, fname in enumerate(fnames[1:]):
                # ignore ds_name if input file has single dataset
                ds_name2read = None if len(ds_names_list[i+1]) == 1 else ds_name
                # read
                data2 = readfile.read(fname, datasetName=ds_name2read)[0]
                # apply operation
                data = add_matrix(data, data2)
            dsDict[ds_name] = data

        # output
        print(f'use metadata from the 1st file: {fnames[0]}')
        writefile.write(dsDict, out_file=out_file, metadata=atr, ref_file=fnames[0])

    return out_file
