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
from mintpy.objects import timeseries
from mintpy.utils import readfile, writefile
from mintpy.diff import check_reference


################################################################################
EXAMPLE = """example:
  add.py mask_1.h5 mask_2.h5 mask_3.h5        -o mask_all.h5
  add.py 081008_100220.unw  100220_110417.unw -o 081008_110417.unw
  add.py timeseries_ERA5.h5 inputs/ERA5.h5    -o timeseries.h5
  add.py timeseriesRg.h5    inputs/TECsub.h5  -o timeseriesRg_TECsub.h5 --force
"""


def create_parser():
    """ Command line parser """
    parser = argparse.ArgumentParser(description='Generate sum of multiple input files.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs='+', help='files (2 or more) to be added')
    parser.add_argument('-o', '--output', dest='outfile', help='output file name')
    parser.add_argument('--force', action='store_true',
                        help='Enforce the adding for the shared dates only for time-series files')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    if len(inps.file) < 2:
        parser.print_usage()
        sys.exit('ERROR: At least 2 input files needed!')

    # only two input files are supported for time-series type
    atr = readfile.read_attribute(inps.file[0])
    if atr['FILE_TYPE'] == 'timeseries' and len(inps.file) != 2:
        raise ValueError('Only TWO files are supported for time-series, input has {}'.format(len(inps.file)))

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


def add_file(fnames, out_file=None, force=False):
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

    atr1 = readfile.read_attribute(fnames[0])
    atr2 = readfile.read_attribute(fnames[1])

    if atr1['FILE_TYPE'] == 'timeseries':
        file1, file2 = fnames[0], fnames[1]

        # check dates shared by two timeseries files
        dateList1 = timeseries(file1).get_date_list()
        dateList2 = timeseries(file2).get_date_list()
        dateListShared = [i for i in dateList1 if i in dateList2]
        dateShared = np.ones((len(dateList1)), dtype=np.bool_)
        if dateListShared != dateList1:
            print('WARNING: {} does not contain all dates in {}'.format(file2, file1))
            if force:
                dateListEx = list(set(dateList1) - set(dateListShared))
                print('Continue and enforce the differencing for their shared dates only.')
                print('\twith following dates are ignored for differencing:\n{}'.format(dateListEx))
                dateShared[np.array([dateList1.index(i) for i in dateListEx])] = 0
            else:
                raise Exception('To enforce the differencing anyway, use --force option.')

        # check reference point
        ref_date, ref_y, ref_x = check_reference(atr1, atr2)

        # read data2 (consider different reference_date/pixel)
        print('read from file: {}'.format(file2))
        data2 = readfile.read(file2, datasetName=dateListShared)[0]

        if ref_y and ref_x:
            print('* referencing data from {} to y/x: {}/{}'.format(os.path.basename(file2), ref_y, ref_x))
            ref_box = (ref_x, ref_y, ref_x + 1, ref_y + 1)
            ref_val = readfile.read(file2, datasetName=dateListShared, box=ref_box)[0]
            data2 -= np.tile(ref_val.reshape(-1, 1, 1), (1, data2.shape[1], data2.shape[2]))

        if ref_date:
            print('* referencing data from {} to date: {}'.format(os.path.basename(file2), ref_date))
            ref_ind = dateListShared.index(ref_date)
            data2 -= np.tile(data2[ref_ind, :, :], (data2.shape[0], 1, 1))

        # read data1
        print('read from file: {}'.format(file1))
        data = readfile.read(file1)[0]

        # apply adding
        mask = data == 0.
        data[dateShared] += data2
        data[mask] = 0.               # Do not change zero phase value
        del data2

        # write file
        writefile.write(data, out_file=out_file, metadata=atr1, ref_file=file1)

    else:
        dsDict = {}
        dsNames = readfile.get_dataset_list(fnames[0])
        for dsName in dsNames:
            # ignore dsName if input file has single dataset
            dsName2read = None if len(dsNames) == 1 else dsName

            print('adding {} ...'.format(dsName))
            data = readfile.read(fnames[0], datasetName=dsName2read)[0]
            for i in range(1, len(fnames)):
                data2 = readfile.read(fnames[i], datasetName=dsName2read)[0]
                data = add_matrix(data, data2)
            dsDict[dsName] = data

        # output
        print('use metadata from the 1st file: {}'.format(fnames[0]))
        writefile.write(dsDict, out_file=out_file, metadata=atr1, ref_file=fnames[0])

    return out_file


################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    print('input files to be added: ({})\n{}'.format(len(inps.file), inps.file))

    inps.outfile = add_file(inps.file, inps.outfile, force=inps.force)

    print('Done.')
    return inps.outfile


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
