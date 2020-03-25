#!/usr/bin/env python
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################

import os
import time
import argparse
import numpy as np
from mintpy.objects import timeseries, giantTimeseries, ifgramStack, ifgramDatasetNames
from mintpy.utils import readfile, writefile


#####################################################################################
EXAMPLE = """example:
  diff.py  velocity.h5    velocity_demErr.h5
  diff.py  timeseries.h5  inputs/ERA5.h5  -o timeseries_ERA5.h5
  diff.py  timeseries.h5  inputs/ERA5.h5  -o timeseries_ERA5.h5  --force
  diff.py  timeseries_ERA5_ramp_demErr.h5  ../GIANT/Stack/LS-PARAMS.h5 -o mintpy_giant.h5
  diff.py  reconUnwrapIfgram.h5  ./inputs/ifgramStack.h5  -o diffUnwrapIfgram.h5

  # multiple files
  diff.py  waterMask.h5  maskSantiago.h5  maskFernandina.h5  -o maskIsabela.h5
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Generates the difference of two input files.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file1', help='file to be substracted.')
    parser.add_argument('file2', nargs='+', help='file used to substract')
    parser.add_argument('-o', '--output', dest='outfile',
                        help='output file name, default is file1_diff_file2.h5')
    parser.add_argument('--force', action='store_true',
                        help='Enforce the differencing for the shared dates only for time-series files')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # for timeseries and ifgramStack, only two files differencing is supported
    atr = readfile.read_attribute(inps.file1)
    if atr['FILE_TYPE'] in ['timeseries', 'ifgramStack']:
        if len(inps.file2) > 1:
            raise SystemExit('ERROR: only one file2 is inputed for {} type'.format(atr['FILE_TYPE']))
    return inps


#####################################################################################
def _check_reference(atr1, atr2):
    """Check reference date and point
    Parameters: atr1/2   - dict, metadata of file1/2
    Returns:    ref_date - str, None for re-referencing in time  is NOT needed
                ref_y/x  - int, None for re-referencing in space is NOT needed
    """
    # reference date
    if atr1['REF_DATE'] == atr2.get('REF_DATE', None):
        ref_date = None
    else:
        ref_date = atr1['REF_DATE']
        #print('consider different reference date')

    # reference point
    ref_y = atr1.get('REF_Y', None)
    ref_x = atr1.get('REF_X', None)
    if ref_x == atr2.get('REF_X', None) or ref_y == atr2.get('REF_Y', None):
        ref_y = None
        ref_x = None
    else:
        ref_y = ref_y
        ref_x = ref_x
        #print('consider different reference pixel')

    if ref_y is not None:
        ref_y = int(ref_y)
    if ref_x is not None:
        ref_x = int(ref_x)
    return ref_date, ref_y, ref_x


def diff_file(file1, file2, outFile=None, force=False):
    """Subtraction/difference of two input files"""
    if not outFile:
        fbase, fext = os.path.splitext(file1)
        if len(file2) > 1:
            raise ValueError('Output file name is needed for more than 2 files input.')
        outFile = '{}_diff_{}{}'.format(fbase, os.path.splitext(os.path.basename(file2[0]))[0], fext)
    print('{} - {} --> {}'.format(file1, file2, outFile))

    # Read basic info
    atr1 = readfile.read_attribute(file1)
    k1 = atr1['FILE_TYPE']
    atr2 = readfile.read_attribute(file2[0])
    k2 = atr2['FILE_TYPE']
    print('input files are: {} and {}'.format(k1, k2))

    if k1 == 'timeseries':
        if k2 not in ['timeseries', 'giantTimeseries']:
            raise Exception('Input multiple dataset files are not the same file type!')
        if len(file2) > 1:
            raise Exception(('Only 2 files substraction is supported for time-series file,'
                             ' {} input.'.format(len(file2)+1)))

        atr1 = readfile.read_attribute(file1)
        atr2 = readfile.read_attribute(file2[0])
        dateList1 = timeseries(file1).get_date_list()
        if k2 == 'timeseries':
            dateList2 = timeseries(file2[0]).get_date_list()
            unit_fac = 1.
        elif k2 == 'giantTimeseries':
            dateList2 = giantTimeseries(file2[0]).get_date_list()
            unit_fac = 0.001
        ref_date, ref_y, ref_x = _check_reference(atr1, atr2)

        # check dates shared by two timeseries files
        dateListShared = [i for i in dateList1 if i in dateList2]
        dateShared = np.ones((len(dateList1)), dtype=np.bool_)
        if dateListShared != dateList1:
            print('WARNING: {} does not contain all dates in {}'.format(file2, file1))
            if force:
                dateExcluded = list(set(dateList1) - set(dateListShared))
                print('Continue and enforce the differencing for their shared dates only.')
                print('\twith following dates are ignored for differencing:\n{}'.format(dateExcluded))
                dateShared[np.array([dateList1.index(i) for i in dateExcluded])] = 0
            else:
                raise Exception('To enforce the differencing anyway, use --force option.')

        # consider different reference_date/pixel
        data2 = readfile.read(file2[0], datasetName=dateListShared)[0] * unit_fac
        if ref_date:
            print('* referencing data from {} to date: {}'.format(os.path.basename(file2[0]), ref_date))
            data2 -= np.tile(data2[dateListShared.index(ref_date), :, :],
                             (data2.shape[0], 1, 1))
        if ref_y and ref_x:
            print('* referencing data from {} to y/x: {}/{}'.format(os.path.basename(file2[0]), ref_y, ref_x))
            data2 -= np.tile(data2[:, ref_y, ref_x].reshape(-1, 1, 1),
                             (1, data2.shape[1], data2.shape[2]))

        data = readfile.read(file1)[0]
        mask = data == 0.
        data[dateShared] -= data2
        data[mask] = 0.               # Do not change zero phase value
        del data2
        writefile.write(data, out_file=outFile, ref_file=file1)

    elif all(i == 'ifgramStack' for i in [k1, k2]):
        obj1 = ifgramStack(file1)
        obj1.open()
        obj2 = ifgramStack(file2[0])
        obj2.open()
        dsNames = list(set(obj1.datasetNames) & set(obj2.datasetNames))
        if len(dsNames) == 0:
            raise ValueError('no common dataset between two files!')
        dsName = [i for i in ifgramDatasetNames if i in dsNames][0]

        # read data
        print('reading {} from file {} ...'.format(dsName, file1))
        data1 = readfile.read(file1, datasetName=dsName)[0]
        print('reading {} from file {} ...'.format(dsName, file2[0]))
        data2 = readfile.read(file2[0], datasetName=dsName)[0]

        # consider reference pixel
        if 'unwrapphase' in dsName.lower():
            print('referencing to pixel ({},{}) ...'.format(obj1.refY, obj1.refX))
            ref1 = data1[:, obj1.refY, obj1.refX]
            ref2 = data2[:, obj2.refY, obj2.refX]
            for i in range(data1.shape[0]):
                data1[i,:][data1[i, :] != 0.] -= ref1[i]
                data2[i,:][data2[i, :] != 0.] -= ref2[i]

        # operation and ignore zero values
        data1[data1 == 0] = np.nan
        data2[data2 == 0] = np.nan
        data = data1 - data2
        del data1, data2
        data[np.isnan(data)] = 0.

        # write to file
        dsDict = {}
        dsDict[dsName] = data
        writefile.write(dsDict, out_file=outFile, ref_file=file1)

    # Sing dataset file
    else:
        data1 = readfile.read(file1)[0]
        data = np.array(data1, data1.dtype)
        for fname in file2:
            data2 = readfile.read(fname)[0]
            data = np.array(data, dtype=np.float32) - np.array(data2, dtype=np.float32)
            data = np.array(data, data1.dtype)
        print('writing >>> '+outFile)
        writefile.write(data, out_file=outFile, metadata=atr1)

    return outFile


def main(iargs=None):
    inps = cmd_line_parse(iargs)
    #start_time = time.time()

    inps.outfile = diff_file(inps.file1, inps.file2, inps.outfile, force=inps.force)

    #m, s = divmod(time.time()-start_time, 60)
    #print('time used: {:02.0f} mins {:02.1f} secs'.format(m, s))
    return inps.outfile


#####################################################################################
if __name__ == '__main__':
    main()
