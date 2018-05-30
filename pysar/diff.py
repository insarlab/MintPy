#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013-2018, Heresh Fattahi, Zhang Yunjun     #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################

import os
import sys
import argparse
import numpy as np
from pysar.objects import timeseries
from pysar.utils import readfile, writefile


#####################################################################################
def diff_data(data1, data2):
    """data1 - data2"""
    data = data1 - data2
    #data[np.isnan(data2)] = data1[np.isnan(data2)];
    return data


def check_reference(atr1, atr2):
    if atr1['REF_DATE'] == atr2['REF_DATE']:
        ref_date = None
    else:
        ref_date = atr1['REF_DATE']
        print('consider different reference date')

    ref_y = int(atr1['REF_Y'])
    ref_x = int(atr1['REF_X'])
    if ref_y == int(atr2['REF_Y']) and ref_x == int(atr2['REF_X']):
        ref_y = None
        ref_x = None
    else:
        print('consider different reference pixel')
    return ref_date, ref_y, ref_x


def diff_file(file1, file2, outFile=None, force=False):
    """Subtraction/difference of two input files"""
    if not outFile:
        fbase, fext = os.path.splitext(file1)
        outFile = '{}_diff_{}{}'.format(fbase, os.path.splitext(os.path.basename(file2))[0], fext)
    print('{} - {} --> {}'.format(os.path.basename(file1),
                                  os.path.basename(file2),
                                  os.path.basename(outFile)))

    # Read basic info
    atr1 = readfile.read_attribute(file1)
    k1 = atr1['FILE_TYPE']
    atr2 = readfile.read_attribute(file2)
    k2 = atr1['FILE_TYPE']
    print('input files are: {} and {}'.format(k1, k2))

    if k1 == 'timeseries':
        if k2 != 'timeseries':
            print('ERROR: input multiple dataset files are not the same file type!')
            sys.exit(1)

        obj1 = timeseries(file1)
        obj1.open()
        obj2 = timeseries(file2)
        obj2.open()
        ref_date, ref_y, ref_x = check_reference(obj1.metadata, obj2.metadata)

        # check dates shared by two timeseries files
        dateListShared = [i for i in obj1.dateList if i in obj2.dateList]
        dateShared = np.ones((obj1.numDate), dtype=np.bool_)
        if dateListShared != obj1.dateList:
            print('WARNING: {} does not contain all dates in {}'.format(file2, file1))
            if force:
                dateExcluded = list(set(obj1.dateList) - set(dateListShared))
                print('Continue and enforce the differencing for their shared dates only.')
                print('\twith following dates are ignored for differencing:\n{}'.format(dateExcluded))
                dateShared[np.array([obj1.dateList.index(i) for i in dateExcluded])] = 0
            else:
                raise Exception('To enforce the differencing anyway, use --force option.')

        # consider different reference_date/pixel
        data2 = obj2.read(dateListShared)
        if ref_date:
            data2 -= np.tile(data2[obj2.dateList.index(ref_date), :, :],
                             (data2.shape[0], 1, 1))
        if ref_y and ref_x:
            data2 -= np.tile(data2[:, ref_y, ref_x].reshape(-1, 1, 1),
                             (1, data2.shape[1], data2.shape[2]))

        data = obj1.read()
        mask = data == 0.
        data[dateShared] -= data2
        data[mask] = 0.               # Do not change zero phase value
        objOut = timeseries(outFile)
        objOut.write2hdf5(data=data, refFile=file1)

    # Sing dataset file
    else:
        data1, atr1 = readfile.read(file1)
        data2, atr2 = readfile.read(file2)
        data = data1 - data2
        print('writing >>> '+outFile)
        writefile.write(data, out_file=outFile, metadata=atr1)

    return outFile


#####################################################################################
EXAMPLE = """example:
  diff.py  velocity.h5      velocity_demCor.h5
  diff.py  timeseries.h5    ECMWF.h5  -o timeseries_ECMWF.h5
  diff.py  timeseries.h5    ECMWF.h5  -o timeseries_ECMWF.h5  --force
  diff.py  unwrapIfgram.h5  reconstruct_unwrapIfgram.h5
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Generates the difference of two input files.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file1', help='file to be substracted.')
    parser.add_argument('file2', help='file used to substract')
    parser.add_argument('-o', '--output', dest='outfile',
                        help='output file name, default is file1_diff_file2.h5')
    parser.add_argument('--force', action='store_true',
                        help='Enforce the differencing for the shared dates only for time-series files')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


def main(iargs=None):
    inps = cmd_line_parse(iargs)
    inps.outfile = diff_file(inps.file1, inps.file2, inps.outfile, force=inps.force)
    return inps.outfile


#####################################################################################
if __name__ == '__main__':
    main()
