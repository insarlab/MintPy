#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013-2018, Zhang Yunjun, Heresh Fattahi     #
# Author:  Zhang Yunjun, Heresh Fattahi                    #
############################################################


import os
import sys
import argparse
import time
import h5py
from numpy import std
from pysar.utils import readfile, ptime
from pysar.objects import timeseries, ifgramStack, geometry, HDFEOS


output = ""

############################################################
EXAMPLE = """example:
  info.py timeseries.h5
  info.py velocity.h5
  info.py ifgramStack.h5

  info.py ifgramStack.h5 --date                  # print master/slave date pairs info of interferograms.
  info.py timeseries.h5 --date                   # print date list of timeseries.
  info.py timeseries.h5 --date > date_list.txt   # print date list of timeseries and save it to txt file.
"""


def create_parser():
    """Create command line parser."""
    parser = argparse.ArgumentParser(description='Display Metadata / Structure information of File',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)
    parser.add_argument('file', type=str, help='File to check')
    parser.add_argument('--date', dest='disp_date', action='store_true',
                        help='Show date/date12 info of input file')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


############################################################
def attributes_string(atr, string=str(), sorting=True):
    ## Print Dictionary of Attributes
    digits = max([len(key) for key in list(atr.keys())] + [0])
    f = '{0:<%d}    {1}' % (digits)
    dictKey = atr.keys()
    if sorting:
        dictKey = sorted(dictKey)
    for key in dictKey:
        string += '  {k:<{d}}    {v}'.format(d=digits, k=key, v=atr[key])
        string += "\n"
    return string


def print_attributes(atr, string=str(), sorting=True):
    print((attributes_string(atr, string, sorting)))


############################################################
def hdf5_structure_string(file):
    global output

    def print_hdf5_structure_obj(name, obj):
        global output
        maxDigit = 25
        if name.startswith('HDFEOS'):
            maxDigit += 30

        if isinstance(obj, h5py.Group):
            output += 'HDF5 group   "/{n}"\n'.format(n=name)
        elif isinstance(obj, h5py.Dataset):
            output += ('HDF5 dataset "/{n:<{w}}": shape {s:<20}, '
                       'dtype <{t}>\n').format(n=name,
                                               w=maxDigit,
                                               s=str(obj.shape),
                                               t=obj.dtype)

    f = h5py.File(file, 'r')
    if len(f.attrs) > 0:
        output += 'Attributes in / level:\n'
        output = attributes_string(f.attrs, output, sorting=True)+"\n"

    f.visititems(print_hdf5_structure_obj)
    f.close()

    local_output = output
    output = ""
    return local_output


## By andrewcollette at https://github.com/h5py/h5py/issues/406
def print_hdf5_structure(file):
    string = hdf5_structure_string(file)
    print(string)


def print_timseries_date_stat(dateList):
    datevector = ptime.date_list2vector(dateList)[1]
    print('Start Date: '+dateList[0])
    print('End   Date: '+dateList[-1])
    print('Number of acquisitions    : %d' % len(dateList))
    print('Std. of acquisition times : %.2f yeras' % std(datevector))
    print('----------------------')
    print('List of dates:')
    print(dateList)
    print('----------------------')
    print('List of dates in years')
    print(datevector)
    return


def get_date_list(fname, print_msg=False):
    atr = readfile.read_attribute(fname)
    k = atr['FILE_TYPE']
    dateList = None
    if k in ['timeseries']:
        obj = timeseries(fname)
        obj.open(print_msg=False)
        dateList = obj.dateList
    elif k in ['ifgramStack']:
        obj = ifgramStack(fname)
        obj.open(print_msg=False)
        dateList = obj.date12List
    else:
        print('--date option can not be applied to {} file, ignore it.'.format(k))
    try:
        obj.close(print_msg=False)
    except:
        pass

    if print_msg and dateList is not None:
        for i in dateList:
            print(i)
    return dateList


def print_pysar_info(fname):
    try:
        atr = readfile.read_attribute(fname)
        k = atr['FILE_TYPE']
        print('{} {:*<40}'.format('*'*20, 'Basic File Info '))
        print('file name: '+atr['FILE_PATH'])
        print('file type: '+atr['FILE_TYPE'])
        if 'Y_FIRST' in atr.keys():
            print('coordinates : GEO')
        else:
            print('coordinates : RADAR')
        if k in ['timeseries']:
            dateList = get_date_list(fname)
            print('\n{} {:*<40}'.format('*'*20, 'Date Stat Info '))
            print_timseries_date_stat(dateList)
    except:
        pass
    return


############################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    if not os.path.isfile(inps.file):
        print('ERROR: input file does not exists: {}'.format(inps.file))
        return
    ext = os.path.splitext(inps.file)[1].lower()

    # --date option
    if inps.disp_date:
        get_date_list(inps.file, print_msg=True)
        return

    # Basic info from PySAR reader
    print_pysar_info(inps.file)

    # Generic Attribute/Structure of all files
    if ext in ['.h5', '.he5']:
        print('\n{} {:*<40}'.format('*'*20, 'HDF5 File Structure '))
        print_hdf5_structure(inps.file)
    else:
        print('\n{} {:*<40}'.format('*'*20, 'Binary File Attributes '))
        atr = readfile.read_attribute(inps.file)
        print_attributes(atr)

    return


############################################################
if __name__ == '__main__':
    main()
