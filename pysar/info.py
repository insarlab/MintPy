#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2013, Heresh Fattahi, Zhang Yunjun          #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import os, sys
import argparse
import time
import h5py
from numpy import std
from pysar.utils import readfile, datetime as ptime
from pysar.objects import timeseries, ifgramStack, geometry, hdfEos5


############################################################
EXAMPLE='''example:
  info.py timeseries.h5
  info.py velocity.h5
  info.py ifgramStack.h5

  info.py ifgramStack.h5 --date                  # print master/slave date pairs info of interferograms.
  info.py timeseries.h5 --date                   # print date list of timeseries.
  info.py timeseries.h5 --date > date_list.txt   # print date list of timeseries and save it to txt file.
'''

def createParser():
    '''Create command line parser.'''
    parser = argparse.ArgumentParser(description='Display Metadata / Structure information of File',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)
    parser.add_argument('file', type=str, help='File to check')
    parser.add_argument('--date', dest='disp_date', action='store_true', help='Show date/date12 info of input file')
    return parser


def cmdLineParse(iargs = None):
    '''Command line parser.'''
    parser = createParser()
    inps = parser.parse_args(args=iargs)
    return inps


############################################################
def print_attributes(atr, sorting=True):
    ## Print Dictionary of Attributes
    digits = digits = max([len(key) for key in atr.keys()]+[0])
    f = '{0:<%d}    {1}'%(digits)
    dictKey = atr.keys()
    if sorting:
        dictKey = sorted(dictKey)
    for key in dictKey:
        print('  {k:<{d}}    {v}'.format(d=digits, k=key, v=atr[key]))
    return


def print_hdf5_structure(File):
    '''Modified from andrewcollette at https://github.com/h5py/h5py/issues/406'''
    def print_hdf5_structure_obj(name, obj):
        if isinstance(obj, h5py.Group):
            print('HDF5 group "/{}"'.format(name))
        elif isinstance(obj, h5py.Dataset):
            print('HDF5 dataset "/{:<30}": shape {:<20}, dtype <{}>'.format(name, str(obj.shape), obj.dtype))
        print_attributes(obj.attrs)
    f=h5py.File(File,'r')
    f.visititems(print_hdf5_structure_obj)
    f.close()
    return


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


def get_date_list(fname, printMsg=False):
    atr = readfile.read_attribute(fname)
    k = atr['FILE_TYPE']
    dateList = None
    if k in ['timeseries']:
        obj = timeseries(fname)
        obj.open(printMsg=False)
        dateList = obj.dateList
    elif k in ['ifgramStack']:
        obj = ifgramStack(fname)
        obj.open(printMsg=False)
        dateList = obj.date12List
    else:
        print('--date option can not be applied to {} file, ignore it.'.format(k))
    try: obj.close(printMsg=False)
    except: pass

    if printMsg and dateList is not None:
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
def main(argv):
    inps = cmdLineParse()
    if not os.path.isfile(inps.file):
        print('ERROR: input file does not exists: {}'.format(inps.file))
        return
    ext = os.path.splitext(inps.file)[1].lower()

    ## --date option
    if inps.disp_date:
        get_date_list(inps.file, printMsg=True)
        return

    ## Basic info from PySAR reader
    print_pysar_info(inps.file)

    ## Generic Attribute/Structure of all files
    if ext in ['.h5','.he5']:
        print('\n{} {:*<40}'.format('*'*20, 'HDF5 File Structure '))
        print_hdf5_structure(inps.file)
    else:
        print('\n{} {:*<40}'.format('*'*20, 'Binary File Attributes '))
        atr = readfile.read_attribute(inps.file)
        print_attributes(atr)

    return


############################################################
if __name__ == '__main__':
    main(sys.argv[1:])

