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


############################################################
def print_attributes(atr, sorting=True):
    ## Print Dictionary of Attributes
    digits = digits = max([len(key) for key in atr.keys()]+[0])
    f = '{0:<%d}    {1}'%(digits)
    dictKey = atr.keys()
    if sorting:
        dictKey = sorted(dictKey)
    for key in dictKey:
        print((f.format(str(key),str(atr[key]))))
    return


## By andrewcollette at https://github.com/h5py/h5py/issues/406
def print_hdf5_structure(File):
    def print_hdf5_structure_obj(name, obj):
        print(name)
        print_attributes(obj.attrs)
    h5file=h5py.File(File,'r')
    h5file.visititems(print_hdf5_structure_obj)
    h5file.close()
    return


def print_timseries_date_info(dateList):
    datevector = ptime.date_list2vector(dateList)[1]
    print('*************** Date Info ***************')
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
def main(argv):
    inps = cmdLineParse()

    ext = os.path.splitext(inps.file)[1].lower()
    if inps.disp_date:
        atr = readfile.read_attribute(inps.file)
        k = atr['FILE_TYPE']
        if k not in ['timeseries','ifgramStack']:
            print('--date option can not be applied to {} file, ignore it.'.format(k))
            inps.disp_date = False
        else:
            f = h5py.File(inps.file,'r')
            dates = f[k].get('date')[:]
            print(dates)
        return

    try:
        atr = readfile.read_attribute(inps.file)
        k = atr['FILE_TYPE']
        print('file name: '+atr['FILE_PATH'])
        print('file type: '+atr['FILE_TYPE'])
        if 'Y_FIRST' in atr.keys():
            print('Coordinates : GEO')
        else:
            print('Coordinates : radar')
    except:
        pass

    if ext in ['.h5','.he5']:
        print('***** HDF5 File Structure *****')
        print_hdf5_structure(inps.file)
    else:
        atr = readfile.read_attribute(inps.file)
        print('***** File Attributes *********')
        print_attributes(atr)
    return


############################################################
if __name__ == '__main__':
    main(sys.argv[1:])

