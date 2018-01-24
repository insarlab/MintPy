#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2013, Heresh Fattahi, Zhang Yunjun          #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import os
import sys
import getopt
import time

import h5py
from numpy import std

import pysar._datetime as ptime
import pysar._readfile as readfile

output = ""
############################################################
def attributes_string(atr, string=str(), sorting=True):
    ## Print Dictionary of Attributes
    digits = digits = max([len(key) for key in atr.keys()] + [0])
    f = '{0:<%d}    {1}' % (digits)

    if sorting:
        keyList = sorted(atr.iterkeys())
    else:
        keyList = atr.iterkeys()
    for key in keyList:
        string += f.format(str(key), str(atr[key]))
        string += "\n"

    return string


def print_attributes(atr, string=str(), sorting=True):
    print(attributes_string(atr, string, sorting))


############################################################
def hdf5_structure_string(file):
    global output
    def print_hdf5_structure_obj(name, obj):
        global output
        output = attributes_string(obj.attrs, output)
        output += "\n"
        output += name

    h5file = h5py.File(file, 'r')
    h5file.visititems(print_hdf5_structure_obj)
    h5file.close()
    return output


## By andrewcollette at https://github.com/h5py/h5py/issues/406
def print_hdf5_structure(file):
    string = hdf5_structure_string(file)
    print(string)


############################################################
def print_timseries_date_info(dateList):
    datevector = ptime.date_list2vector(dateList)[1]
    print('*************** Date Info ***************')
    print('Start Date: ' + dateList[0])
    print('End   Date: ' + dateList[-1])
    print('Number of acquisitions      : %d' % len(dateList))
    print('Std.   of acquisition times : %.2f yeras' % std(datevector))
    print('----------------------')
    print('List of dates:')
    print(dateList)
    print('----------------------')
    print('List of dates in years')
    print(datevector)
    return


############################################################
def usage():
    print('usage: info.py file [eNum] [--tree/structure]\n'
          '\n'
          'Display the general information of File\n'
          '\n'
          'arguments:\n'
          '  file : HDF5 file, support all .h5 files\n'
          '  eNum : number of interferogram/coherence in the group\n'
          '         (1 as the first)\n'
          '  --struct/structure/tree : show the structure tree\n'
          '\n'
          'example:\n'
          '  info.py timeseries.h5\n'
          '  info.py velocity.h5\n'
          '  info.py unwrapIfgram.h5\n'
          '  info.py unwrapIfgram.h5    3\n'
          '\n'
          '  info.py timeseries.h5 --tree\n'
          '  info.py timeseries.h5 --date   # print out date list of timeseries HDF5 file\n'
          '    ')
    return


############################################################
def main(argv):
    ##### Check Inputs
    try:
        File = argv[0]
    except:
        usage();sys.exit(1)
    ext = os.path.splitext(File)[1].lower()

    #################### File Structure #####################
    try:
        print(argv[1])
        if argv[1] in ['--struct', '--structure', '--tree'] and ext in ['.h5', '.he5']:
            print('***** HDF5 File Structure *****')
            print_hdf5_structure(File)
            return
    except:
        pass

    #################### Basic Info #####################
    try:
        atr = readfile.read_attribute(File)
    except:
        print('Can not read file: ' + File)
        sys.exit(1)
    k = atr['FILE_TYPE']

    # Print out date list for timeseries HDF5 file
    try:
        if k in ['timeseries'] and argv[1] in ['--date']:
            h5 = h5py.File(File, 'r')
            dateList = h5[k].keys()
            for date in dateList:
                print(date)
            h5.close()
            return
    except: pass

    print('\n************************ File Info *****************************')
    print('File name   : '+os.path.basename(File))
    print('File type   : '+atr['PROCESSOR']+' '+atr['FILE_TYPE'])
    try:
        atr['X_FIRST']
        print('Coordinates : GEO')
    except:
        print('Coordinates : radar')

    print('\n************************ File Info *****************************')
    print('File name   : ' + os.path.basename(File))
    print('File type   : ' + atr['PROCESSOR'] + ' ' + atr['FILE_TYPE'])
    try:
        atr['X_FIRST']
        print('Coordinates : GEO')
    except:
        print('Coordinates : radar')

    #################### HDF5 File Info #####################
    if ext in ['.h5', '.he5']:
        h5file = h5py.File(File, 'r')
        ##### Group Info
        print('All groups in this file:')
        print(h5file.keys())

        ##### DateList / IgramList
        if k in ['interferograms', 'coherence', 'wrapped', 'timeseries']:
            epochList = sorted(h5file[k].keys())

    if k == 'timeseries':
        try:
            print_timseries_date_info(epochList)
        except:
            pass
        print('*************** Attributes **************')
        print_attributes(atr)

    elif k in ['interferograms', 'coherence', 'wrapped']:
        ##### Plot Attributes of One Epoch
        try:
            epochNum = int(argv[1])
            epochAtr = h5file[k][epochList[epochNum - 1]].attrs
            print('*****************************************')
            print(epochList[epochNum - 1])
            print('*************** Attributes **************')
            print_attributes(epochAtr)
            print('*****************************************')
            print(epochList[epochNum - 1])
        ##### Plot Epoch List Info
        except:

            print('*****************************************')
            print('Number of '+k+': '+str(len(epochList)) )
            print('*****************************************')
            print('List of the '+k+':             number')
            for i in range(len(epochList)):
                print(epochList[i] + '    ' + str(i + 1))
            print('*****************************************')
            print('Number of ' + k + ': ' + str(len(epochList)))

    ##### All other file types, except for timeseries/interferograms/coherence/wrapped
    else:
        print('*************** Attributes **************')
        print_attributes(atr)

    try:
        h5file.close()
    except:
        pass
    print('****************************************************************')
    return


############################################################
if __name__ == '__main__':
    main(sys.argv[1:])
