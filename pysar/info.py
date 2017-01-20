#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#
# Yunjun, Jul 2015: add 'coherence','wrapped' option, eNum
# Yunjun, Oct 2015: merge 'interferograms','coherence','wrapped' into one
# Yunjun, Oct 2016: add --tree option to check hdf5 tree/structure


import os
import sys
import getopt
import time

import h5py
from numpy import std

import pysar._readfile as readfile
import pysar._datetime as ptime


############################################################
def print_attributes(atr, sort=1):
    ## Print Dictionary of Attributes
    keyDigitList   = []
    for key in atr.iterkeys():
        keyDigitList.append(len(key))
    digits = max(keyDigitList+[0])
    f = '{0:<%d}    {1}'%(digits)

    if sort == 1: keyList = sorted(atr.iterkeys())
    else:         keyList = atr.iterkeys()
    for key in keyList:
        print(f.format(str(key),str(atr[key])))

    return


############################################################
## By andrewcollette at https://github.com/h5py/h5py/issues/406
def print_hdf5_structure(File):
    def print_hdf5_structure_obj(name, obj):
        print name
        print_attributes(obj.attrs)
    h5file=h5py.File(File,'r')
    h5file.visititems(print_hdf5_structure_obj)
    h5file.close()
    return


############################################################
def print_timseries_date_info(dateList):
    datevector = ptime.date_list2vector(dateList)[1]
    print '*************** Date Info ***************'
    print 'Start       Date: '+dateList[0]
    print 'End         Date: '+dateList[-1]
    print 'Number             of acquisitions      : '+str(len(dateList))
    print 'Standard deviation of acquisition times : '+str(std(datevector))+' years'
    print '----------------------'
    print 'List of dates:'
    print dateList
    print '----------------------'
    print 'List of dates in years'
    print datevector
    return

############################################################
def usage():
    print '''
***************************************************************
  Displayes the general information of the PySAR product h5 file.

  Usage:
      info.py hdf5File  [eNum]

      file : HDF5 file, support all .h5 files
      eNum : number of interferogram/coherence in the group
             (1 as the first)
      --struct/structure/tree : show the structure tree

  Example:

      info.py timeseries.h5
      info.py velocity_demCor_masked.h5
      info.py LoadedData.h5
      info.py LoadedData.h5    3

      info.py timeseries.h5 --tree

***************************************************************
    '''

############################################################
def main(argv):

    ##### Check Inputs
    try:    File = argv[0]
    except: usage();sys.exit(1)
    print '\n************************ File Info *****************************'

    #################### Basic Info #####################
    try: atr = readfile.read_attribute(File)
    except: print 'Can not read file: '+File; sys.exit(1)
    ext = os.path.splitext(File)[1].lower()
    k = atr['FILE_TYPE']
    print 'File name   : '+os.path.basename(File)
    print 'File type   : '+atr['PROCESSOR']+' '+atr['FILE_TYPE']
    try:  atr['X_FIRST'];  print 'Coordinates : GEO'
    except:                print 'Coordinates : radar'


    #################### File Structure #####################
    try:
        argv[1]
        if argv[1] in ['--struct','--structure','--tree'] and ext in ['.h5','.he5']:
            print '***** HDF5 File Structure *****'
            print_hdf5_structure(File)
            return
    except: pass

    #################### HDF5 File Info #####################
    if ext in ['.h5','.he5']:
        h5file=h5py.File(File,'r')
        ##### Group Info
        print 'Al groups in this file:'
        print h5file.keys()

        ##### DateList / IgramList
        if k in ['interferograms','coherence','wrapped','timeseries']:
            epochList = h5file[k].keys()
            epochList = sorted(epochList)

    if k == 'timeseries':
        try: print_timseries_date_info(epochList)
        except: pass

        print '*************** Attributes **************'
        print_attributes(atr)

    elif k in ['interferograms', 'coherence', 'wrapped']:
        ##### Plot Attributes of One Epoch
        try: 
            epochNum = int(argv[1])
            epochAtr = h5file[k][epochList[epochNum-1]].attrs
            print '*****************************************'
            print epochList[epochNum-1]
            print '*************** Attributes **************'
            print_attributes(epochAtr)
            print '*****************************************'
            print epochList[epochNum-1]
        ##### Plot Epoch List Info
        except: 
            print '*****************************************'
            print 'Number of '+k+': '+str(len(epochList)) 
            print '*****************************************'
            print 'List of the '+k+':             number'
            for i in range(len(epochList)):
                print epochList[i] + '    ' + str(i+1)
            print '*****************************************'
            print 'Number of '+k+': '+str(len(epochList))

    ##### All other file types, except for timeseries/interferograms/coherence/wrapped
    else:
        print '*************** Attributes **************'
        print_attributes(atr)

    try: h5file.close()
    except: pass
    print '****************************************************************'
    return


############################################################
if __name__ == '__main__':
    main(sys.argv[1:])



