#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2016, Yunjun Zhang                          #
# Author:  Yunjun Zhang                                    #
############################################################
#


import sys
import os

import h5py
import numpy as np

import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar.info as info



################################################################################

def usage():
    print '''
***************************************************************
  Add / Update attributes of hdf5 file manually

  Usage:
      add_attribute.py hdf5_file key1=value1 key2=value2 ...
      add_attribute.py hdf5_file attribute_list_file
      
      attribute_list_file is a text file like below, similar to PySAR template file
      track        =    422
      frame        =    620-650
      phase_ramp   =    quadratic

  Example:
      add_attribute.py timeseries.h5 add_attribute.txt
      add_attribute.py timeseries.h5 track=422 frame=650

***************************************************************
    '''
    return


def main(argv):

    ##### Check Inputs
    if not argv or argv[0] in ['-h','--help']:
        usage()
        sys.exit(1)
    if len(argv) < 2:  print('\nAt lease 2 inputs are needed.\n'); sys.exit(1)

    ##### Read Original Attributes
    print '************ Add / Update HDF5 File Attributes *************'
    File = argv[0]
    atr  = readfile.read_attributes(File)
    print 'Input file is '+atr['PROCESSOR']+' '+atr['FILE_TYPE']+': '+File

    ##### Read New Attributes
    atr_new = dict()
    for i in range(1,len(argv)):
        if os.path.isfile(argv[i]):
            atr_tmp = readfile.read_template(argv[i])
            atr_new.update(atr_tmp)
        else:
            atr_tmp = argv[i].split('=')
            atr_new[atr_tmp[0].strip()] = atr_tmp[1].strip()
    print 'The following attributes will be added/updated:'
    info.print_attributes(atr_new)

    ##### Update h5 File
    k = atr['FILE_TYPE']
    h5 = h5py.File(File,'r+')
    for key, value in atr_new.iteritems():
        h5[k].attrs[key] = value
    h5.close
    print 'Done.'

    return


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])  

