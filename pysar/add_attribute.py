#! /usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2016, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################
#


import os
import sys

import h5py
import numpy as np

import pysar.utils.readfile as readfile
import pysar.utils.writefile as writefile
import pysar.utils.utils as ut
import pysar.info as info



################################################################################

def usage():
    print('''
***************************************************************
  Add / Update attributes of hdf5 file manually

  Usage:
      add_attribute.py hdf5_file key1=value1 key2=value2 ...
      add_attribute.py hdf5_file attribute_list_file
      
      attribute_list_file is a text file like below, similar to PySAR template file
      track        =    422
      frame        =    620-650
      phase_ramp   =    quadratic
      
      ref_y=None - use None value to delete item from HDF5 file attributes

  Example:
      add_attribute.py timeseries.h5 add_attribute.txt
      add_attribute.py timeseries.h5 track=422 frame=650
      add_attribute.py unwrapIfgram.h5  ref_y=None  ref_x=None

***************************************************************
    ''')
    return


def main(argv):

    ##### Check Inputs
    if not argv or argv[0] in ['-h','--help']:
        usage();  sys.exit(1)
    if len(argv) < 2 or not argv[1]:
        raise Exception('\nAt lease 2 inputs are needed.\n')

    ##### Read Original Attributes
    #print '************ Add / Update HDF5 File Attributes *************'
    File = argv[0]
    atr  = readfile.read_attribute(File)
    print(('Input file is '+atr['PROCESSOR']+' '+atr['FILE_TYPE']+': '+File))

    ##### Read New Attributes
    atr_new = dict()
    for i in range(1,len(argv)):
        if os.path.isfile(argv[i]):
            atr_tmp = readfile.read_template(argv[i])
            atr_new.update(atr_tmp)
        else:
            atr_tmp = argv[i].split('=')
            atr_new[atr_tmp[0].strip()] = atr_tmp[1].strip()
    print("The following attributes will be added/updated, or removed if new value is 'None':")
    info.print_attributes(atr_new)

    ext = os.path.splitext(File)[1]
    ##### Update h5 File
    if ext in ['.h5','.he5']:
        File = ut.add_attribute(File, atr_new)
    else:
        if not ut.update_attribute_or_not(atr_new, atr):
            print('All updated (removed) attributes already exists (do not exists) and have the same value, skip update.')
        else:
            for key, value in iter(atr_new.items()):
                # delete the item is new value is None
                if value == 'None':
                    try: atr.pop(key)
                    except: pass
                else:
                    atr[key] = value
            if atr['PROCESSOR'] == 'roipac':
                print(('writing >>> '+File+'.rsc'))
                writefile.write_roipac_rsc(atr, File+'.rsc')

    return File


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])  

