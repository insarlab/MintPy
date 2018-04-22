#! /usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2016, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################
#


import os, sys
import h5py
import numpy as np
from pysar.utils import readfile, writefile, utils as ut
import pysar.info as info


################################################################################
def usage():
    print('''
usage: add_attribute.py file metadataFile
       add_attribute.py file key1=value1 [key2=value2 [...]]

Add/Update attributes to file.

Example:
  add_attribute.py timeseries.h5 unavco_attribute.txt
  add_attribute.py timeseries.h5 track=422 frame=650

  Use None value to delete attribute:
  add_attribute.py unwrapIfgram.h5  ref_y=None  ref_x=None
    ''')
    return


def main(argv):
    ##### Check Inputs
    if not argv or argv[0] in ['-h','--help']:
        usage();  sys.exit(1)
    if len(argv) < 2 or not argv[1]:
        raise Exception('\nAt lease 2 inputs are needed.\n')

    ##### Read Original Attributes
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
    if ext in ['.h5','.he5']:
        File = ut.add_attribute(File, atr_new)
    else:
        if not ut.update_attribute_or_not(atr_new, atr):
            print('All updated (removed) attributes already exists (do not exists) and have the same value, skip update.')
        else:
            for key, value in iter(atr_new.items()):
                if value == 'None':
                    try: atr.pop(key)
                    except: pass
                else:
                    atr[key] = value
            print(('writing >>> '+File+'.rsc'))
            writefile.write_roipac_rsc(atr, out_file=File+'.rsc')
    return File


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])  

